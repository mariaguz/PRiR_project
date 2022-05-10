#include <omp.h>
#include "mpi.h" 
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <locale.h>

#define THREADS 7
#define SIZE 1e6
#define RANGE 10000
#define LB 2000
#define UB 5000



double* generate_data(int n, int low, int high) {

	double* data = new double[n];

	for (int i = 0; i < n; i++) {
		data[i] = ((double)(rand() % 10000) / 10000.0) * (high - low) + low;
	}

	return data;
}

int main(int argc, char** argv) {

    setlocale (LC_ALL,"");
	int rank, size, tag, i;
    double* data;
	tag = 0;
	rank = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//
    if (rank == 0) {
        std::cout << "Watki " << THREADS << "\tProcesy " << size;
        ////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = sequence = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        data = generate_data(SIZE, 0, RANGE);

        double avg_seq, sum_seq = 0;
        int count_seq = 0;

        auto beginTime = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < SIZE; i++) {
            if (data[i] >= LB && data[i] <= UB) {
                sum_seq += data[i];
                count_seq++;
            }
        }
        avg_seq = sum_seq / count_seq;
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - beginTime);
        //Wypisanie rezultatów
        printf("\nCzas Sekwencyjne: %f", time.count() * 1e-6);
        printf("\nAVG Sekwencyjne: %.2f", avg_seq);


		////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = OpenMP = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		int array_size = SIZE;
        double sum_OpenMP = 0, avg_OpenMP;
        int count_OpenMP = 0;
        auto beginTime_OMP = std::chrono::high_resolution_clock::now();

#pragma omp parallel for num_threads(THREADS) reduction(+:sum_OpenMP) reduction(+:count_OpenMP)
        for (i = 0; i < array_size; ++i) {
            if (data[i] >= LB && data[i] <= UB) {
                sum_OpenMP += data[i];
                count_OpenMP++;
            }
        }
        avg_OpenMP = sum_OpenMP / count_OpenMP;
        time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - beginTime_OMP);
		//Wypisanie rezultatów
		printf("\nCzas OpenMP: %f", time.count() * 1e-6);
		printf("\nAVG OpenMP: %.2f", avg_OpenMP);
	}
	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = MPI = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

	MPI_Barrier(MPI_COMM_WORLD);
	int count = SIZE/(size), rest = SIZE - count * (size);
	double sum_MPI = 0, cnt_MPI = 0, avg_MPI;
	double * local_array;
    double sum = 0;
    double count_in_tread = 0;
    double* results = new double[2];
    double end;
    int cnt;

//    double start = MPI_Wtime();
    auto s = std::chrono::high_resolution_clock::now();

	if (rank == 0) {
        for (i = 0; i < count; ++i) {
            if (data[i] >= LB && data[i] <= UB) {
                sum += data[i];
                count_in_tread++;
            }
        }
        for(int dest = 1; dest < size; ++dest)
        {
            if(dest < size - 1) {
                cnt = count;
            } else {
                cnt = count + rest;
            }

            MPI_Send(&data[(dest) * cnt], cnt, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            //printf("\nP0 sent a %d elements to P%d.", cnt, dest);
        }
	}

	if (rank != 0) {
	    if(rank == size - 1){
            cnt = count + rest;
	        local_array = new double [cnt];
	    } else {
	        cnt = count;
	        local_array = new double [cnt];
	    }
        MPI_Recv(local_array, cnt, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (i = 0; i < cnt; ++i) {
            if (local_array[i] >= LB && local_array[i] <= UB) {
                sum_MPI += local_array[i];
                cnt_MPI++;
            }
        }
        results[0] = sum_MPI;
        results[1] = cnt_MPI;
        MPI_Send(&results[0], 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

	}
	MPI_Barrier(MPI_COMM_WORLD);


	if (rank == 0) {
        for(int src = 1; src < size; ++src)
        {
            MPI_Recv(results, 2, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += results[0];
            count_in_tread += results[1];
        }

        avg_MPI = sum / count_in_tread;
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - s);
//        end = MPI_Wtime();

        printf("\nCzas MPI: %f", time.count() * 1e-6);
        printf("\nAVG MPI: %.2f", avg_MPI);
    }

    ////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = HYBRID = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    MPI_Barrier(MPI_COMM_WORLD);

	double sum_hybrid = 0, cnt_hybrid = 0;

//    start = MPI_Wtime();
    s = std::chrono::high_resolution_clock::now();
    if (rank == 0) {
#pragma omp parallel for num_threads(THREADS) reduction(+:sum) reduction(+:count_in_tread)
        for (i = 0; i < count; ++i) {
            if (data[i] >= LB && data[i] <= UB) {
                sum += data[i];
                count_in_tread++;
            }
        }
        for(int dest = 1; dest < size; ++dest)
        {
            if(dest < size - 1) {
                cnt = count;
            } else {
                cnt = count + rest;
            }

            MPI_Send(&data[(dest) * cnt], cnt, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            //printf("\nP0 sent a %d elements to P%d.", cnt, dest);
        }
    }

    if (rank != 0) {
        if(rank == size - 1){
            cnt = count + rest;
            local_array = new double [cnt];
        } else {
            cnt = count;
            local_array = new double [cnt];
        }
        MPI_Recv(local_array, cnt, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for num_threads(THREADS) reduction(+:sum_hybrid) reduction(+:cnt_hybrid)
        for (i = 0; i < cnt; ++i) {
            if (local_array[i] >= LB && local_array[i] <= UB) {
                sum_hybrid += local_array[i];
                cnt_hybrid++;
            }
        }

        results[0] = sum_hybrid;
        results[1] = cnt_hybrid;
        MPI_Send(&results[0], 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

    }
    MPI_Barrier(MPI_COMM_WORLD);


    if (rank == 0) {
        for(int src = 1; src < size; ++src)
        {

            MPI_Recv(results, 2, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += results[0];
            count_in_tread += results[1];
        }

        avg_MPI = sum / count_in_tread;
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - s);
//        end = MPI_Wtime();

        printf("\nCzas HYBRID: %f", time.count() * 1e-6);
        printf("\nAVG HYBRID: %.2f", avg_MPI);
    }

    MPI_Finalize();
	return 0;
}
