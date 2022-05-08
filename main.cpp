#include <omp.h>
#include "mpi.h" 
#include <iostream>
#include <ctime>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>

omp_lock_t lock;

#define THREADS 4
#define SIZE 1000000
#define RANGE 10000


double* computeAVG_values(int n, double* data, int val_min, int val_max) {
	auto beginTime = std::chrono::high_resolution_clock::now();
	
	double avg = 0, sum = 0;
	int counter = 0;

	float min = *std::min_element(data, data + n);
	float max = *std::max_element(data, data + n);
	//std::cout << "\nMin = " << min << "\nMax = " << max;

	if (val_min < min || val_max > max) {
		std::cout << "\n(Srednia z zakresu wartosci) ERROR: Wartosci poza zakresem.\n";
		return NULL;
	}
	else {
		if (val_min > val_max) {
			std::cout << "\n(Srednia z zakresu wartosci) ERROR: Wartosc min > max.\n";
			return NULL;
		}
		else {
			for (int i = 0; i < n; i++) {
				if (data[i] >= val_min && data[i] <= val_max) {


					sum += data[i];
					counter++;
				}
			}

			avg = sum / counter;
			//std::cout << "\nSrednia (zakres wartosci) = " << avg;

			auto endTime = std::chrono::high_resolution_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);

			double* results = new double[4];

			results[0] = avg;
			results[1] = sum;
			results[2] = counter;
			results[3] = time.count() * 1e-3;

			return results;
		}
	}

}

double* generate_data(int n, int low, int high) {
	double* data = new double[n];
	for (int i = 0; i < n; i++) {
		data[i] = ((double)(rand() % 10000) / 10000.0) * (high - low) + low;
	}
	return data;
}

double* split_data(double* data, int begin, int end) {
	int size = end - begin;
	double* subset = new double[size];
	for (int i = 0; i < size; i++)
	{
		subset[i] = data[i + begin];
	}
	return subset;
}

int main(int argc, char** argv) {
	int rank, ranksent, size, source, tag, i, len, src, number = 5;
    double* data;
	std::string num;
	len = 1000;
	tag = 0;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = sequence = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		data = generate_data(SIZE, 0, RANGE);


		double* results_seq = computeAVG_values(SIZE, data, 1000, 9000);

		//Wypisanie rezultatów
		printf("\nCzas Sekwencyjne: %f", results_seq[3]);
		printf("\nAVG Sekwencyjne: %f", results_seq[0]);





		////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = OpenMP = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
		//double* data = generate_data(SIZE, 0, 1000);
		double* sums_OpenMP = new double[THREADS];
		double* counts_OpenMP = new double[THREADS];



		auto beginTime = std::chrono::high_resolution_clock::now();

		omp_init_lock(&lock);
		omp_set_num_threads(THREADS);

#pragma omp parallel
		{
			//pobranie ID wątku
			int tread_num = omp_get_thread_num();

			// obliczenie przedziałów dla których srednie beda liczone rownolegle
			int begin = round(double(SIZE / THREADS) * tread_num);
			int end = round(double(SIZE / THREADS) * tread_num + SIZE / THREADS);
			int subsize = end - begin;

			//podział zbioru
			double* splitted_data = new double[subsize];
			splitted_data = split_data(data, begin, end);

			//obliczenie wyników w każdym wątku
			double* results = computeAVG_values(subsize, splitted_data, 1000, 9000);

			sums_OpenMP[tread_num] = results[1];
			counts_OpenMP[tread_num] = results[2];

		}

		//obliczenie średniej ogólnej
		double sum_OpenMP = 0;
		double count_OpenMP = 0;

		for (int i = 0; i < THREADS; i++) {
			sum_OpenMP = sum_OpenMP + sums_OpenMP[i];
			count_OpenMP = count_OpenMP + counts_OpenMP[i];
		}

		double avg_OpenMP = sum_OpenMP / count_OpenMP;


		auto endTime = std::chrono::high_resolution_clock::now();
		auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);

		//Wypisanie rezultatów
		printf("\nCzas OpenMP: %f", time.count() * 1e-3);
		printf("\nAVG OpenMP: %f", avg_OpenMP);
	}
	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = MPI = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	int* array = new int[len];
	int count = SIZE/size;
	int rest = SIZE - count * size;
	double * local_array;
	if (rank == 0) {
	    int cnt;
        for(int dest = 1; dest < size; ++dest)
        {
            if(dest < size - 1) {
                cnt = count;
            } else {
                cnt = count + rest;
            }

            MPI_Send(&data[(dest - 1) * cnt], cnt, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            //printf("\nP0 sent a %d elements to P%d.", cnt, dest);
        }
	}

	if (rank != 0) {
	    int cnt;
	    if(rank == size - 1){
	        local_array = new double [count + rest];
	        cnt = count + rest;
	    } else {
	        cnt = count;
	        local_array = new double [count];
	    }
        MPI_Recv(local_array, cnt, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    double* results = computeAVG_values(cnt, local_array, 1000, 9000);
        MPI_Send(&results[0], 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

	}
	MPI_Barrier(MPI_COMM_WORLD);


	if (rank == 0) {
	    double sum = 0;
	    double count = 0;
        for(int src = 1; src < size; ++src)
        {
            double* results = new double[4];
            MPI_Recv(results, 4, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += results[1];
            count += results[2];
        }

        double avg_MPI = sum/count;
        double end = MPI_Wtime();
        printf("\nCzas MPI: %f", end - start);
        printf("\nAVG MPI: %f", avg_MPI);
    }
    MPI_Finalize();
	return 0;
}
