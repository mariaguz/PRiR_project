#include <omp.h>
#include "mpi.h" 
#include <iostream>
#include <ctime>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>
#include <locale.h>

omp_lock_t lock;

#define THREADS 4
#define SIZE 3e6
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

double* split_data(double* data, int begin, int end) {

	int range = end - begin;

	if (begin > end) {
		std::cout << "\n(Split_data) ERROR: Begin>end.\n";
		return NULL;
	}
	else {
		double* splitted_data = new double[range];
		int j = 0;

		for (int i = begin; i < end; i++) {
			splitted_data[j] = data[i];
			j++;
		}

		return splitted_data;
	}
	
}

/*
double computeAVG_indexes(int n, double* data, int first, int last) {
	auto beginTime = std::chrono::high_resolution_clock::now();

	float avg = 0;
	int counter = 0;

	if (first < 0 || last >= n) {
		std::cout << "\n(Srednia z zakresu indeksow) ERROR: Indeksy poza zakresem.\n";
		return -1;
	}
	else {
		if (first > last) {
			std::cout << "\n(Srednia z zakresu indeksow) ERROR: Indeks 'last' jest mniejszy od 'first'\n";
			return -1;
		}
		else {
			for (int i = first; i <= last; i++) {
				avg += data[i];
				counter++;
			}

			avg = avg / counter;

			std::cout << "\nSrednia (zakres indeksow) = " << avg;
			auto endTime = std::chrono::high_resolution_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);

			printf("\nCzas: %f\n", time.count() * 1e-3);
			return avg;
		}
	}

	


}
*/

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
//			std::cout << "\nSrednia (zakres wartosci) = " << avg;

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

    setlocale (LC_ALL,"");
    omp_init_lock(&lock);
	int rank, size, tag;
    double* data;
	tag = 0;
	rank = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (rank == 0) {
        std::cout << "Watki " << THREADS << "\tProcesy " << size;
		////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = sequence = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		data = generate_data(SIZE, 0, RANGE);

		double* results_seq = computeAVG_values(SIZE, data, LB, UB);

		//Wypisanie rezultatów
		printf("\nCzas Sekwencyjne: %f", results_seq[3]);
		printf("\nAVG Sekwencyjne: %.2f", results_seq[0]);





		////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = OpenMP = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
		//double* data = generate_data(SIZE, 0, 1000);
		double* sums_OpenMP = new double[THREADS];
		double* counts_OpenMP = new double[THREADS];



		auto beginTime = std::chrono::high_resolution_clock::now();

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
			double* results = computeAVG_values(subsize, splitted_data, LB, UB);

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
		printf("\nAVG OpenMP: %.2f", avg_OpenMP);
	}
	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = MPI = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	int count = SIZE/(size-1);
	int rest = SIZE - count * (size-1);
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
            cnt = count + rest;
	        local_array = new double [cnt];
	    } else {
	        cnt = count;
	        local_array = new double [cnt];
	    }
        MPI_Recv(local_array, cnt, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    double* results = computeAVG_values(cnt, local_array, LB, UB);
        MPI_Send(&results[0], 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

	}
	MPI_Barrier(MPI_COMM_WORLD);


	if (rank == 0) {
	    double sum = 0;
	    double count_in_tread = 0;
        for(int src = 1; src < size; ++src)
        {
            double* results = new double[4];
            MPI_Recv(results, 4, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += results[1];
            count_in_tread += results[2];
        }

        double avg_MPI = sum / count_in_tread;
        double end = MPI_Wtime();
        printf("\nCzas MPI: %f", end - start);
        printf("\nAVG MPI: %.2f", avg_MPI);
    }

    ////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = HYBRID = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
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
            cnt = count + rest;
            local_array = new double [cnt];
        } else {
            cnt = count;
            local_array = new double [cnt];
        }
        MPI_Recv(local_array, cnt, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        omp_set_num_threads(THREADS);
        double* sums_hyb = new double[THREADS];
        double* counts_hyb = new double[THREADS];
#pragma omp parallel
        {
            //pobranie ID wątku
            int tread_num = omp_get_thread_num();

            // obliczenie przedziałów dla których srednie beda liczone rownolegle
            int begin = round(double(cnt / THREADS) * tread_num);
            int end = round(double(cnt / THREADS) * tread_num + cnt / THREADS);
            int subsize = end - begin;

            //podział zbioru
            double* splitted_data = split_data(local_array, begin, end);

            //obliczenie wyników w każdym wątku
            double* results = computeAVG_values(subsize, splitted_data, LB, UB);

            sums_hyb[tread_num] = results[1];
            counts_hyb[tread_num] = results[2];

        }

        //obliczenie średniej ogólnej
        double results[2] = {0};
        for (int i = 0; i < THREADS; i++) {
            results[0] += sums_hyb[i];
            results[1] += counts_hyb[i];
        }

        MPI_Send(&results[0], 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

    }
    MPI_Barrier(MPI_COMM_WORLD)

    if (rank == 0) {
        double sum = 0;
        double count_in_tread = 0;
        for(int src = 1; src < size; ++src)
        {
            double* results = new double[4];
            MPI_Recv(results, 2, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += results[0];
            count_in_tread += results[1];
        }

        double avg_hyb = sum / count_in_tread;
        double end = MPI_Wtime();
        printf("\nCzas HYBRID: %f", end - start);
        printf("\nAVG HYBRID: %.2f", avg_hyb);
    }

    MPI_Finalize();
	return 0;
}
