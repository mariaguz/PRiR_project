#include <omp.h>
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

int main() {


	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = sequence = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	double* data = generate_data(SIZE, 0, RANGE);


	double* results_seq = computeAVG_values(SIZE, data, 1000, 9000);

	//Wypisanie rezultatów
	printf("\nCzas Sekwencyjne: %f", results_seq[3]);
	printf("\nAVG Sekwencyjne: %f", results_seq[0]);
	
	
	
	
	
	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = OpenMP = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
	//double* data = generate_data(SIZE, 0, 1000);
	double* sums_OpenMP = new double [THREADS];
	double* counts_OpenMP = new double[THREADS];


	
	auto beginTime = std::chrono::high_resolution_clock::now();

	omp_init_lock(&lock);
	omp_set_num_threads(THREADS);

#pragma omp parallel
	{
		//pobranie ID wątku
		int rank = omp_get_thread_num();
		double size = SIZE;
		double threads = THREADS;



		// obliczenie przedziałów dla których srednie beda liczone rownolegle
		int begin = round(double(size / threads) * rank);
		int end = round(double(size / threads) * rank + size / threads);
		int subsize = end - begin;

		//podział zbioru
		double* splitted_data = new double[subsize];
		splitted_data = split_data(data, begin, end);

		//obliczenie wyników w każdym wątku
		double* results = computeAVG_values(subsize, splitted_data, 1000, 9000);
		
		sums_OpenMP[rank] = results[1];
		counts_OpenMP[rank] = results[2];

	}

	//obliczenie średniej ogólnej
	double sum_OpenMP=0;
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



	return 0;
}
