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



double* generate_data(int n, int low, int high) {

	double* data = new double[n];

	for (int i = 0; i < n; i++) {
		data[i] = ((double)(rand() % 10000) / 10000.0) * (high - low) + low;
	}

	return data;
}

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

double computeAVG_values(int n, double* data, int val_min, int val_max) {
	auto beginTime = std::chrono::high_resolution_clock::now();
	
	float avg = 0;
	int counter = 0;

	float min = *std::min_element(data, data + n);
	float max = *std::max_element(data, data + n);
	//std::cout << "\nMin = " << min << "\nMax = " << max;

	if (val_min < min || val_max > max) {
		std::cout << "\n(Srednia z zakresu wartosci) ERROR: Wartosci poza zakresem.\n";
		return -1;
	}
	else {
		if (val_min > val_max) {
			std::cout << "\n(Srednia z zakresu wartosci) ERROR: Wartosc min > max.\n";
			return -1;
		}
		else {
			for (int i = 0; i < n; i++) {
				if (data[i] >= min && data[i] <= max) {


					avg += data[i];
					counter++;
				}
			}

			avg = avg / counter;

			std::cout << "\nSrednia (zakres wartosci) = " << avg;
			auto endTime = std::chrono::high_resolution_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);

			printf("\nCzas: %f\n", time.count() * 1e-3);
			
			return avg;
		}
	}

}

int main() {


	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = sequence = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	double* data = generate_data(SIZE, 0, 1000);


	//float avg1 = computeAVG_indexes(SIZE, data, 10000, 10);
	//float avg2 = computeAVG_values(SIZE, data, 9, 0);


	
	
	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = OpenMP = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
	double* avgs_OpenMP = new double [THREADS];
	
	
	
	auto beginTime = std::chrono::high_resolution_clock::now();

	omp_init_lock(&lock);
	omp_set_num_threads(THREADS);

#pragma omp parallel
	{
		//pobranie ID wątku
		int rank = omp_get_thread_num();
		//printf("thread %d", omp_get_thread_num());

		// obliczenie przedziałów dla których srednie beda liczone rownolegle
		
		int begin = round((SIZE / THREADS) * rank);
		int end = round(((SIZE / THREADS) * rank) + SIZE/THREADS);
		//printf("thread: %d, begin: %d, end: %d \n",rank,begin, end);

		for(int i = begin; i <end; i++)
		

	}






	return 0;
}
