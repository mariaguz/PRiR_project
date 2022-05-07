#include <omp.h>
#include <iostream>
#include <ctime>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>



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



	double* data1 = generate_data(1000000, 0, 1000);
	double* data2 = generate_data(2000000, 0, 2000);
	double* data3 = generate_data(3000000, 0 ,3000);


	float avg11 = computeAVG_indexes(1000000, data1, 10000, 10);
	float avg12 = computeAVG_values(1000000, data1, 9, 0);

	std::cout << "=====================\n";
	float avg21 = computeAVG_indexes(2000000, data2, 0, 400000);
	float avg22 = computeAVG_values(2000000, data2, 6, 9);

	std::cout << "=====================\n";
	float avg31 = computeAVG_indexes(3000000, data3, 0, 1700000);
	float avg32 = computeAVG_values(3000000, data3, 6, 9);








	return 0;
}
