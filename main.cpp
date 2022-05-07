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

//void zad_1()
//{
//    omp_init_lock(&lock);
//    srand(time(nullptr));
//    int a=0,i;
//    omp_set_num_threads(4);
//
////---------------------------------------------------------------------
//    printf("Petla firstprivate:\n");
//#pragma omp parallel for firstprivate(a)
//    for(i=0;i<10;i++)
//    {
//        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
//        a++;
//    }
//    printf("Po petli firstprivate a=%d\n\n",a);
////-------------------------------------------------------------------
//    printf("Petla private:\n");
//#pragma omp parallel for private(a)
//    for(a=0;a<10;a++)
//    {
//        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
//    }
//    printf("Po petli private a=%d\n\n",a);
////---------------------------------------------------------------------
//    printf("Petla lastprivate:\n");
//#pragma omp parallel for lastprivate(a)
//    for(i=0;i<10;i++)
//    {
//        //1. spróbuj zmienić tą wartość na zmienną losową i zobacz jak to działa
//        a = omp_get_thread_num();
//
//        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
//    }
//    printf("Po petli lastprivate a=%d\n\n",a);
////---------------------------------------------------------------------
//    printf("Petla shared:\n");
//    a=0;
//#pragma omp parallel for shared(a)
//    for(i=0;i<10;i++)
//    {
//        //2. Co się stanie gdy wyłączymy zamek?
//        omp_set_lock(&lock);
//        a=omp_get_thread_num();
//        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
//        omp_unset_lock(&lock);
//    }
//
//    //3. Jaka bedzie wartosc "a" po kilkukrotnym wywołaniu programu?
//    printf("Po petli shared a=%d\n\n",a);
////---------------------------------------------------------------------
//    printf("Petla bez zadnej klauzuli:\n");
//    a=0;
//#pragma omp parallel for
//    for(i=0;i<10;i++)
//    {
//        a++;
//        printf("Thread %d Iteration %d a=%d\n",omp_get_thread_num(),i,a);
//    }
//    //4. Jaka jest domysla klauzula?
//    printf("Po petli bez klauzuli a=%d\n",a);
//}
//
//double zad_2(int n, int m, int p){
//
//    const size_t N = 100;
//    const size_t M = 500;
//    const size_t P = 140;
//
//    int m_1[N][M];
//    int m_2[M][P];
//    int m_3[N][P];
//
//
//    for (auto & i : m_1){
//        for (int & j : i){
//            j = 2;
//        }
//    }
//
//    for (auto & i : m_2){
//        for (int & j : i){
//            j = 3;
//        }
//    }
//
//
//    double t_0 = omp_get_wtime();
//    for(int i=0;i<N;i++){
//        for(int j=0;j<P;j++){
//            for(int k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    double t_1 = omp_get_wtime();
//    double t = t_1 - t_0;
//    printf("Bez zrownoleglenia:\n czas: %f\n\n", t);
//
//    int m_check[N][P];
//
//    for(int i=0;i<N;i++){
//        for(int j=0;j<P;j++){
////            std::cout << m_3[i][j];
//            m_check[i][j] = m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    std:: cout << "\n\n";
//
//
//
//
//
//
//
////    ================================ zrownoleglenie petli zewnetrznej
//// -----------------shared static
//
//
//
//
//    t_0 = omp_get_wtime();
//    #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static)
//    for(int i=0;i<N;i++){
//        for(int j=0;j<P;j++){
//            for(int k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla zewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static):\n czas: %f\n\n", t);
//
////    for(auto & i : m_3){
////        for(int j : i){
////            std::cout << j;
////        }
////        std:: cout << "\n";
////    }
////
////    for(int i=0;i<N;i++){
////        for(int j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//
//    std:: cout << "\n\n";
////    ================================ zrownoleglenie petli zewnetrznej
//// -----------------shared dynamic
//
//
//
//
//    t_0 = omp_get_wtime();
//    #pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic)
//    for(int i=0;i<N;i++){
//        for(int j=0;j<P;j++){
//            for(int k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla zewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic):\n czas: %f\n\n", t);
//
//
////    for(auto & i : m_3){
////        for(int j : i){
////            std::cout << j;
////        }
////        std:: cout << "\n";
////    }
////
////    for(int i=0;i<N;i++){
////        for(int j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//
//    std:: cout << "\n\n";
//
////    ================================ zrownoleglenie petli zewnetrznej
//// -----------------private dynamic
//int i,j,k;
//
//
//
//    t_0 = omp_get_wtime();
//#pragma omp parallel for private(i,j,k) schedule(static)
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            for( k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla zewnetrzna \n #pragma omp parallel for private(i,j,k) schedule(static):\n czas: %f\n\n", t);
//
//
////    for(i=0;i<N;i++){
////        for(j=0;j<P;j++){
////            std::cout << m_3[i][j];
////        }
////        std:: cout << "\n";
////    }
////
////    for( i=0;i<N;i++){
////        for( j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//    std:: cout << "\n\n";
//
////    ================================ zrownoleglenie petli zewnetrznej
//// -----------------private dynamic
//
//
//
//
//    t_0 = omp_get_wtime();
//#pragma omp parallel for private(i,j,k) schedule(dynamic)
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            for( k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla zewnetrzna \n #pragma omp parallel for private(i,j,k) schedule(dynamic):\n czas: %f\n\n", t);
//
//
////    for(i=0;i<N;i++){
////        for(j=0;j<P;j++){
////            std::cout << m_3[i][j];
////        }
////        std:: cout << "\n";
////    }
////
////    for( i=0;i<N;i++){
////        for( j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//
//    std:: cout << "\n\n";
//
//////    ================================ zrownoleglenie petli wewnetrznej
////
////
////
////// -----------------shared static
//
//    t_0 = omp_get_wtime();
//
//    for( i=0;i<N;i++){
//        #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static)
//        for( j=0;j<P;j++){
//            for( k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla wewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static):\n czas: %f\n\n", t);
//
////    for(i=0;i<N;i++){
////        for(j=0;j<P;j++){
////            std::cout << m_3[i][j];
////        }
////        std:: cout << "\n";
////    }
////
////    for( i=0;i<N;i++){
////        for( j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//
//    std:: cout << "\n\n";
//
//
//
//    //    ================================ zrownoleglenie petli wewnetrznej
//
//
//
//// -----------------shared dynamic
//
//    t_0 = omp_get_wtime();
//
//    for( i=0;i<N;i++){
//#pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic)
//        for( j=0;j<P;j++){
//            for( k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla wewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic):\n czas: %f\n\n", t);
//
////    for(i=0;i<N;i++){
////        for(j=0;j<P;j++){
////            std::cout << m_3[i][j];
////        }
////        std:: cout << "\n";
////    }
////
////    for( i=0;i<N;i++){
////        for( j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//
//    std:: cout << "\n\n";
//
//
//    ////    ================================ zrownoleglenie petli wewnetrznej
////
////
////
////// -----------------private static
//
//    t_0 = omp_get_wtime();
//
//    for( i=0;i<N;i++){
//#pragma omp parallel for private(j,k) schedule(static)
//        for( j=0;j<P;j++){
//            for( k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla wewnetrzna \n #pragma omp parallel for private(j,k) schedule(static):\n czas: %f\n\n", t);
//
////    for(i=0;i<N;i++){
////        for(j=0;j<P;j++){
////            std::cout << m_3[i][j];
////        }
////        std:: cout << "\n";
////    }
////
////    for( i=0;i<N;i++){
////        for( j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//    std:: cout << "\n\n";
//
//
//    ////    ================================ zrownoleglenie petli wewnetrznej
////
////
////
////// -----------------private dynamic
//
//    t_0 = omp_get_wtime();
//
//    for( i=0;i<N;i++){
//#pragma omp parallel for private(j,k) schedule(dynamic)
//        for( j=0;j<P;j++){
//            for( k=0; k<M; k++){
//                m_3[i][j] = m_1[i][k] * m_2[k][j];
//            }
//        }
//    }
//    t_1 = omp_get_wtime();
//    t = t_1 - t_0;
//    printf("Petla wewnetrzna \n #pragma omp parallel for private(j,k) schedule(dynamic):\n czas: %f\n\n", t);
//
////    for(i=0;i<N;i++){
////        for(j=0;j<P;j++){
////            std::cout << m_3[i][j];
////        }
////        std:: cout << "\n";
////    }
////
////    for( i=0;i<N;i++){
////        for( j=0;j<P;j++){
////            if(m_check[i][j] == m_3[i][j]){
////                printf(" | ok | ");
////            }
////        }
////        std:: cout << "\n";
////    }
//
//    std:: cout << "\n\n";
//    return 0;
//////    return ts;
//
//}

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

double* computeAVG_values(int n, double* data, int val_min, int val_max) {
	auto beginTime = std::chrono::high_resolution_clock::now();
	
	double avg = 0, sum = 0;
	int counter = 0;

	float min = *std::min_element(data, data + n);
	float max = *std::max_element(data, data + n);
	std::cout << "\nMin = " << min << "\nMax = " << max;

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
			std::cout<<"\nsum: " << sum;
			std::cout << "\ncount: " << counter;
			avg = sum / counter;

			std::cout << "\nSrednia (zakres wartosci) = " << avg;
			auto endTime = std::chrono::high_resolution_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);

			printf("\nCzas: %f\n", time.count() * 1e-3);
			
			double results[4] = { avg, sum, counter, time.count() * 1e-3 };
			
			return results;
		}
	}

}

int main() {


	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = sequence = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	double* data = generate_data(SIZE, 0, 1000);


	//double avg1 = computeAVG_indexes(SIZE, data, 10, 100)[0];
	//double avg2 = computeAVG_values(SIZE, data, 50, 100)[0];
	
	double* results = computeAVG_values(SIZE, data, 0, 10);

	std::cout << results[0] << " " << results[1] << " " << results[2];

	
	
	////= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = OpenMP = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
//	double* avgs_OpenMP = new double [THREADS];
//	
//	
//	
//	auto beginTime = std::chrono::high_resolution_clock::now();
//
//	omp_init_lock(&lock);
//	omp_set_num_threads(THREADS);
//
//#pragma omp parallel
//	{
//		//pobranie ID wątku
//		int rank = omp_get_thread_num();
//		//printf("thread %d", omp_get_thread_num());
//
//		// obliczenie przedziałów dla których srednie beda liczone rownolegle
//		
//		int begin = round((SIZE / THREADS) * rank);
//		int end = round(((SIZE / THREADS) * rank) + SIZE/THREADS);
//		//printf("thread: %d, begin: %d, end: %d \n",rank,begin, end);
//
//		for (int i = begin; i < end; i++) {
//
//		}
//		
//
//	}






	return 0;
}
