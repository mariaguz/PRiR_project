#include <omp.h>
#include <iostream>
#include <ctime>
#include <stdexcept>
#include <cstdlib>
#include <string>



omp_lock_t lock;

void zad_1()
{
    omp_init_lock(&lock);
    srand(time(nullptr));
    int a=0,i;
    omp_set_num_threads(4);

//---------------------------------------------------------------------
    printf("Petla firstprivate:\n");
#pragma omp parallel for firstprivate(a)
    for(i=0;i<10;i++)
    {
        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
        a++;
    }
    printf("Po petli firstprivate a=%d\n\n",a);
//-------------------------------------------------------------------
    printf("Petla private:\n");
#pragma omp parallel for private(a)
    for(a=0;a<10;a++)
    {
        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
    }
    printf("Po petli private a=%d\n\n",a);
//---------------------------------------------------------------------
    printf("Petla lastprivate:\n");
#pragma omp parallel for lastprivate(a)
    for(i=0;i<10;i++)
    {
        //1. spróbuj zmienić tą wartość na zmienną losową i zobacz jak to działa
        a = omp_get_thread_num();

        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
    }
    printf("Po petli lastprivate a=%d\n\n",a);
//---------------------------------------------------------------------
    printf("Petla shared:\n");
    a=0;
#pragma omp parallel for shared(a)
    for(i=0;i<10;i++)
    {
        //2. Co się stanie gdy wyłączymy zamek?
        omp_set_lock(&lock);
        a=omp_get_thread_num();
        printf("Watek %d a=%d\n",omp_get_thread_num(),a);
        omp_unset_lock(&lock);
    }

    //3. Jaka bedzie wartosc "a" po kilkukrotnym wywołaniu programu?
    printf("Po petli shared a=%d\n\n",a);
//---------------------------------------------------------------------
    printf("Petla bez zadnej klauzuli:\n");
    a=0;
#pragma omp parallel for
    for(i=0;i<10;i++)
    {
        a++;
        printf("Thread %d Iteration %d a=%d\n",omp_get_thread_num(),i,a);
    }
    //4. Jaka jest domysla klauzula?
    printf("Po petli bez klauzuli a=%d\n",a);
}

double zad_2(int n, int m, int p){

    const size_t N = 100;
    const size_t M = 500;
    const size_t P = 140;

    int m_1[N][M];
    int m_2[M][P];
    int m_3[N][P];


    for (auto & i : m_1){
        for (int & j : i){
            j = 2;
        }
    }

    for (auto & i : m_2){
        for (int & j : i){
            j = 3;
        }
    }


    double t_0 = omp_get_wtime();
    for(int i=0;i<N;i++){
        for(int j=0;j<P;j++){
            for(int k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    double t_1 = omp_get_wtime();
    double t = t_1 - t_0;
    printf("Bez zrownoleglenia:\n czas: %f\n\n", t);

    int m_check[N][P];

    for(int i=0;i<N;i++){
        for(int j=0;j<P;j++){
//            std::cout << m_3[i][j];
            m_check[i][j] = m_3[i][j];
        }
        std:: cout << "\n";
    }

    std:: cout << "\n\n";







//    ================================ zrownoleglenie petli zewnetrznej
// -----------------shared static




    t_0 = omp_get_wtime();
    #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static)
    for(int i=0;i<N;i++){
        for(int j=0;j<P;j++){
            for(int k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla zewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static):\n czas: %f\n\n", t);

//    for(auto & i : m_3){
//        for(int j : i){
//            std::cout << j;
//        }
//        std:: cout << "\n";
//    }
//
//    for(int i=0;i<N;i++){
//        for(int j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }

    std:: cout << "\n\n";
//    ================================ zrownoleglenie petli zewnetrznej
// -----------------shared dynamic




    t_0 = omp_get_wtime();
    #pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic)
    for(int i=0;i<N;i++){
        for(int j=0;j<P;j++){
            for(int k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla zewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic):\n czas: %f\n\n", t);


//    for(auto & i : m_3){
//        for(int j : i){
//            std::cout << j;
//        }
//        std:: cout << "\n";
//    }
//
//    for(int i=0;i<N;i++){
//        for(int j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }

    std:: cout << "\n\n";

//    ================================ zrownoleglenie petli zewnetrznej
// -----------------private dynamic
int i,j,k;



    t_0 = omp_get_wtime();
#pragma omp parallel for private(i,j,k) schedule(static)
    for( i=0;i<N;i++){
        for( j=0;j<P;j++){
            for( k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla zewnetrzna \n #pragma omp parallel for private(i,j,k) schedule(static):\n czas: %f\n\n", t);


//    for(i=0;i<N;i++){
//        for(j=0;j<P;j++){
//            std::cout << m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }
    std:: cout << "\n\n";

//    ================================ zrownoleglenie petli zewnetrznej
// -----------------private dynamic




    t_0 = omp_get_wtime();
#pragma omp parallel for private(i,j,k) schedule(dynamic)
    for( i=0;i<N;i++){
        for( j=0;j<P;j++){
            for( k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla zewnetrzna \n #pragma omp parallel for private(i,j,k) schedule(dynamic):\n czas: %f\n\n", t);


//    for(i=0;i<N;i++){
//        for(j=0;j<P;j++){
//            std::cout << m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }

    std:: cout << "\n\n";

////    ================================ zrownoleglenie petli wewnetrznej
//
//
//
//// -----------------shared static

    t_0 = omp_get_wtime();

    for( i=0;i<N;i++){
        #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static)
        for( j=0;j<P;j++){
            for( k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla wewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(static):\n czas: %f\n\n", t);

//    for(i=0;i<N;i++){
//        for(j=0;j<P;j++){
//            std::cout << m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }

    std:: cout << "\n\n";



    //    ================================ zrownoleglenie petli wewnetrznej



// -----------------shared dynamic

    t_0 = omp_get_wtime();

    for( i=0;i<N;i++){
#pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic)
        for( j=0;j<P;j++){
            for( k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla wewnetrzna \n #pragma omp parallel for shared(m_1, m_2, m_3) schedule(dynamic):\n czas: %f\n\n", t);

//    for(i=0;i<N;i++){
//        for(j=0;j<P;j++){
//            std::cout << m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }

    std:: cout << "\n\n";


    ////    ================================ zrownoleglenie petli wewnetrznej
//
//
//
//// -----------------private static

    t_0 = omp_get_wtime();

    for( i=0;i<N;i++){
#pragma omp parallel for private(j,k) schedule(static)
        for( j=0;j<P;j++){
            for( k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla wewnetrzna \n #pragma omp parallel for private(j,k) schedule(static):\n czas: %f\n\n", t);

//    for(i=0;i<N;i++){
//        for(j=0;j<P;j++){
//            std::cout << m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }
    std:: cout << "\n\n";


    ////    ================================ zrownoleglenie petli wewnetrznej
//
//
//
//// -----------------private dynamic

    t_0 = omp_get_wtime();

    for( i=0;i<N;i++){
#pragma omp parallel for private(j,k) schedule(dynamic)
        for( j=0;j<P;j++){
            for( k=0; k<M; k++){
                m_3[i][j] = m_1[i][k] * m_2[k][j];
            }
        }
    }
    t_1 = omp_get_wtime();
    t = t_1 - t_0;
    printf("Petla wewnetrzna \n #pragma omp parallel for private(j,k) schedule(dynamic):\n czas: %f\n\n", t);

//    for(i=0;i<N;i++){
//        for(j=0;j<P;j++){
//            std::cout << m_3[i][j];
//        }
//        std:: cout << "\n";
//    }
//
//    for( i=0;i<N;i++){
//        for( j=0;j<P;j++){
//            if(m_check[i][j] == m_3[i][j]){
//                printf(" | ok | ");
//            }
//        }
//        std:: cout << "\n";
//    }

    std:: cout << "\n\n";
    return 0;
////    return ts;

}

float* generate(int n){
    auto* data = new float[n];

    for (int i = 0; i < n; i++) {
        data[i] = float(rand() % 100) / 5;
    }

    return data;
}

int main(){

//    zad_2(2, 3, 4);
//    double t1[] = zad_2(2, 3, 4);
//    double t3[] = zad_2(3, 4, 5);
//    double t5[] = zad_2(5, 6, 7);
//
//    std::cout << t1[0];
//    std::cout << t1[1];
//    std::cout << t3[0];
//    std::cout << t3[1];
//    std::cout << t5[0];
//    std::cout << t5[1];

float* data1 = generate(1000000);
float* data2 = generate(2000000);
float* data3 = generate(3000000);
//for(int i = 0; i < 10; i++){
//    std::cout << data[i] << " ";
//}



    return 0;
}
