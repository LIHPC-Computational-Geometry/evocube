#include <iostream>
#include <chrono>
#include <omp.h>

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds; 

static long num_steps = 10000000;
double step;

int main(){
    
    auto main_start = steady_clock::now();
    //omp_set_num_threads(4);

    int i;
    double x, pi=0.0;
    double sum = 0;
    step = 1.0/(double) num_steps;

    int NUM_THREADS = omp_get_num_procs();

    omp_set_num_threads(NUM_THREADS);
    std::cout << "# of threads: " << omp_get_num_threads() << std::endl;

    //int id = omp_get_thread_num();
    //printf("hello(%d) ", id);
    //printf("world(%d)\n", id);
    //int nthrds = omp_get_num_threads();


#pragma omp parallel for reduction(+:sum)
    for (i = 0; i<num_steps; i+=1){
        x = (i + 0.5) * step;
        sum += 4.0 / (1.0 + x * x);

    }

//#pragma omp atomic
    

    pi = step * sum;


    printf("pi(%f) ", pi);

    auto main_end = steady_clock::now();
    std::cout << "Total time (μs): "<< duration_cast<microseconds>(main_end - main_start).count() << std::endl;
    return 0;
}