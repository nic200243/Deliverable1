#include <iostream>
#include <vector>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#define N_ITERATIONS 10

using namespace std;


float random (float min, float max);

bool checkSym(const vector<float> & mat, const unsigned int _size);
void matTranspose(const vector<float> & mat, vector<float> & MatTransposed, const unsigned int _size);

bool checkSymImp(const vector<float> & mat, const unsigned int _size);
void matTransposeImp(const vector<float> & mat, vector<float> & MatTransposed, const unsigned int _size);

bool checkSymOMP(const vector<float> & mat, const unsigned int _size);
void matTransposeOMP(const vector<float> & mat, vector<float> & MatTransposed, const unsigned int _size);

void print_details_serial (const int _iteration, const double _time_checkSym_seq, const double _time_transpose_seq, const double _time_checkSym_imp, const double _time_transpose_imp);
void print_details_parallel (const int _iteration, const double _time_checkSym_OMP, const double _time_transpose_OMP);


int main(){
#ifndef _OPENMP
    cout << "This program is not compiled with OpenMP." << endl;
#else
	srand(time(NULL));
	const unsigned int SIZE[4] = {16, 512, 1024, 4096};
	double start, end;
	for (int size=0; size<4; size++){
		cout << "--- MATRIX SIZE N="<< SIZE[size] << " ---" << endl;
		vector<float> M = vector<float>((SIZE[size]*SIZE[size]), 0.0);
	    vector<float> T = vector<float>((SIZE[size]*SIZE[size]), 0.0);
	    double sum_time_checkSym_seq = 0.0;
	    double sum_time_transpose_seq = 0.0;
	    double sum_time_checkSym_imp = 0.0;
	    double sum_time_transpose_imp = 0.0;
	    
	    for (int i=0; i<SIZE[size]; i++) {
	        for (int j=0; j<SIZE[size]; j++) {
	            M[(i*SIZE[size])+j] = random(-100.0, 100.0);
	        }
	    }
	    
	    for(int iteration=0; iteration<N_ITERATIONS; iteration++){
	    	start = omp_get_wtime();
			checkSym(M, SIZE[size]);
			end = omp_get_wtime();
			double time_checkSym_seq = end-start;
			
			start = omp_get_wtime();
		    matTranspose(M, T, SIZE[size]);
		    end = omp_get_wtime();
			double time_transpose_seq = end-start;
			
			start = omp_get_wtime();
			checkSymImp(M, SIZE[size]);
			end = omp_get_wtime();
			double time_checkSym_imp = end-start;
		    
			start = omp_get_wtime();
		    matTransposeImp(M, T, SIZE[size]);
		    end = omp_get_wtime();
			double time_transpose_imp = end-start;
			
			//print_details_serial(iteration, time_checkSym_seq, time_transpose_seq, time_checkSym_imp, time_transpose_imp);
			
			sum_time_checkSym_seq += time_checkSym_seq;
			sum_time_transpose_seq += time_transpose_seq;
			sum_time_checkSym_imp += time_checkSym_imp;
			sum_time_transpose_imp += time_transpose_imp;
		}
		
		double avg_time_checkSym_seq = sum_time_checkSym_seq/N_ITERATIONS;
		double avg_time_transpose_seq = sum_time_transpose_seq/N_ITERATIONS;
		double avg_time_checkSym_imp = sum_time_checkSym_imp/N_ITERATIONS;
		double avg_time_transpose_imp = sum_time_transpose_imp/N_ITERATIONS;
		cout << "---" << endl;
		cout << "Average elapsed execution times of serial implementations evaluated over " << N_ITERATIONS << " iterations." << endl;
		cout << "Avg time_checkSym_seq = " << avg_time_checkSym_seq << "	| " << "Avg time_transpose_seq = " << avg_time_transpose_seq << endl;
		cout << "Avg time_checkSym_imp = " << avg_time_checkSym_imp << "	| " << "Avg time_transpose_imp = " << avg_time_transpose_imp << endl;
		cout << "---" << endl << endl << endl;
		
    	for (int num_threads=1; num_threads<=64; num_threads*=2) {
    	double sum_time_checkSym_OMP = 0.0;
	    double sum_time_transpose_OMP = 0.0;
	    omp_set_num_threads(num_threads);
	    
	        for (int iteration=0; iteration<N_ITERATIONS; iteration++){
	        	start = omp_get_wtime();
				checkSymOMP(M, SIZE[size]);
				end = omp_get_wtime();
				double time_checkSym_OMP = end-start;
			    
				start = omp_get_wtime();
			    matTransposeOMP(M, T, SIZE[size]);
			    end = omp_get_wtime();
				double time_transpose_OMP = end-start;
	        	
	        	//print_details_parallel(iteration, time_checkSym_OMP, time_transpose_OMP);
	        	
	        	sum_time_checkSym_OMP += time_checkSym_OMP;
				sum_time_transpose_OMP += time_transpose_OMP;
			}
			
			double avg_time_checkSym_OMP = sum_time_checkSym_OMP/N_ITERATIONS;
			double avg_time_transpose_OMP = sum_time_transpose_OMP/N_ITERATIONS;
			
			double avg_speedup_checkSym = avg_time_checkSym_seq / avg_time_checkSym_OMP;
			double avg_speedup_transpose = avg_time_transpose_seq / avg_time_transpose_OMP;
        	double avg_efficiency_checkSym = (avg_speedup_checkSym / num_threads) * 100;
        	double avg_efficiency_transpose = (avg_speedup_transpose / num_threads) * 100;
        	
			cout << "---" << endl;
			cout << "Average elapsed execution times of parallel implementations evaluated over " << N_ITERATIONS << " iterations using " << num_threads << " threads." << endl;
			cout << "Avg time_checkSym_OMP = " << avg_time_checkSym_OMP << "	| " << "Avg time_transpose_OMP = " << avg_time_transpose_OMP << endl;
			cout << "Speedup checkSym = " << avg_speedup_checkSym << "		| " << "Speedup transpose = " << avg_speedup_transpose << endl;
			cout << "Efficiency checkSym = " << avg_efficiency_checkSym << "%		| " << "Efficiency transpose = " << avg_efficiency_transpose << "%" << endl;
			cout << "---" << endl << endl;
		}
		
		cout << endl << endl << endl <<endl;
	}
	
#endif

    return 0;
}

float random (float min, float max) {
    return (float) rand() / (float) RAND_MAX * (max - min) + min;
}

bool checkSym(const vector<float> & mat, const unsigned int _size) {
	bool isSym = true;
    for (int i=0; i<_size; i++) {
        for (int j=0; j<_size; j++) {
            if (mat[(i*_size)+j] != mat[(j*_size)+i])
                isSym = false;
        }
    }
    return isSym;
}

void matTranspose(const vector<float> & mat, vector<float> & MatTransposed, const unsigned int _size) {
    for (int i=0; i<_size; i++) {
        for (int j=0; j<_size; j++) {	
            MatTransposed[(i*_size)+j] = mat[(j*_size)+i];
        }
    }
}

#pragma GCC push_options
#pragma GCC optimize("O2", "unroll-loops", "tree-vectorize")
bool checkSymImp(const vector<float> & mat, const unsigned int _size) {
	bool isSym = true;
	#pragma GCC ivdep
	#pragma GCC simd
    for (int i=0; i<_size; i++) {
        for (int j=0; j<_size; j++) {
        	isSym = ( mat[(i*_size)+j] == mat[(j*_size)+i] );
        }
    }
    return isSym;
}

void matTransposeImp(const vector<float> & mat, vector<float> & MatTransposed, const unsigned int _size) {
	#pragma GCC ivdep
	#pragma GCC simd
    for (int i=0; i<_size; i++) {
        for (int j=0; j<_size; j++) {
            MatTransposed[(i*_size)+j] = mat[(j*_size)+i];	
        }
    }
}
#pragma GCC pop_options

bool checkSymOMP(const vector<float> & mat, const unsigned int _size){
	bool isSym = true;
	#pragma omp parallel shared(isSym)
	{
		#pragma omp for collapse(2) reduction(&&:isSym) 
		for (int i=0; i<_size; i++) {
	        for (int j=0; j<_size; j++) {
	            if (mat[(i*_size)+j] != mat[(j*_size)+i])	
	                isSym = false;
	        }
	    }
	}
	return isSym;
}

void matTransposeOMP(const vector<float> & mat, vector<float> & MatTransposed, const unsigned int _size){
	#pragma omp parallel
	{
		#pragma omp for collapse(2) 
		for (int i=0; i<_size; i++) {
	        for (int j=0; j<_size; j++) {	
	            MatTransposed[(i*_size)+j] = mat[(j*_size)+i];
	        }
		}
	}
}

void print_details_serial (const int _iteration, const double _time_checkSym_seq, const double _time_transpose_seq, const double _time_checkSym_imp, const double _time_transpose_imp){
	cout << "-----------------------------------------------------------------" << endl;
	cout << "ITERATION " << _iteration+1 << " OF " << N_ITERATIONS << endl;
	cout << "time_checkSym_seq = " << _time_checkSym_seq << "	|";
	cout << " time_transpose_seq = " << _time_transpose_seq << endl;
	cout << "time_checkSym_imp = " << _time_checkSym_imp << "	|";
	cout << " time_transpose_imp = " << _time_transpose_imp << endl;
}

void print_details_parallel (const int _iteration, const double _time_checkSym_OMP, const double _time_transpose_OMP){
	cout << "-----------------------------------------------------------------" << endl;
	cout << "[ITERATION " << _iteration+1 << " OF " << N_ITERATIONS << "]"<< endl;
	cout << "time_checkSym_OMP = " << _time_checkSym_OMP << "	| ";
	cout << "time_transpose_OMP = " << _time_transpose_OMP << endl;
}

