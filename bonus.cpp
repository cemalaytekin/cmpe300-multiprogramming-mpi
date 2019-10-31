/*
Student Name: Cemal Aytekin
Student Number: 2015400126
Compile Status: Compiling
Program Status: Working
Notes: This code written according to second approach.(BONUS) Explanation of variables, functions and algorithms are commented in the code.
*/

#include <stdio.h>
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <iterator>
#include <vector>

using namespace std;

// It splits the line to the words. Helper method for reading data from the input file
template<class Container>
void split1(const string &str, Container &cont) {
    istringstream iss(str);
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter(cont));
}


double GAMMA;			// keeps the value of gamma
double BETA;			// keeps the value of beta
double PI;			// keeps the value of pi
int N;				// number of slaves
int n;
string INPUT_FILE;		// name of the input file
string OUTPUT_FILE;		// name of the output file
int ARRAY_SIZE = 200;		// one-size of the 2D array read from input file
int MAX_ITERATION = 500000;	// maximum iteration

int main(int argc, char* argv[])
{

    srandom((unsigned)time(NULL));

    // take the arguments and assign the values
    INPUT_FILE = argv[1];
    OUTPUT_FILE = argv[2];
    BETA = stod(argv[3]);
    PI = stod(argv[4]);
    GAMMA = log((1-PI)/PI)/2.0;
   
    // create a communication
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    N = size-1;
    n = sqrt(N);
   

    MAX_ITERATION /= N;
        
    // if master process
    if(rank == 0){
	
	ifstream infile(INPUT_FILE);
	ofstream ofile;
	ofile.open(OUTPUT_FILE);

        int arr[ARRAY_SIZE][ARRAY_SIZE];	// array read from input file
	int newarr[ARRAY_SIZE][ARRAY_SIZE];	// array received from slave processes

	// read input file
        for(int i=0; i<ARRAY_SIZE; i++){
            string line;
	    getline(infile, line);
	    vector<string> words;
	    split1(line, words);
	 
	    for(int j=0; j<ARRAY_SIZE; j++)
 		arr[i][j] = stoi(words[j]);    
         
	}
	
	int s_subarr[ARRAY_SIZE/n][ARRAY_SIZE/n];
	int r_subarr[ARRAY_SIZE/n][ARRAY_SIZE/n];

	// send the sub-arrays to the slaves
        for(int i=0; i<N; i++){
	    for(int m=0; m<ARRAY_SIZE/n; m++){
		for( int k=0; k<ARRAY_SIZE/n; k++){
		   s_subarr[m][k] = arr[m+(i/n)*(ARRAY_SIZE/n)][k+(i%n)*(ARRAY_SIZE/n)];
		}
	    }
	    
	    MPI_Send(s_subarr, (ARRAY_SIZE/n)*(ARRAY_SIZE/n), MPI_INT, i+1, 0, MPI_COMM_WORLD);
	}

	
	// receive the final version of the subarrays from slaves
	for(int i=0; i<N; i++){
	    MPI_Recv(r_subarr, (ARRAY_SIZE/n)*(ARRAY_SIZE/n), MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    for(int m=0; m<ARRAY_SIZE/n; m++){
		for( int k=0; k<ARRAY_SIZE/n; k++){
		   newarr[m+(i/n)*(ARRAY_SIZE/n)][k+(i%n)*(ARRAY_SIZE/n)] = r_subarr[m][k];
		}
	    }
	}
	
	// print final array to output file
	for(int i=0; i<ARRAY_SIZE; i++){
	    for(int j=0; j<ARRAY_SIZE; j++){
		ofile<<newarr[i][j]<<" ";
	    }
	    ofile<<endl;
	}
	
    }

    // if slave process 
    else{

        int X[ARRAY_SIZE/n][ARRAY_SIZE/n];		// sub data, received from master
	int Z[ARRAY_SIZE/n][ARRAY_SIZE/n];		// copy of X
	int r_upper[ARRAY_SIZE/n];			// received from upper adjacent process
	int r_lower[ARRAY_SIZE/n];			// received from lower adjacent process
	int r_right[ARRAY_SIZE/n];			// received from right adjacent process
	int r_left[ARRAY_SIZE/n];			// received from left adjacent process
	int s_upper[ARRAY_SIZE/n];			// send to upper adjacent process
	int s_lower[ARRAY_SIZE/n];			// send to lower adjacent process
	int s_right[ARRAY_SIZE/n];			// send to right adjacent process
	int s_left[ARRAY_SIZE/n];			// send to left adjacent process
	int r_left_top;					// received from left-top adjacent process
	int r_left_down;				// received from left-down adjacent process
	int r_right_top;				// received from right-top adjacent process
	int r_right_down;				// received from right-down adjacent process
	int s_left_top;					// send to left-top adjacent process
	int s_left_down;				// send to left-down adjacent process
	int s_right_top;				// send to right-top adjacent process
	int s_right_down;				// send to right-down adjacent process

	// receive the subarray from master process

	MPI_Recv(X, (ARRAY_SIZE/n)*(ARRAY_SIZE/n), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


	// copy array X to Z
	for(int i=0; i<ARRAY_SIZE/n; i++){
	    for(int j=0; j<ARRAY_SIZE/n; j++){
		Z[i][j] = X[i][j];
	    }
	}
   
	// start to iteration
        for(int l=0; l<MAX_ITERATION; l++){
	  
		
	    // send UPPER
	    if(rank>n){
		for(int s=0; s<ARRAY_SIZE/n; s++)
		    s_upper[s] = Z[0][s];		
	    	MPI_Send(s_upper, ARRAY_SIZE/n, MPI_INT, rank-n, 0, MPI_COMM_WORLD);
	    }
	
	    // send DOWN
	    if(rank<=N-n){
		for(int s=0; s<ARRAY_SIZE/n; s++)
		    s_lower[s] = Z[(ARRAY_SIZE/n)-1][s];	
	        MPI_Send(s_lower, ARRAY_SIZE/n, MPI_INT, rank+n, 0, MPI_COMM_WORLD);
	    }

	    // send RIGHT
	    if(rank%n!=0){
		for(int s=0; s<ARRAY_SIZE/n; s++)
		    s_right[s] = Z[s][ARRAY_SIZE/n-1];	
	    	MPI_Send(s_right, ARRAY_SIZE/n, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	    }

	    // send LEFT
	    if(rank%n!=1){
		for(int s=0; s<ARRAY_SIZE/n; s++)
		    s_left[s] = Z[s][0];		
	    	MPI_Send(s_left, ARRAY_SIZE/n, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
	    }

	    // send left-top
	    if(rank%n!=1 && rank>n){
		s_left_top = Z[0][0];
		MPI_Send(&s_left_top, 1, MPI_INT, rank-n-1, 0, MPI_COMM_WORLD);
	    }
		
	    // send left-down
	    if(rank%n!=1 && rank<=N-n){
		s_left_down = Z[(ARRAY_SIZE/n)-1][0];
		MPI_Send(&s_left_down,1, MPI_INT, rank+n-1, 0, MPI_COMM_WORLD);
	    }		

	    // send right-top
	    if(rank%n!=0 && rank>n){
		s_right_top = Z[(ARRAY_SIZE/n)-1][0];
		MPI_Send(&s_right_top,1, MPI_INT, rank-n+1,  0, MPI_COMM_WORLD);
	    }	

	    // send right-down
	    if(rank%n!=0 && rank<=N-n){
		s_right_down = Z[(ARRAY_SIZE/n)-1][(ARRAY_SIZE/n)-1];
		MPI_Send(&s_right_down,1, MPI_INT, rank+n+1, 0, MPI_COMM_WORLD);
	    }	



	    // receive upper
	    if(rank>n)
		MPI_Recv(r_upper, ARRAY_SIZE/n, MPI_INT, rank-n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // receive lower
	    if(rank <= N-n)
		MPI_Recv(r_lower, ARRAY_SIZE/n, MPI_INT, rank+n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // receive right
	    if(rank%n!=0)
		MPI_Recv(r_right, ARRAY_SIZE/n, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // receive left
	    if(rank%n!=1)
		MPI_Recv(r_left, ARRAY_SIZE/n, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	    // receive left-top
	    if(rank%n!=1 && rank>n)
		MPI_Recv(&r_left_top, 1, MPI_INT, rank-n-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // receive left-bottom
	    if(rank%n!=1 && rank<=N-n){
		MPI_Recv(&r_left_down, 1, MPI_INT, rank+n-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }		

	    // receive right-top
	    if(rank%n!=0 && rank>n){
		MPI_Recv(&r_right_top, 1,  MPI_INT, rank-n+1,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }	

	    // receive right-down
	    if(rank%n!=0 && rank<=N-n){
		MPI_Recv(&r_right_down, 1, MPI_INT, rank+n+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }	



 	    int i = rand()%(ARRAY_SIZE/n)+0;
            int j = rand()%(ARRAY_SIZE/n)+0;
	

	    double random = ((double) rand() / (double)RAND_MAX);

	    int sum = 0;			// keeps the sum of the indexes surraounding
	    int value = 0;			// value of that index
            for(int k=-1; k<=1; k++){
		for(int m=-1; m<=1; m++){
			
		    value = 0;
 		    if(m==0 && k==0)
			continue;
		
		    // left-top
		    else if( j+m == -1 && i+k == -1 && rank%n!=1 && rank>n)
			value = r_left_top;

		    // left-bottom
		    else if( j+m == -1 && i+k == ARRAY_SIZE/n && rank%n!=1 && rank<=N-n)
			value = r_left_down;

		    // right-top
		    else if( j+m == ARRAY_SIZE/n && i+k == -1 && rank%n!=0 && rank>n)
			value = r_right_top;

		    // right-bottom
		    else if( j+m == ARRAY_SIZE/n && i+k == ARRAY_SIZE/n && rank%n!=0 && rank<=N-n)
			value = r_right_down;

		    // left bound
		    else if(j+m == -1 && rank%n!=1)
			value = r_left[i+k];

		    // right bound
		    else if(j+m == ARRAY_SIZE/n && rank%n!=0)
			value = r_right[i+k];
		    
    		    // upper bound
		    else if(i+k == -1 && rank>n)
			    value = r_upper[j+m];

		    //lower bound
		    else if(i+k == ARRAY_SIZE/n && rank<= N-n )
			    value = r_lower[j+m]; 

		    else
			value = Z[i+k][j+m];
    
		    sum+=value;	

		}
	    }
	    // calculate the delta_E
	    double delta_E = ((-2*GAMMA*X[i][j]*Z[i][j]) + (-2*BETA*Z[i][j]*sum));
		
	    // flip if necessary
	    if(delta_E > log(random))
		Z[i][j] = -Z[i][j];
	      
	    
	}

	// send the new subarray to the master
	MPI_Send(Z, (ARRAY_SIZE/n)*(ARRAY_SIZE/n), MPI_INT, 0, 0, MPI_COMM_WORLD);	
    		
    }

    // finish the communication
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
