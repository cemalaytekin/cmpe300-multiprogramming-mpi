/*
Student Name: Cemal Aytekin
Student Number: 2015400126
Compile Status: Compiling
Program Status: Working
Notes: This code written as first approach. Explanation of variables, functions and algorithms are commented in the code.
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

	// send the sub-arrays to the slaves
        for(int i=0; i<N; i++)
	    MPI_Send(arr[i*(ARRAY_SIZE/N)], ((ARRAY_SIZE/N)*ARRAY_SIZE), MPI_INT, i+1, 0, MPI_COMM_WORLD);	
	
	// receive the final version of the subarrays from slaves
	for(int i=0; i<N; i++)
	    MPI_Recv(newarr[i*(ARRAY_SIZE/N)], ((ARRAY_SIZE/N)*ARRAY_SIZE), MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

        int X[ARRAY_SIZE/N][ARRAY_SIZE];	// subarray receieved from master
	int Z[ARRAY_SIZE/N][ARRAY_SIZE];	// copy of the X
	int r_upper[ARRAY_SIZE];		// line received from above 
	int r_lower[ARRAY_SIZE];		// line received from below 
	int s_upper[ARRAY_SIZE];		// line sent to above 
	int s_lower[ARRAY_SIZE];		// line sent to below 

	// receive the subarray from master process
	MPI_Recv(X, (ARRAY_SIZE/N)*ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// copy array X to Z
	for(int i=0; i<ARRAY_SIZE/N; i++){
	    for(int j=0; j<ARRAY_SIZE; j++){
		Z[i][j] = X[i][j];
	    }
	}
   
	// start to iteration
        for(int l=0; l<MAX_ITERATION; l++){
	  
            // if rank is not 1 (the top process), then send the top line to the above 
	    if(rank!=1){
	        // copy the top line to the new array
		for(int s=0; s<ARRAY_SIZE; s++){
		    s_upper[s] = Z[0][s];	
		}
	    	MPI_Send(s_upper, ARRAY_SIZE, MPI_INT, rank-1, 0, MPI_COMM_WORLD);

	    }
            // if rank is not N (the bottom process), then send the bottom line to the belove 
	    if(rank!=N){
	        // copy the bottom line to the new array
		for(int s=0; s<ARRAY_SIZE; s++){
		    s_lower[s] = Z[(ARRAY_SIZE/N)-1][s];	
		}
	        MPI_Send(s_lower, ARRAY_SIZE, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

	    }
		
	    // if rank is not 1 then receive the upper line from above.
	    if(rank!=1)
		MPI_Recv(r_upper, ARRAY_SIZE, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    // if rank is not N then receive the lower line from below.
	    if(rank!=N)
		MPI_Recv(r_lower, ARRAY_SIZE, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // generate the random index
 	    int i = rand()%(ARRAY_SIZE/N)+0;
            int j = rand()%ARRAY_SIZE+0;
	
	    // generate a number double between 0 and 1
	    double random = ((double) rand() / (double)RAND_MAX);

	    int sum = 0;			// keeps the sum of the indexes surraounding
	    int value = 0;			// value of that index
            for(int k=-1; k<=1; k++){
		for(int m=-1; m<=1; m++){
			
		    value = 0;
 		    if(m==0 && k==0)
			continue;

		    // left bound
		    else if(j+m == -1)
			value = 0;
		    // right bound
		    else if(j+m == ARRAY_SIZE)
			value = 0;
		    
    		    // upper bound
		    else if(i+k == -1){
			if(rank != 1)
			    value = r_upper[j+m];
			else
			    value = 0;	
		    }
		    //lower bound
		    else if(i+k == (ARRAY_SIZE/N)){
			if(rank!=N)
			    value = r_lower[j+m]; 
			else
			    value = 0;   
		    }
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
	MPI_Send(Z, (ARRAY_SIZE/N)*ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);	
    		
    }

    // finish the communication
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
