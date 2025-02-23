#ifndef MPI_ENV_H_
#define MPI_ENV_H_

#include "mpi.h"
#include "ga.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/resource.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <ctype.h>
#include <cstdlib>
#include <set>
#include <map>
#include <queue> 
#include <deque>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

#define MPI_CHECK(call) { \
    int mpi_err = (call); \
    if (mpi_err != MPI_SUCCESS) { \
        char error_string[256]; \
        int length_of_error_string; \
        MPI_Error_string(mpi_err, error_string, &length_of_error_string); \
        fprintf(stderr, "MPI error in function %s: %s\n", #call, error_string); \
        MPI_Abort(MPI_COMM_WORLD, mpi_err); \
    } \
}

class MPIEnviroment
{
public:
	//process variables
    int rank, nprocs;
	char processor_name[256];
	char hostname[255];
	int  memusage;
    int namelen, rc;
	
	//file variables
    MPI_File cFile;
    MPI_Offset size;
	unsigned long long start_pos, end_pos, datasize;
	unsigned long long read_offset;
	
	//process status
    MPI_Status status;

    //ga param
    int ga_rank, ga_nprocs;

public:
    void init(int argc, char **argv);
    void init_ga(int argc, char **argv);
    void finalize();
    void finalize_ga();
    void File_open(char *File_name);
	void print(const char *message);
};

#endif
