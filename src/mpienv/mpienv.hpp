#ifndef MPI_ENV_H_
#define MPI_ENV_H_

#include "mpi.h"
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

public:
    void init(int argc, char **argv);
    void finalize();
    void File_open(char *File_name);
	void print(const char *message);
};

#endif
