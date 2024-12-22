#include "mpi.h"
#include "mpienv.hpp"


void MPIEnviroment::init(int argc, char **argv)
{
	MPI_CHECK(MPI_Init(&argc, &argv));
	MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
	MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
	MPI_CHECK(MPI_Get_processor_name(processor_name, &namelen));
}

void MPIEnviroment::print(const char *message)
{
	int who = RUSAGE_SELF;
	struct rusage usage;
	int usageErrno;
	usageErrno = getrusage(who, &usage);
	if(usageErrno == EFAULT)  		printf("Error: EFAULT\n");
	else if(usageErrno == EINVAL) 		printf("Error: EINVAL\n");

	printf("Proc:%d (%s)[%ld] -> %s\n", rank, processor_name, usage.ru_maxrss, message); 
}

void MPIEnviroment::finalize()
{
	MPI_CHECK(MPI_Finalize());
}

void MPIEnviroment::File_open(char *File_name)
{
        rc = MPI_File_open(MPI_COMM_WORLD, File_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile);

        if(rc) {
                printf("Proc:%d Unable to open file \"%s\"\n", rank, File_name); fflush(stdout);
                exit(0);
        }
        else    {
                printf("Proc:%d file\"%s\" opened ...\n", rank, File_name);
        }
}