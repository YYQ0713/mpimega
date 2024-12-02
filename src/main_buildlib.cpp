#include "sequence/io/sequence_lib.h"
#include "utils/utils.h"
#include "mpienv/mpienv.hpp"

void DisplayHelp(const char *program) {
  pfprintf(stderr, "Usage {s} <read_lib_file> <out_prefix>\n", program);
}

int main_build_lib(int argc, char **argv) {

  MPIEnviroment mpienv;
  mpienv.init(argc, argv);

  if (mpienv.rank == 0)
  {
    AutoMaxRssRecorder recorder;

    if (argc < 3) {
      DisplayHelp(argv[0]);
      exit(1);
    }

    SequenceLibCollection::Build(argv[1], argv[2]);
  }

  mpienv.finalize();

  return 0;
}