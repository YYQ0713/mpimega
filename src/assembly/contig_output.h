//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_OUTPUT_H
#define MEGAHIT_CONTIG_OUTPUT_H

#include <cassert>
#include <cstdint>
#include <string>

#include "sequence/io/contig/contig_writer.h"
#include "mpienv/mpienv.hpp"

class UnitigGraph;

void OutputContigs(UnitigGraph &graph, ContigWriter *contig_writer,
                   ContigWriter *final_contig_writer, bool change_only,
                   uint32_t min_standalone);
void MPIOutputContigs(UnitigGraph &graph, MPIContigWriter *contig_writer,
                   MPIContigWriter *final_contig_writer, bool change_only,
                   uint32_t min_standalone, MPIEnviroment &mpienv);
#endif  // MEGAHIT_CONTIG_OUTPUT_H
