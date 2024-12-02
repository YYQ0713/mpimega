/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <omp.h>
#include <stdio.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include "definitions.h"
#include "sorting/kmer_counter.h"
#include "sorting/read_to_sdbg.h"
#include "sorting/seq_to_sdbg.h"
#include "utils/options_description.h"
#include "utils/utils.h"
#include "mpienv/mpienv.hpp"

int main_kmer_count(int argc, char **argv) {

  MPIEnviroment mpienv;
  mpienv.init(argc, argv);
  
  AutoMaxRssRecorder recorder;
  // parse option
  OptionsDescription desc;
  KmerCounterOption opt;

  desc.AddOption("kmer_k", "k", opt.k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.solid_threshold,
                 "min frequency to output an edge");
  desc.AddOption(
      "host_mem", "", opt.host_mem,
      "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file,
                 "read library configuration file.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file.empty()) {
      throw std::logic_error("No read library configuration file!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder count --input_file fastx_file -o out"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  KmerCounter runner(opt, mpienv);
  runner.Run();

  mpienv.finalize();

  return 0;
}