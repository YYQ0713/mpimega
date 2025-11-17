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

int main_kmer_count(int argc, char **argv, MPIEnviroment &mpienv) {
  
  //AutoMaxRssRecorder recorder;
  // parse option
  OptionsDescription desc;
  //KmerCounterOption opt;//todo
  KmerCounterOption opt;//todo

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
  // KmerCounter runner(opt);
  runner.Run();

  return 0;
}


int main_seq2sdbg(int argc, char **argv, MPIEnviroment &mpienv) {

  AutoMaxRssRecorder recorder;
  Seq2SdbgOption opt;

  for (int i = 1; i < argc; ++i) {
    std::string option = argv[i];

    auto has_next = [&](void) { return (i + 1) < argc; };
    auto next_arg = [&](void) -> const char* { return has_next() ? argv[i + 1] : nullptr; };

    if (option == "--host_mem") {
      if (!has_next()) { std::cerr << "--host_mem needs a value\n"; exit(1); }
      try { opt.host_mem = std::stod(argv[++i]); }
      catch (...) { std::cerr << "invalid --host_mem value\n"; exit(1); }
    }
    else if (option == "-k" || option == "--kmer_size") {
      if (!has_next()) { std::cerr << option << " needs a value\n"; exit(1); }
      try { opt.k = std::stoi(argv[++i]); }
      catch (...) { std::cerr << "invalid " << option << " value\n"; exit(1); }
    }
    else if (option == "--kmer_from") {
      if (!has_next()) { std::cerr << "--kmer_from needs a value\n"; exit(1); }
      try { opt.k_from = std::stoi(argv[++i]); }
      catch (...) { std::cerr << "invalid --kmer_from value\n"; exit(1); }
    }
    else if (option == "-t" || option == "--num_cpu_threads") {
      if (!has_next()) { std::cerr << option << " needs a value\n"; exit(1); }
      try { opt.n_threads = std::stoi(argv[++i]); }
      catch (...) { std::cerr << "invalid " << option << " value\n"; exit(1); }
    }
    else if (option == "--contig") {
      if (!has_next()) { std::cerr << "--contig needs a value\n"; exit(1); }
      opt.contig = argv[++i];
    }
    else if (option == "--bubble") {
      if (!has_next()) { std::cerr << "--bubble needs a value\n"; exit(1); }
      opt.bubble_seq = argv[++i];
    }
    else if (option == "--addi_contig") {
      if (!has_next()) { std::cerr << "--addi_contig needs a value\n"; exit(1); }
      opt.addi_contig = argv[++i];
    }
    else if (option == "--local_contig") {
      if (!has_next()) { std::cerr << "--local_contig needs a value\n"; exit(1); }
      opt.local_contig = argv[++i];
    }
    else if (option == "--input_prefix") {
      if (!has_next()) { std::cerr << "--input_prefix needs a value\n"; exit(1); }
      opt.input_prefix = argv[++i];
    }
    else if (option == "-o" || option == "--output_prefix") {
      if (!has_next()) { std::cerr << option << " needs a value\n"; exit(1); }
      opt.output_prefix = argv[++i];
    }
    else if (option == "--need_mercy") {
      opt.need_mercy = true;
    }
    else if (option == "--mem_flag") {
      if (!has_next()) { std::cerr << "--mem_flag needs a value\n"; exit(1); }
      try { opt.mem_flag = std::stoi(argv[++i]); }
      catch (...) { std::cerr << "invalid --mem_flag value\n"; exit(1); }
    }
    else {
      std::cerr << "Unknown option: " << option << "\n";
    }
  }


  SeqToSdbg runner(opt, mpienv);
  // SeqToSdbg runner(opt);
  runner.Run();

  return 0;
}

int main_read2sdbg(int argc, char **argv, MPIEnviroment &mpienv) {
  AutoMaxRssRecorder recorder;

  // parse option the same as kmer_count
  OptionsDescription desc;
  Read2SdbgOption opt;

  desc.AddOption("kmer_k", "k", opt.k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.solid_threshold,
                 "min frequency to output an edge");
  desc.AddOption(
      "host_mem", "", opt.host_mem,
      "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file,
                 "input fast[aq] file, can be gzip'ed. \"-\" for stdin.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");
  desc.AddOption("need_mercy", "", opt.need_mercy, "to add mercy edges.");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file.empty()) {
      throw std::logic_error("No input file!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr
        << "Usage: sdbg_builder read2sdbg --read_lib_file fastx_file -o out"
        << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  SeqPkgWithSolidMarker pkg;

  {
    // stage 1
    Read2SdbgS1 runner(opt, &pkg, mpienv);
    // Read2SdbgS1 runner(opt, &pkg);
    if (opt.solid_threshold > 1) {
      runner.Run();
    } else {
      runner.Initialize();
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  {
    // stage 2
    Read2SdbgS2 runner(opt, &pkg, mpienv);
    // Read2SdbgS2 runner(opt, &pkg);
    runner.Run();
  }

  return 0;
}