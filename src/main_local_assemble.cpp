/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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
#include <iostream>
#include <string>
#include <stdexcept>

#include "localasm/local_assemble.h"
#include "utils/options_description.h"
#include "utils/utils.h"
#include "mpienv/mpienv.hpp"

namespace {

LocalAsmOption ParseLocalAsmOptions(int argc, char *argv[]) {
  LocalAsmOption opt;
  OptionsDescription desc;

  desc.AddOption("contig_file", "c", opt.contig_file, "contig file");
  desc.AddOption("lib_file_prefix", "l", opt.lib_file_prefix,
                 "lib file prefix");
  desc.AddOption("kmin", "", opt.kmin, "");
  desc.AddOption("kmax", "", opt.kmax, "");
  desc.AddOption("step", "", opt.step, "");
  desc.AddOption("seed_kmer", "", opt.seed_kmer,
                 "kmer size for seeding alignments");
  desc.AddOption("min_contig_len", "", opt.min_contig_len, "");
  desc.AddOption("min_mapping_len", "", opt.min_mapping_len, "");
  desc.AddOption("sparsity", "", opt.sparsity, "sparsity of hash mapper");
  desc.AddOption("similarity", "", opt.similarity,
                 "alignment similarity threshold");
  desc.AddOption("num_threads", "t", opt.num_threads, "");
  desc.AddOption("output_file", "o", opt.output_file, "");

  try {
    desc.Parse(argc, argv);
    if (opt.contig_file == "") {
      throw std::logic_error("no contig file!");
    }
    if (opt.lib_file_prefix == "") {
      throw std::logic_error("no read file!");
    }
    if (opt.output_file == "") {
      throw std::logic_error("no output file!");
    }
    if (opt.num_threads == 0) {
      opt.num_threads = omp_get_max_threads();
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0]
              << " -c contigs.fa -r reads.fq -o out.local_contig.fa"
              << std::endl;
    std::cerr << "options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
  return opt;
}
}  // namespace

void ParseLAOptions(int argc, char **argv, LocalAsmOption &opt) {
  for (int i = 1; i < argc; ++i) {
    std::string option = argv[i];
    if (option == "-c" || option == "--contig_file") {
      if (i + 1 <= argc) {
        opt.contig_file = argv[++i];
      }
    }
    else if (option == "-l" || option == "--lib_file_prefix") {
      if (i + 1 <= argc) {
        opt.lib_file_prefix = argv[++i];
      }
    }
    else if (option == "--kmin") {
      if (i + 1 <= argc) {
        opt.kmin = std::stoi(argv[++i]);
      }
    }
    else if (option == "--kmax") {
      if (i + 1 <= argc) {
        opt.kmax = std::stoi(argv[++i]);
      }
    }
    else if (option == "--step") {
      if (i + 1 <= argc) {
        opt.step = std::stoi(argv[++i]);
      }
    }
    else if (option == "--seed_kmer") {
      if (i + 1 <= argc) {
        opt.seed_kmer = std::stoi(argv[++i]);
      }
    }
    else if (option == "--min_contig_len") {
      if (i + 1 <= argc) {
        opt.min_contig_len = std::stoi(argv[++i]);
      }
    }
    else if (option == "--min_mapping_len") {
      if (i + 1 <= argc) {
        opt.min_mapping_len = std::stoi(argv[++i]);
      }
    }
    else if (option == "--sparsity") {
      if (i + 1 <= argc) {
        opt.sparsity = std::stoi(argv[++i]);
      }
    }
    else if (option == "--similarity") {
      if (i + 1 <= argc) {
        opt.similarity = std::stod(argv[++i]);
      }
    }
    else if (option == "-t" || option == "--num_threads") {
      if (i + 1 <= argc) {
        opt.num_threads = std::stoi(argv[++i]);
      }
    }
    else if (option == "-o" || option == "--output_file") {
      if (i + 1 <= argc) {
        opt.output_file = argv[++i];
      }
    }
  }

  try {
    if (opt.contig_file == "") {
      throw std::logic_error("no contig file!");
    }
    if (opt.lib_file_prefix == "") {
      throw std::logic_error("no read file!");
    }
    if (opt.output_file == "") {
      throw std::logic_error("no output file!");
    }
    if (opt.num_threads == 0) {
      opt.num_threads = omp_get_max_threads();
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0]
              << " -c contigs.fa -r reads.fq -o out.local_contig.fa"
              << std::endl;
    //std::cerr << "options:" << std::endl;
    //std::cerr << desc << std::endl;
    exit(1);
  }
};

int main_local(int argc, char **argv, MPIEnviroment &mpienv) {
  AutoMaxRssRecorder recorder;
  LocalAsmOption opt;
  ParseLAOptions(argc, argv, opt);
  //auto opt = ParseLocalAsmOptions(argc, argv);
  omp_set_num_threads(opt.num_threads);
  RunLocalAssembly(opt, mpienv);
  return 0;
}