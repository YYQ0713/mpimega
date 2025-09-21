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
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

#include "assembly/all_algo.h"
#include "assembly/contig_output.h"
#include "assembly/contig_stat.h"
#include "utils/histgram.h"
#include "utils/options_description.h"
#include "utils/utils.h"
#include "mpienv/mpienv.hpp"
#include "assembly/unitig_graph.h"
using std::string;

namespace {
/*
struct LocalAsmOption {
  string sdbg_name;
  string output_prefix{"out"};
  int num_cpu_threads{0};

  int local_width{1000};
  int max_tip_len{-1};
  int min_standalone{200};
  double min_depth{-1};
  bool is_final_round{false};
  int bubble_level{2};
  int merge_len{20};
  double merge_similar{0.98};
  int prune_level{2};
  double disconnect_ratio{0.1};
  double low_local_ratio{0.2};
  int cleaning_rounds{5};
  bool output_standalone{false};
  bool careful_bubble{false};

  string contig_file() { return output_prefix + ".contigs.fa"; }
  string standalone_file() { return output_prefix + ".final.contigs.fa"; }
  string addi_contig_file() { return output_prefix + ".addi.fa"; }
  string bubble_file() { return output_prefix + ".bubble_seq.fa"; }
} opt;


void ParseAsmOption(int argc, char *argv[]) {
  OptionsDescription desc;

  desc.AddOption("sdbg_name", "s", opt.sdbg_name,
                 "succinct de Bruijn graph name");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads,
                 "number of cpu threads");
  desc.AddOption("max_tip_len", "", opt.max_tip_len,
                 "max length for tips to be removed. -1 for 2k");
  desc.AddOption(
      "min_standalone", "", opt.min_standalone,
      "min length of a standalone contig to output to final.contigs.fa");
  desc.AddOption("bubble_level", "", opt.bubble_level, "bubbles level 0-3");
  desc.AddOption("merge_len", "", opt.merge_len,
                 "merge complex bubbles of length <= merge_len * k");
  desc.AddOption("merge_similar", "", opt.merge_similar,
                 "min similarity of complex bubble merging");
  desc.AddOption("prune_level", "", opt.prune_level,
                 "strength of low local depth contig pruning (0-3)");
  desc.AddOption("disconnect_ratio", "", opt.disconnect_ratio,
                 "ratio threshold for disconnecting contigs");
  desc.AddOption("low_local_ratio", "", opt.low_local_ratio,
                 "ratio to define low depth contigs");
  desc.AddOption("cleaning_rounds", "", opt.cleaning_rounds,
                 "number of rounds of graphs cleaning");
  desc.AddOption("min_depth", "", opt.min_depth,
                 "if prune_level >= 2, permanently remove low local coverage "
                 "unitigs under this threshold");
  desc.AddOption("is_final_round", "", opt.is_final_round,
                 "this is the last iteration");
  desc.AddOption("output_standalone", "", opt.output_standalone,
                 "output standalone contigs to *.final.contigs.fa");
  desc.AddOption("careful_bubble", "", opt.careful_bubble,
                 "remove bubble carefully");

  try {
    desc.Parse(argc, argv);
    std::cout << opt.sdbg_name << std::endl;
    std::cout << argv[2] << std::endl;
    if (opt.sdbg_name.empty()) {
      throw std::logic_error("no succinct de Bruijn graph name!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " -s sdbg_name -o output_prefix"
              << std::endl;
    std::cerr << "options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}
*/

}  // namespace

// class AsmOptions {
//   public:
//       string sdbg_name = "";
//       string output_prefix = "out";
//       int num_cpu_threads = 0;
//       int local_width = 1000;
//       int max_tip_len = -1;
//       int min_standalone = 200;
//       double min_depth = -1;
//       bool is_final_round = false;
//       int bubble_level = 2;
//       int merge_len = 20;
//       double merge_similar = 0.98;
//       int prune_level = 2;
//       double disconnect_ratio = 0.1;
//       double low_local_ratio = 0.2;
//       int cleaning_rounds = 5;
//       bool output_standalone = false;
//       bool careful_bubble = false;

//       string contig_file() { return output_prefix + ".contigs.fa"; }
//       string standalone_file() { return output_prefix + ".final.contigs.fa"; }
//       string addi_contig_file() { return output_prefix + ".addi.fa"; }
//       string bubble_file() { return output_prefix + ".bubble_seq.fa"; }
// };

void ParseAsmOption(int argc, char **argv, AsmOptions &opt) {
  for (int i = 1; i < argc; ++i) {
    std::string option = argv[i];
    if (option == "-s" || option == "--sdbg_name") {
      if (i + 1 <= argc) {
        opt.sdbg_name = argv[++i];
      }
    }
    else if (option == "-o" || option == "--output_prefix") {
      if (i + 1 <= argc) {
        opt.output_prefix = argv[++i];
      }
    }
    else if (option == "-t" || option == "--num_cpu_threads") {
      if (i + 1 <= argc) {
        opt.num_cpu_threads = std::stoi(argv[++i]);
      }
    }
    else if (option == "--max_tip_len") {
      if (i + 1 <= argc) {
        opt.max_tip_len = std::stoi(argv[++i]);
      }
    }
    else if (option == "--min_standalone") {
      if (i + 1 <= argc) {
        opt.min_standalone = std::stoi(argv[++i]);
      }
    }
    else if (option == "--bubble_level") {
      if (i + 1 <= argc) {
        opt.bubble_level = std::stoi(argv[++i]);
      }
    }
    else if (option == "--merge_len") {
      if (i + 1 <= argc) {
        opt.merge_len = std::stoi(argv[++i]);
      }
    }
    else if (option == "--merge_similar") {
      if (i + 1 <= argc) {
        opt.merge_similar = std::stod(argv[++i]);
      }
    }
    else if (option == "--prune_level") {
      if (i + 1 <= argc) {
        opt.prune_level = std::stoi(argv[++i]);
      }
    }
    else if (option == "--disconnect_ratio") {
      if (i + 1 <= argc) {
        opt.disconnect_ratio = std::stod(argv[++i]);
      }
    }
    else if (option == "--low_local_ratio") {
      if (i + 1 <= argc) {
        opt.low_local_ratio = std::stod(argv[++i]);
      }
    }
    else if (option == "--cleaning_rounds") {
      if (i + 1 <= argc) {
        opt.cleaning_rounds = std::stoi(argv[++i]);
      }
    }
    else if (option == "--min_depth") {
      if (i + 1 <= argc) {
        opt.min_depth = std::stod(argv[++i]);
      }
    }
    else if (option == "--is_final_round") {
      if (i + 1 <= argc) {
        opt.is_final_round = true;
      }
    }
    else if (option == "--output_standalone") {
      if (i + 1 <= argc) {
        opt.output_standalone = true;
      }
    }
    else if (option == "--careful_bubble") {
      if (i + 1 <= argc) {
        opt.careful_bubble = true;
      }
    }
  }

  try {
    if (opt.sdbg_name.empty()) {
      throw std::logic_error("no succinct de Bruijn graph name!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " -s sdbg_name -o output_prefix"
              << std::endl;
    //std::cerr << "options:" << std::endl;
    //std::cerr << desc << std::endl;
    exit(1);
  }
};

int main_assemble(int argc, char **argv, MPIEnviroment &mpienv) {
  AutoMaxRssRecorder recorder;
  AsmOptions opt;
  ParseAsmOption(argc, argv, opt);
  SDBG dbg;
  SimpleTimer timer;
  size_t vmrss_kb;
  
  // graph loading
  timer.reset();
  timer.start();
  xinfo("Loading succinct de Bruijn graph: {s}", opt.sdbg_name.c_str());
  dbg.LoadFromFile(opt.sdbg_name.c_str(), opt.num_cpu_threads);
  timer.stop();
  xinfoc("Done. Time elapsed: {}\n", timer.elapsed());
  xinfo("Number of Edges: {}; K value: {}\n", dbg.size(), dbg.k());

  // set cpu threads
  if (opt.num_cpu_threads == 0) {
    opt.num_cpu_threads = omp_get_max_threads();
  }
  omp_set_num_threads(opt.num_cpu_threads);
  xinfo("Number of CPU threads: {}\n", opt.num_cpu_threads);

  // set tip len
  if (opt.max_tip_len == -1) {
    opt.max_tip_len = dbg.k() * 2;
  }
  // set min depth
  if (opt.min_depth <= 0) {
    opt.min_depth = sdbg_pruning::InferMinDepth(dbg);
    xinfo("min depth set to {.3}\n", opt.min_depth);
  }

  // tips removal before building unitig graph
  if (opt.max_tip_len > 0) {
    timer.reset();
    timer.start();
    sdbg_pruning::RemoveTips(dbg, opt.max_tip_len, mpienv);
    timer.stop();
    xinfo("Tips removal done! Time elapsed(sec): {.3}\n", timer.elapsed());
  }

  // construct unitig graph
  timer.reset();
  timer.start();
  UnitigGraph graph(&dbg, mpienv);
  timer.stop();
  xinfo("unitig graph size: {}, time for building: {.3}\n", graph.size(),
  timer.elapsed());

  //CalcAndPrintStat(graph);
  xinfo("sizeof(UnitigGraphVertex): {}\n", sizeof(UnitigGraphVertex));

  // set up bubble
  // ContigWriter bubble_writer(opt.bubble_file());
  MPIContigWriter bubble_writer(opt.bubble_file(), mpienv.rank);
  NaiveBubbleRemover naiver_bubble_remover;
  ComplexBubbleRemover complex_bubble_remover;
  complex_bubble_remover.SetMergeSimilarity(opt.merge_similar)
      .SetMergeLevel(opt.merge_len);
  Histgram<int64_t> bubble_hist;
  if (opt.careful_bubble) {
    naiver_bubble_remover.SetCarefulThreshold(0.2).SetWriter(&bubble_writer);
    complex_bubble_remover.SetCarefulThreshold(0.2).SetWriter(&bubble_writer);
  }

  MPI_Barrier(MPI_COMM_WORLD); // Barrier before graph cleaning
  // graph cleaning
  for (int round = 1; round <= opt.cleaning_rounds; ++round) {
    xinfo("Graph cleaning round {}\n", round);
    bool changed = false;
    if (round > 1) {
      timer.reset();
      timer.start();
      uint32_t num_tips = RemoveTips(graph, opt.max_tip_len, mpienv);
      changed |= num_tips > 0;
      timer.stop();
      xinfo("Tips removed: {}, time: {.3}\n", num_tips, timer.elapsed());
    }
    // remove bubbles
    if (opt.bubble_level >= 1) {
      timer.reset();
      timer.start();
      uint32_t num_bubbles = naiver_bubble_remover.PopBubbles(graph, true, mpienv);
      timer.stop();
      xinfo("Number of bubbles removed: {}, Time elapsed(sec): {.3}\n",
            num_bubbles, timer.elapsed());
      changed |= num_bubbles > 0;
    }

    // remove complex bubbles
    if (opt.bubble_level >= 2) {
      timer.reset();
      timer.start();
      uint32_t num_bubbles = complex_bubble_remover.PopBubbles(graph, true, mpienv);
      timer.stop();
      xinfo("Number of complex bubbles removed: {}, Time elapsed(sec): {}\n",
            num_bubbles, timer.elapsed());
      changed |= num_bubbles > 0;
    }

    // disconnect
    timer.reset();
    timer.start();
    uint32_t num_disconnected =
        DisconnectWeakLinks(graph, mpienv, opt.disconnect_ratio);
    timer.stop();
    xinfo("Number unitigs disconnected: {}, time: {.3}\n", num_disconnected,
          timer.elapsed());
    changed |= num_disconnected > 0;

    // excessive pruning
    uint32_t num_excessive_pruned = 0;
    if (opt.prune_level >= 3) {
      timer.reset();
      timer.start();
      num_excessive_pruned = RemoveLowDepth(graph, opt.min_depth);
      num_excessive_pruned += naiver_bubble_remover.PopBubbles(graph, true, mpienv);
      if (opt.bubble_level >= 2 && opt.merge_len > 0) {
        num_excessive_pruned += complex_bubble_remover.PopBubbles(graph, true, mpienv);
      }
      timer.stop();
      xinfo("Unitigs removed in (more-)excessive pruning: {}, time: {.3}\n",
            num_excessive_pruned, timer.elapsed());
    } else if (opt.prune_level >= 2) {
      timer.reset();
      timer.start();
      RemoveLocalLowDepth(graph, opt.min_depth, opt.max_tip_len,
                          opt.local_width, std::min(opt.low_local_ratio, 0.1),
                          true, &num_excessive_pruned, mpienv);
      timer.stop();
      xinfo("Unitigs removed in excessive pruning: {}, time: {.3}\n",
            num_excessive_pruned, timer.elapsed());
    }
    if (!changed) break;
  }

  //ContigStat stat = CalcAndPrintStat(graph);
  
  // output contigs
  //ContigWriter contig_writer(opt.contig_file());
  MPIContigWriter mpi_contig_writer(opt.contig_file(), mpienv.rank);
  //ContigWriter standalone_writer(opt.standalone_file());
  MPIContigWriter mpi_standalone_writer(opt.standalone_file(), mpienv.rank);

  if (!(opt.is_final_round &&
        opt.prune_level >=
            1)) {  // otherwise output after local low depth pruning
    timer.reset();
    timer.start();
    //OutputContigs(graph, &contig_writer,
    //              opt.output_standalone ? &standalone_writer : nullptr, false,
    //              opt.min_standalone);
    MPIOutputContigs(graph, &mpi_contig_writer,
                  opt.output_standalone ? &mpi_standalone_writer : nullptr, false,
                  opt.min_standalone, mpienv);
    timer.stop();
    xinfo("Time to output: {}\n", timer.elapsed());
  }

  // remove local low depth & output as contigs
  if (opt.prune_level >= 1) {
    //ContigWriter addi_contig_writer(opt.addi_contig_file());
    MPIContigWriter mpi_addi_contig_writer(opt.addi_contig_file(), mpienv.rank);

    timer.reset();
    timer.start();
    uint32_t num_removed = IterateLocalLowDepth(
        graph, opt.min_depth, opt.max_tip_len, opt.local_width,
        opt.low_local_ratio, mpienv, opt.is_final_round);

    uint32_t n_bubbles = 0;
    if (opt.bubble_level >= 2 && opt.merge_len > 0) {
      complex_bubble_remover.SetWriter(nullptr);
      n_bubbles = complex_bubble_remover.PopBubbles(graph, false, mpienv);
      timer.stop();
    }
    xinfo(
        "Number of local low depth unitigs removed: {}, complex bubbles "
        "removed: {}, time: {}\n",
        num_removed, n_bubbles, timer.elapsed());
    CalcAndPrintStat(graph);

    //timer.reset();
    //timer.start();
    //MPI_Bcast(&vtx_size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    //if (mpienv.rank != 0) {
    //  graph.vertices_resize(vtx_size);
    //}
    //graph.Mpi_Bcast_vertices();
    //timer.stop();
    //xinfo("Time to MPI communication after local low depth unitigs removed: {}\n", timer.elapsed());

    timer.reset();
    timer.start();
    if (!opt.is_final_round) {
      //OutputContigs(graph, &addi_contig_writer, nullptr, true, 0);
      MPIOutputContigs(graph, &mpi_addi_contig_writer, nullptr, true, 0, mpienv);
    } else {
      //OutputContigs(graph, &contig_writer,
      //              opt.output_standalone ? &standalone_writer : nullptr, false,
      //              opt.min_standalone);
      MPIOutputContigs(graph, &mpi_contig_writer,
                    opt.output_standalone ? &mpi_standalone_writer : nullptr, false,
                    opt.min_standalone, mpienv);
    }
    timer.stop();
    xinfo("Time to output after local low depth unitigs removed: {}\n", timer.elapsed());
    //auto stat_changed = CalcAndPrintStat(graph, false, true);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  vmrss_kb = getCurrentRSS_kb();
  xinfo("End of asm currentRSS: {} KB\n", vmrss_kb);

  return 0;
}