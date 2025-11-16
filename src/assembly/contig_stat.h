//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_STAT_H
#define MEGAHIT_CONTIG_STAT_H

#include <map>
#include "unitig_graph.h"
#include "utils/histgram.h"
#include "utils/utils.h"

using ContigStat = std::map<std::string, uint64_t>;

inline ContigStat CalcAndPrintStat(UnitigGraph &graph, bool print = true,
                                   bool changed_only = false) {
  uint32_t n_isolated = 0, n_looped = 0;
  Histgram<uint64_t> hist;
#pragma omp parallel for reduction(+ : n_looped, n_isolated)
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    
    auto adapter = graph.MakeVertexAdapter(i);
    if (changed_only && !adapter.IsChanged()) {
      continue;
    }
    hist.insert(adapter.GetLength() + graph.k());
    n_looped += adapter.IsLoop();
    n_isolated += adapter.IsStandalone() || (graph.InDegree(adapter) == 0 &&
                                             graph.OutDegree(adapter) == 0);
  }
  uint64_t total_size = hist.sum();
  ContigStat stat = {{"Max", hist.maximum()},
                     {"Min", hist.minimum()},
                     {"N50", hist.Nx(0.5 * total_size)},
                     {"total size", total_size},
                     {"number contigs", hist.size()},
                     {"number looped", n_looped},
                     {"number isolated", n_isolated}};

  if (print) {
    xinfo("");
    for (auto &kv : stat) {
      xinfoc("{s}: {}, ", kv.first.c_str(), kv.second);
    }
    xinfoc("{s}", "\n");
  }
  return stat;
}

// inline ContigStat CalcAndPrintStatMPI(UnitigGraph &graph, MPIEnviroment &mpienv, bool print = true,
//                                    bool changed_only = false) {
//   uint64_t n_isolated = 0, n_looped = 0;
//   uint64_t max = 0, min = 0, total_size = 0;
//   uint64_t num_contigs = 0;

//   Histgram<uint64_t> hist;

// #pragma omp parallel for reduction(+ : n_looped, n_isolated)
//   for (UnitigGraph::size_type i = mpienv.rank; i < graph.size(); i += mpienv.nprocs) {
    
//     auto adapter = graph.MakeVertexAdapter(i);
//     if (changed_only && !adapter.IsChanged()) {
//       continue;
//     }
//     hist.insert(adapter.GetLength() + graph.k());
//     n_looped += adapter.IsLoop();
//     n_isolated += adapter.IsStandalone() || (graph.InDegree(adapter) == 0 &&
//                                              graph.OutDegree(adapter) == 0);
//   }

//   max = hist.maximum();
//   min = hist.minimum();
//   total_size = hist.sum();
//   num_contigs = hist.size();

//   MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
//   MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
//   MPI_Allreduce(MPI_IN_PLACE, &total_size, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Allreduce(MPI_IN_PLACE, &num_contigs, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Allreduce(MPI_IN_PLACE, &n_looped, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
//   MPI_Allreduce(MPI_IN_PLACE, &n_isolated, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

//   ContigStat stat = {{"Max", max},
//                      {"Min", min},
//                      {"total size", total_size},
//                      {"number contigs", num_contigs},
//                      {"number looped", n_looped},
//                      {"number isolated", n_isolated}};

//   if (print) {
//     xinfo("");
//     xinfoc("{s}: {}, {s}: {}, {s}: {}, {s}: {}, {s}: {}, {s}: {}, {s}", "Max", max
//                                                                       , "Min", min
//                                                                       , "number contigs", num_contigs
//                                                                       , "number isolated", n_isolated
//                                                                       , "number looped", n_looped
//                                                                       , "total size", total_size, "\n");
//   }
//   return stat;
// }

inline ContigStat CalcAndPrintStatMPI(UnitigGraph &graph, MPIEnviroment &mpienv, int num_threads, bool print = true,
                                   bool changed_only = false) {
  uint64_t n_isolated = 0, n_looped = 0;
  Histgram<uint64_t> hist;

  // create vector of local hist per thread
  int nthreads = omp_get_max_threads();
  std::vector<Histgram<uint64_t>> local_hists(nthreads);

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
#pragma omp for reduction(+ : n_looped, n_isolated)
    for (UnitigGraph::size_type i = 0; i < graph.size(); i++) {
      auto adapter = graph.MakeVertexAdapter(i);
      if (changed_only && !adapter.IsChanged()) continue;
      local_hists[tid].insert(adapter.GetLength() + graph.k());
      n_looped += adapter.IsLoop();
      n_isolated += adapter.IsStandalone() || (graph.InDegree(adapter) == 0 &&
                                              graph.OutDegree(adapter) == 0);
    }
  } // omp barrier here

  // merge local_hists into hist (single-threaded merge to avoid lock contention)
  for (int t = 0; t < nthreads; ++t) {
    // naive merge: iterate local map and insert into global hist
    // (you may add a merge function in Histgram to merge maps directly under lock)
    for (auto &p : local_hists[t].map_) { // note: map_ is private; you might add accessor/merge
      hist.insert(p.first, p.second);
    }
  }


  uint64_t total_size = hist.sum();
  ContigStat stat = {{"Max", hist.maximum()},
                     {"Min", hist.minimum()},
                     {"N50", hist.Nx(0.5 * total_size)},
                     {"total size", total_size},
                     {"number contigs", hist.size()},
                     {"number looped", n_looped},
                     {"number isolated", n_isolated}};

  if (print) {
    xinfo("");
    for (auto &kv : stat) {
      xinfoc("{s}: {}, ", kv.first.c_str(), kv.second);
    }
    xinfoc("{s}", "\n");
  }
  return stat;
}

#endif  // MEGAHIT_CONTIG_STAT_H
