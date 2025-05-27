//
// Created by vout on 11/21/18.
//

#include "weak_link_remover.h"
#include "unitig_graph.h"

uint32_t DisconnectWeakLinks(UnitigGraph &graph, MPIEnviroment &mpienv, double local_ratio = 0.1) {
  uint32_t num_disconnected = 0;
  int64_t num_edges_mean = graph.size() / mpienv.nprocs;
  int64_t remain = graph.size() % mpienv.nprocs;
  int64_t start_index = mpienv.rank * num_edges_mean + (mpienv.rank < remain ? mpienv.rank : remain);
  int64_t end_index = start_index + num_edges_mean + (mpienv.rank < remain ? 1 : 0);

#pragma omp parallel for reduction(+ : num_disconnected)
  for (UnitigGraph::size_type i = start_index; i < end_index; ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    if (adapter.IsStandalone() || adapter.IsPalindrome()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      UnitigGraph::VertexAdapter nexts[4];
      double depths[4];
      double total_depth = 0;
      int degree = graph.GetNextAdapters(adapter, nexts);
      if (degree <= 1) {
        continue;
      }
      for (int j = 0; j < degree; ++j) {
        depths[j] = nexts[j].GetAvgDepth();
        total_depth += depths[j];
      }
      for (int j = 0; j < degree; ++j) {
        if (depths[j] <= local_ratio * total_depth) {
          num_disconnected += nexts[j].SetToDisconnect();
        }
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &num_disconnected, 1, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);
  graph.Mpi_Allreduce_vertices();
  graph.MPIRefresh(false, mpienv.rank);
  return num_disconnected;
}