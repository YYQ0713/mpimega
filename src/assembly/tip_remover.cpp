//
// Created by vout on 11/21/18.
//

#include "tip_remover.h"
#include "unitig_graph.h"

uint32_t RemoveTips(UnitigGraph &graph, uint32_t max_tip_len, MPIEnviroment &mpienv) {
  uint32_t num_removed = 0;
  int64_t num_edges_mean = graph.size() / mpienv.nprocs;
  int64_t remain = graph.size() % mpienv.nprocs;
  int64_t start_index = mpienv.rank * num_edges_mean + (mpienv.rank < remain ? mpienv.rank : remain);
  int64_t end_index = start_index + num_edges_mean + (mpienv.rank < remain ? 1 : 0);
  for (uint32_t thre = 2; thre < max_tip_len;
       thre = std::min(thre * 2, max_tip_len)) {
#pragma omp parallel for reduction(+ : num_removed)
    for (UnitigGraph::size_type i = start_index; i < end_index; ++i) {
      auto adapter = graph.MakeVertexAdapter(i);
      if (adapter.GetLength() >= thre) {
        continue;
      }
      if (adapter.IsStandalone()) {
        bool success = adapter.SetToDelete();
        assert(success);
        num_removed += success;
      } else {
        UnitigGraph::VertexAdapter nexts[4], prevs[4];
        int outd = graph.GetNextAdapters(adapter, nexts);
        int ind = graph.GetPrevAdapters(adapter, prevs);

        if (ind + outd == 0) {
          bool success = adapter.SetToDelete();
          assert(success);
          num_removed += success;
        } else if (outd == 1 && ind == 0) {
          if (nexts[0].GetAvgDepth() > 8 * adapter.GetAvgDepth()) {
            bool success = adapter.SetToDelete();
            assert(success);
            num_removed += success;
          }
        } else if (outd == 0 && ind == 1) {
          if (prevs[0].GetAvgDepth() > 8 * adapter.GetAvgDepth()) {
            bool success = adapter.SetToDelete();
            assert(success);
            num_removed += success;
          }
        }
      }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, &num_removed, 1, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);
    graph.Mpi_Allreduce_vertices();
    graph.MPIRefresh(false, mpienv.rank);
    if (thre >= max_tip_len) {
      break;
    }
  }
  return num_removed;
}
