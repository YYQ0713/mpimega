//
// Created by vout on 11/21/18.
//

#include "tip_remover.h"
#include "unitig_graph.h"

uint32_t RemoveTips(UnitigGraph &graph, uint32_t max_tip_len, MPIEnviroment &mpienv) {
  uint32_t num_removed = 0;
  // kmlib::AtomicBitVector<uint8_t> to_delete;

  for (uint32_t thre = 2; thre < max_tip_len;
       thre = std::min(thre * 2, max_tip_len)) {
    // to_delete.reset(graph.size());//reset to_delet
#pragma omp parallel for reduction(+ : num_removed)
    for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
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
    
    //MPI_Allreduce(MPI_IN_PLACE, &num_removed, 1, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);
    // MPI_Allreduce(MPI_IN_PLACE, to_delete.data_array_.data(), to_delete.data_array_.size(), MPI_UINT8_T, MPI_BOR, MPI_COMM_WORLD);

// #pragma omp parallel for reduction(+ : num_removed)
//    for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
//      if (to_delete.at(i)) {
//        auto adapter = graph.MakeVertexAdapter(i);
//        bool success = adapter.SetToDelete();
//        assert(success);
//        num_removed++;
//      }
//    }

    graph.Refresh(false);
    if (thre >= max_tip_len) {
      break;
    }
  }
  return num_removed;
  
//   uint32_t num_removed = 0;
//   for (uint32_t thre = 2; thre < max_tip_len;
//        thre = std::min(thre * 2, max_tip_len)) {
// #pragma omp parallel for reduction(+ : num_removed)
//     for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
//       auto adapter = graph.MakeVertexAdapter(i);
//       if (adapter.GetLength() >= thre) {
//         continue;
//       }
//       if (adapter.IsStandalone()) {
//         bool success = adapter.SetToDelete();
//         assert(success);
//         num_removed += success;
//       } else {
//         UnitigGraph::VertexAdapter nexts[4], prevs[4];
//         int outd = graph.GetNextAdapters(adapter, nexts);
//         int ind = graph.GetPrevAdapters(adapter, prevs);

//         if (ind + outd == 0) {
//           bool success = adapter.SetToDelete();
//           assert(success);
//           num_removed += success;
//         } else if (outd == 1 && ind == 0) {
//           if (nexts[0].GetAvgDepth() > 8 * adapter.GetAvgDepth()) {
//             bool success = adapter.SetToDelete();
//             assert(success);
//             num_removed += success;
//           }
//         } else if (outd == 0 && ind == 1) {
//           if (prevs[0].GetAvgDepth() > 8 * adapter.GetAvgDepth()) {
//             bool success = adapter.SetToDelete();
//             assert(success);
//             num_removed += success;
//           }
//         }
//       }
//     }

//     graph.Refresh(false);
//     if (thre >= max_tip_len) {
//       break;
//     }
//   }
//   return num_removed;
}
