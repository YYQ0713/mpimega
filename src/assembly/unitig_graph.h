//
// Created by vout on 11/10/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_H
#define MEGAHIT_UNITIG_GRAPH_H

#include <deque>
#include <limits>
#include "parallel_hashmap/phmap.h"
#include "sparsepp/spp.h"
#include "sdbg/sdbg.h"
#include "unitig_graph_vertex.h"
#include <algorithm>

class AsmOptions {
  public:
      string sdbg_name = "";
      string output_prefix = "out";
      int num_cpu_threads = 0;
      int local_width = 1000;
      int max_tip_len = -1;
      int min_standalone = 200;
      double min_depth = -1;
      bool is_final_round = false;
      int bubble_level = 2;
      int merge_len = 20;
      double merge_similar = 0.98;
      int prune_level = 2;
      double disconnect_ratio = 0.1;
      double low_local_ratio = 0.2;
      int cleaning_rounds = 5;
      bool output_standalone = false;
      bool careful_bubble = false;

      string contig_file() { return output_prefix + ".contigs.fa"; }
      string standalone_file() { return output_prefix + ".final.contigs.fa"; }
      string addi_contig_file() { return output_prefix + ".addi.fa"; }
      string bubble_file() { return output_prefix + ".bubble_seq.fa"; }
};

class UnitigGraph {
 public:
  using Vertex = UnitigGraphVertex;
  using VertexAdapter = UnitigGraphVertex::Adapter;
  using size_type = VertexAdapter::size_type;
  static const size_type kMaxNumVertices =
      std::numeric_limits<size_type>::max() - 1;
  static const size_type kNullVertexID = kMaxNumVertices + 1;

 public:
  explicit UnitigGraph(SDBG *sdbg, MPIEnviroment &mpienv);
  UnitigGraph(const UnitigGraph &) = delete;
  UnitigGraph(const UnitigGraph &&) = delete;
  ~UnitigGraph() = default;
  size_type size() const { return vertices_.size(); }
  size_t k() const { return sdbg_->k(); }

 public:
  void Refresh(bool mark_changed = false);
  void Mpi_Allreduce_vertices();
  void Mpi_Bcast_vertices();
  void UniGather();
  void show_info(int rank);
  void vertices_resize(size_t size);
  // void vertices_sort();
  size_t vertices_size();
  std::string VertexToDNAString(VertexAdapter adapter);
  uint32_t VerticesIndexWithSdbgId(uint64_t sdbg_id); 
  // void init_db(AsmOptions &opt);
  // void OpenReadWrite_db();
  // void OpenReadOnly_db();
  // void Delete_db();
  // bool is_del(uint32_t vtx_id);
  // void Destroy_db();

 public:
  /*
   * Function for VertexAdapter obtaining & traversal
   */
  VertexAdapter MakeVertexAdapter(size_type id, int strand = 0) {
    return adapter_impl_.MakeVertexAdapter(id, strand);
  }
  int GetNextAdapters(VertexAdapter &adapter, VertexAdapter *out) {
    return adapter_impl_.GetNextAdapters(adapter, out);
  }
  int GetPrevAdapters(VertexAdapter &adapter, VertexAdapter *out) {
    return adapter_impl_.GetPrevAdapters(adapter, out);
  }
  int OutDegree(VertexAdapter &adapter) {
    return adapter_impl_.OutDegree(adapter);
  }
  int InDegree(VertexAdapter &adapter) {
    return adapter_impl_.InDegree(adapter);
  }

 private:
  /*
   * Function for SudoVertexAdapter obtaining & traversal
   */
  using SudoVertexAdapter = UnitigGraphVertex::SudoAdapter;
  SudoVertexAdapter MakeSudoAdapter(size_type id, int strand = 0) {
    return sudo_adapter_impl_.MakeVertexAdapter(id, strand);
  }
  int GetNextAdapters(SudoVertexAdapter &adapter, SudoVertexAdapter *out) {
    return sudo_adapter_impl_.GetNextAdapters(adapter, out);
  }
  int GetPrevAdapters(SudoVertexAdapter &adapter, SudoVertexAdapter *out) {
    return sudo_adapter_impl_.GetPrevAdapters(adapter, out);
  }
  int OutDegree(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.OutDegree(adapter);
  }
  int InDegree(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.InDegree(adapter);
  }
  SudoVertexAdapter NextSimplePathAdapter(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.NextSimplePathAdapter(adapter);
  }
  SudoVertexAdapter PrevSimplePathAdapter(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.PrevSimplePathAdapter(adapter);
  }

 private:
  /**
   * A wrapper for operating different types of adapters
   * @tparam AdapterType type of the vertex adapter
   */
  template <class AdapterType>
  class AdapterImpl {
   public:
    AdapterImpl(UnitigGraph *graph) : graph_(graph) {}

   public:
    AdapterType MakeVertexAdapter(size_type id, int strand = 0) {
      return {graph_->vertices_[id], strand, id};
    }
    int GetNextAdapters(AdapterType &adapter, AdapterType *out) {
      uint64_t next_starts[4];
      int degree = graph_->sdbg_->OutgoingEdges(adapter.e(), next_starts);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i] = MakeVertexAdapterWithSdbgId(next_starts[i]);
        }
      }
      return degree;
    }
    int GetPrevAdapters(AdapterType &adapter, AdapterType *out) {
      adapter.ReverseComplement();
      int degree = GetNextAdapters(adapter, out);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i].ReverseComplement();
        }
      }
      adapter.ReverseComplement();
      return degree;
    }
    int OutDegree(AdapterType &adapter) {
      return GetNextAdapters(adapter, nullptr);
    }
    int InDegree(AdapterType &adapter) {
      adapter.ReverseComplement();
      int degree = OutDegree(adapter);
      adapter.ReverseComplement();
      return degree;
    }
    AdapterType NextSimplePathAdapter(AdapterType &adapter) {
      uint64_t next_sdbg_id = graph_->sdbg_->NextSimplePathEdge(adapter.e());
      if (next_sdbg_id != SDBG::kNullID) {
        return MakeVertexAdapterWithSdbgId(next_sdbg_id);
      } else {
        return AdapterType{};
      }
    }
    AdapterType PrevSimplePathAdapter(AdapterType &adapter) {
      adapter.ReverseComplement();
      AdapterType ret = NextSimplePathAdapter(adapter);
      ret.ReverseComplement();
      adapter.ReverseComplement();
      return ret;
    }

   private:
    AdapterType MakeVertexAdapterWithSdbgId(uint64_t sdbg_id) {
      uint32_t id = graph_->id_map_.at(sdbg_id);
      AdapterType adapter(graph_->vertices_[id], 0, id);
      if (adapter.b() != sdbg_id) {
        adapter.ReverseComplement();
      }
      return adapter;
    }

   private:
    UnitigGraph *graph_;
  };

  void RefreshDisconnected();

 private:
  SDBG *sdbg_{};
  MPIEnviroment mpienv_;
  //std::deque<UnitigGraphVertex> vertices_;
  std::vector<UnitigGraphVertex> vertices_;
  //kmlib::AtomicBitVector<uint8_t> vtx_del_flag_;
  //std::vector<UnitigGraphVertex> loop_vertices_;
  phmap::flat_hash_map<uint64_t, size_type> id_map_;
  //spp::sparse_hash_map<uint64_t, size_type> id_map_;
  //spp::sparse_hash_map<uint64_t, uint64_t> strand_map_;
  AdapterImpl<VertexAdapter> adapter_impl_;
  AdapterImpl<SudoVertexAdapter> sudo_adapter_impl_;
  // rocksdb::DB* db_;
  // rocksdb::Options options_;
  // rocksdb::WriteOptions write_options_;
  // std::string db_path_;
};

class PackedRecord {
  public:
    PackedRecord(uint32_t id_,
                uint64_t ns, uint64_t ne,
                uint64_t nrcs, uint64_t nrce,
                uint64_t td, uint32_t len,
                bool loop, bool change)
        : new_start(ns), new_end(ne),
          new_rc_start(nrcs), new_rc_end(nrce),
          total_depth(td), id(id_), length(len),
          flags((loop ? 1 : 0) | (change ? 2 : 0)) {}
    // 用两个 64bit 的字段容纳 4 个 48bit 的成员（共 192bit）
    uint64_t new_start   : 48;
    uint64_t new_end     : 48;
    uint64_t new_rc_start: 48;
    uint64_t new_rc_end  : 48;

    uint64_t total_depth;

    uint32_t id;
    uint32_t length;

    uint8_t flags; // bit 0: is_loop, bit 1: is_change

    bool is_loop()   const { return flags & 1; }
    bool is_change() const { return flags & 2; }
};

#endif  // MEGAHIT_UNITIG_GRAPH_H
