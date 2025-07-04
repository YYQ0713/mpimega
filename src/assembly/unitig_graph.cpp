//
// Created by vout on 11/19/18.
//

#include "unitig_graph.h"
#include <omp.h>
#include <cmath>

#include "kmlib/kmbitvector.h"
#include "utils/mutex.h"
#include "utils/utils.h"

UnitigGraph::UnitigGraph(SDBG *sdbg, MPIEnviroment &mpienv)
    : sdbg_(sdbg), mpienv_(mpienv), adapter_impl_(this), sudo_adapter_impl_(this) {
  id_map_.clear();
  vertices_.clear();
  loop_vertices_.clear();
  SpinLock path_lock;
  AtomicBitVector locks(sdbg_->size());
  size_t count_palindrome = 0;

// assemble simple paths
#pragma omp parallel for reduction(+ : count_palindrome)
  for (uint64_t edge_idx = mpienv.rank; edge_idx < sdbg_->size(); edge_idx += mpienv.nprocs) {
    if (sdbg_->IsValidEdge(edge_idx) &&
        sdbg_->NextSimplePathEdge(edge_idx) == SDBG::kNullID &&
        locks.try_lock(edge_idx)) {
      bool will_be_added = true;
      uint64_t cur_edge = edge_idx;
      uint64_t prev_edge;
      int64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
      uint32_t length = 1;

      while ((prev_edge = sdbg_->PrevSimplePathEdge(cur_edge)) !=
             SDBG::kNullID) {
        cur_edge = prev_edge;
        if (!locks.try_lock(cur_edge)) {
          will_be_added = false;
          break;
        }
        depth += sdbg_->EdgeMultiplicity(cur_edge);
        ++length;
      }

      if (!will_be_added) {
        continue;
      }

      uint64_t rc_start = sdbg_->EdgeReverseComplement(edge_idx);
      uint64_t rc_end;
      assert(rc_start != SDBG::kNullID);

      auto rc_tmp = sdbg_->EdgeReverseComplement(cur_edge);
      //if (!locks.try_lock(rc_start)) {
      if (std::max(edge_idx, cur_edge) < std::max(rc_start, rc_tmp)) {
          will_be_added = false;
      } else {
        // lock through the rc path
        locks.try_lock(rc_start);
        uint64_t rc_cur_edge = rc_start;
        rc_end = rc_cur_edge;
        bool extend_full = true;
        while ((rc_cur_edge = sdbg_->NextSimplePathEdge(rc_cur_edge)) !=
               SDBG::kNullID) {
          rc_end = rc_cur_edge;
          if (!locks.try_lock(rc_cur_edge)) {
            extend_full = false;
            break;
          }
        }
        if (!extend_full) {
          rc_end = sdbg_->EdgeReverseComplement(cur_edge);
          assert(rc_end != SDBG::kNullID);
        }
      }

      if (will_be_added) {
        std::lock_guard<SpinLock> lk(path_lock);
        vertices_.emplace_back(cur_edge, edge_idx, rc_start, rc_end, depth,
                               length);
        count_palindrome += cur_edge == rc_start;
      }
    }
  }
  xinfo("Graph size without loops: {}, palindrome: {}\n", vertices_.size(),
        count_palindrome);

  MPI_Allreduce(MPI_IN_PLACE, locks.data_array_.data(), locks.data_array_.size(), MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);

  // assemble looped paths
  std::mutex loop_lock;
  size_t count_loop = 0;
#pragma omp parallel for
  for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
    if (!locks.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
      std::lock_guard<std::mutex> lk(loop_lock);
      if (!locks.at(edge_idx)) {
        uint64_t cur_edge = edge_idx;
        uint64_t rc_edge = sdbg_->EdgeReverseComplement(edge_idx);
        uint64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
        uint32_t length = 0;
        // whether it is marked before entering the loop
        bool rc_marked = locks.at(rc_edge);

        while (!locks.at(cur_edge)) {
          locks.set(cur_edge);
          depth += sdbg_->EdgeMultiplicity(cur_edge);
          ++length;
          cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
          assert(cur_edge != SDBG::kNullID);
        }
        assert(cur_edge == edge_idx);

        if (!rc_marked) {
          uint64_t start = sdbg_->NextSimplePathEdge(edge_idx);
          uint64_t end = edge_idx;
          loop_vertices_.emplace_back(start, end, sdbg_->EdgeReverseComplement(end),
                                 sdbg_->EdgeReverseComplement(start), depth,
                                 length, true);
          count_loop += 1;
        }
      }
    }
  }
  xinfo("Graph size of loops: {}, count_loop: {}\n", loop_vertices_.size(), count_loop);
  sdbg_->FreeMultiplicity();
  UniGather();

  vertices_.reserve(vertices_.size() + loop_vertices_.size());
  // 移动拼接，避免拷贝
  vertices_.insert(vertices_.end(), std::make_move_iterator(loop_vertices_.begin()), std::make_move_iterator(loop_vertices_.end()));
  // swap 释放 v2 的容量（缩容技巧）
  std::vector<UnitigGraphVertex>().swap(loop_vertices_);  // v2 清空且 capacity = 0

  xinfo("Graph size without loops: {}, palindrome: {}\n", vertices_.size(), count_palindrome);

  //vertices_sort();
  //show_info(mpienv.rank);
  //MPI_Barrier(MPI_COMM_WORLD); // Barrier
  //exit(0);

  if (vertices_.size() >= kMaxNumVertices) {
    xfatal(
        "Too many vertices in the unitig graph ({} >= {}), "
        "you may increase the kmer size to remove tons of erroneous kmers.\n",
        vertices_.size(), kMaxNumVertices);
  }
  
  id_map_.reserve(vertices_.size() * 2 - count_palindrome);

  for (size_type i = 0; i < vertices_.size(); ++i) {
    VertexAdapter adapter(vertices_[i]);
    id_map_[adapter.b()] = i;
    id_map_[adapter.rb()] = i;
  }
  assert(vertices_.size() * 2 - count_palindrome >= id_map_.size());
}

void UnitigGraph::RefreshDisconnected() {
  SpinLock mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.IsToDelete() || adapter.IsPalindrome() || adapter.IsLoop()) {
      continue;
    }

    uint8_t to_disconnect = adapter.IsToDisconnect();
    adapter.ReverseComplement();
    uint8_t rc_to_disconnect = adapter.IsToDisconnect();
    adapter.ReverseComplement();

    if (!to_disconnect && !rc_to_disconnect) {
      continue;
    }

    if (adapter.GetLength() <= to_disconnect + rc_to_disconnect) {
      adapter.SetToDelete();
      continue;
    }

    auto old_start = adapter.b();
    auto old_end = adapter.e();
    auto old_rc_start = adapter.rb();
    auto old_rc_end = adapter.re();
    uint64_t new_start, new_end, new_rc_start, new_rc_end;

    if (to_disconnect) {
      new_start = sdbg_->NextSimplePathEdge(old_start);
      new_rc_end = sdbg_->PrevSimplePathEdge(old_rc_end);
      assert(new_start != SDBG::kNullID && new_rc_end != SDBG::kNullID);
      sdbg_->SetInvalidEdge(old_start);
      sdbg_->SetInvalidEdge(old_rc_end);
    } else {
      new_start = old_start;
      new_rc_end = old_rc_end;
    }

    if (rc_to_disconnect) {
      new_rc_start = sdbg_->NextSimplePathEdge(old_rc_start);
      new_end = sdbg_->PrevSimplePathEdge(old_end);
      assert(new_rc_start != SDBG::kNullID && new_end != SDBG::kNullID);
      sdbg_->SetInvalidEdge(old_rc_start);
      sdbg_->SetInvalidEdge(old_end);
    } else {
      new_rc_start = old_rc_start;
      new_end = old_end;
    }

    uint32_t new_length =
        adapter.GetLength() - to_disconnect - rc_to_disconnect;
    uint64_t new_total_depth = lround(adapter.GetAvgDepth() * new_length);
    adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
    adapter.SetLength(new_length);
    adapter.SetTotalDepth(new_total_depth);

    std::lock_guard<SpinLock> lk(mutex);
    if (to_disconnect) {
      id_map_.erase(old_start);
      id_map_[new_start] = i;
    }
    if (rc_to_disconnect) {
      id_map_.erase(old_rc_start);
      id_map_[new_rc_start] = i;
    }
  }
}

void UnitigGraph::MPIRefresh(bool set_changed, int rank) {
  static const uint8_t kDeleted = 0x1;
  static const uint8_t kVisited = 0x2;
  size_t vertices_size = 0;
  RefreshDisconnected();
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.IsToDelete()) {
      continue;
    }
    adapter.SetFlag(kDeleted);
    if (adapter.IsStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.e();
      for (size_t j = 1; j < adapter.GetLength(); ++j) {
        auto prev = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev;
        assert(cur_edge != SDBG::kNullID);
      }
      assert(cur_edge == adapter.b());
      sdbg_->SetInvalidEdge(cur_edge);
      if (adapter.IsPalindrome()) {
        break;
      }
    }
  }

  if (rank == 0) {
    AtomicBitVector locks(size());
  #pragma omp parallel for
    for (size_type i = 0; i < vertices_.size(); ++i) {
      auto adapter = MakeSudoAdapter(i);
      if (adapter.IsStandalone() || (adapter.GetFlag() & kDeleted)) {
        continue;
      }
      for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
        if (PrevSimplePathAdapter(adapter).IsValid()) {
          continue;
        }
        if (!locks.try_lock(i)) {
          break;
        }
        std::vector<SudoVertexAdapter> linear_path;
        for (auto cur = NextSimplePathAdapter(adapter); cur.IsValid();
            cur = NextSimplePathAdapter(cur)) {
          linear_path.emplace_back(cur);
        }

        if (linear_path.empty()) {
          adapter.SetFlag(kVisited);
          break;
        }

        size_type back_id = linear_path.back().UnitigId();
        if (back_id != i && !locks.try_lock(back_id)) {
          if (back_id > i) {
            locks.unlock(i);
            break;
          } else {
            locks.lock(back_id);
          }
        }

        auto new_length = adapter.GetLength();
        auto new_total_depth = adapter.GetTotalDepth();
        adapter.SetFlag(kVisited);

        for (auto &v : linear_path) {
          new_length += v.GetLength();
          new_total_depth += v.GetTotalDepth();
          if (v.canonical_id() != adapter.canonical_id()) v.SetFlag(kDeleted);
        }

        auto new_start = adapter.b();
        auto new_rc_end = adapter.re();
        auto new_rc_start = linear_path.back().rb();
        auto new_end = linear_path.back().e();

        adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
        adapter.SetLength(new_length);
        adapter.SetTotalDepth(new_total_depth);
        if (set_changed) adapter.SetChanged();
        break;
      }
    }

    // looped path
    std::mutex mutex;
  #pragma omp parallel for
    for (size_type i = 0; i < vertices_.size(); ++i) {
      auto adapter = MakeSudoAdapter(i);
      if (!adapter.IsStandalone() && !adapter.GetFlag()) {
        std::lock_guard<std::mutex> lk(mutex);
        if (adapter.GetFlag()) {
          continue;
        }

        uint32_t length = adapter.GetLength();
        uint64_t total_depth = adapter.GetTotalDepth();
        SudoVertexAdapter next_adapter = adapter;
        while (true) {
          next_adapter = NextSimplePathAdapter(next_adapter);
          assert(next_adapter.IsValid());
          if (next_adapter.b() == adapter.b()) {
            break;
          }
          next_adapter.SetFlag(kDeleted);
          length += next_adapter.GetLength();
          total_depth += next_adapter.GetTotalDepth();
        }

        auto new_start = adapter.b();
        auto new_end = sdbg_->PrevSimplePathEdge(new_start);
        auto new_rc_end = adapter.re();
        auto new_rc_start = sdbg_->NextSimplePathEdge(new_rc_end);
        assert(new_start == sdbg_->EdgeReverseComplement(new_rc_end));
        assert(new_end == sdbg_->EdgeReverseComplement(new_rc_start));

        adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
        adapter.SetLength(length);
        adapter.SetTotalDepth(total_depth);
        adapter.SetLooped();
        if (set_changed) adapter.SetChanged();
      }
    }

    vertices_.resize(std::remove_if(vertices_.begin(), vertices_.end(),
                                    [](UnitigGraphVertex &a) {
                                      return SudoVertexAdapter(a).GetFlag() &
                                            kDeleted;
                                    }) - vertices_.begin());
    
    vertices_size = vertices_.size();
  } // rank == 0

  double start_time = 0.0, end_time = 0.0;
  if (rank == 0) {
      start_time = MPI_Wtime();  // 记录开始时间（仅 Rank 0）
  }

  MPI_Bcast(&vertices_size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    vertices_.resize(vertices_size);
  }
  UnitigGraph::Mpi_Bcast_vertices();

  if (rank == 0) {
      end_time = MPI_Wtime();  // 记录结束时间（仅 Rank 0）
      xinfo("MPI Mpi_Bcast_vertices Elapsed time: {}\n", (end_time - start_time));
  }
  size_type num_changed = 0;
#pragma omp parallel for reduction(+ : num_changed)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    assert(adapter.IsStandalone() || adapter.GetFlag());
    adapter.SetFlag(0);
    id_map_.at(adapter.b()) = i;
    id_map_.at(adapter.rb()) = i;
    num_changed += adapter.IsChanged();
  }
}

void UnitigGraph::Refresh(bool set_changed) {
  static const uint8_t kDeleted = 0x1;
  static const uint8_t kVisited = 0x2;
  RefreshDisconnected();
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.IsToDelete()) {
      continue;
    }
    adapter.SetFlag(kDeleted);
    if (adapter.IsStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.e();
      for (size_t j = 1; j < adapter.GetLength(); ++j) {
        auto prev = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev;
        assert(cur_edge != SDBG::kNullID);
      }
      assert(cur_edge == adapter.b());
      sdbg_->SetInvalidEdge(cur_edge);
      if (adapter.IsPalindrome()) {
        break;
      }
    }
  }

  //vertices_sort();
  //show_info(mpienv_.rank);
  AtomicBitVector locks(size());
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (adapter.IsStandalone() || (adapter.GetFlag() & kDeleted)) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      if (PrevSimplePathAdapter(adapter).IsValid()) {
        continue;
      }
      if (!locks.try_lock(i)) {
        break;
      }
      std::vector<SudoVertexAdapter> linear_path;
      for (auto cur = NextSimplePathAdapter(adapter); cur.IsValid();
           cur = NextSimplePathAdapter(cur)) {
        linear_path.emplace_back(cur);
      }

      if (linear_path.empty()) {
        adapter.SetFlag(kVisited);
        break;
      }

      size_type back_id = linear_path.back().UnitigId();
      if (back_id != i && !locks.try_lock(back_id)) {
        if (back_id > i) {
          locks.unlock(i);
          break;
        } else {
          locks.lock(back_id);
        }
      }

      auto new_length = adapter.GetLength();
      auto new_total_depth = adapter.GetTotalDepth();
      adapter.SetFlag(kVisited);

      for (auto &v : linear_path) {
        new_length += v.GetLength();
        new_total_depth += v.GetTotalDepth();
        if (v.canonical_id() != adapter.canonical_id()) v.SetFlag(kDeleted);
      }

      auto new_start = adapter.b();
      auto new_rc_end = adapter.re();
      auto new_rc_start = linear_path.back().rb();
      auto new_end = linear_path.back().e();

      adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
      adapter.SetLength(new_length);
      adapter.SetTotalDepth(new_total_depth);
      if (set_changed) adapter.SetChanged();
      break;
    }
  }

  // looped path
  std::mutex mutex;
#pragma omp parallel for
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.IsStandalone() && !adapter.GetFlag()) {
      std::lock_guard<std::mutex> lk(mutex);
      if (adapter.GetFlag()) {
        continue;
      }

      uint32_t length = adapter.GetLength();
      uint64_t total_depth = adapter.GetTotalDepth();
      SudoVertexAdapter next_adapter = adapter;
      while (true) {
        next_adapter = NextSimplePathAdapter(next_adapter);
        assert(next_adapter.IsValid());
        if (next_adapter.b() == adapter.b()) {
          break;
        }
        next_adapter.SetFlag(kDeleted);
        length += next_adapter.GetLength();
        total_depth += next_adapter.GetTotalDepth();
      }

      auto new_start = adapter.b();
      auto new_end = sdbg_->PrevSimplePathEdge(new_start);
      auto new_rc_end = adapter.re();
      auto new_rc_start = sdbg_->NextSimplePathEdge(new_rc_end);
      assert(new_start == sdbg_->EdgeReverseComplement(new_rc_end));
      assert(new_end == sdbg_->EdgeReverseComplement(new_rc_start));

      adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
      adapter.SetLength(length);
      adapter.SetTotalDepth(total_depth);
      adapter.SetLooped();
      if (set_changed) adapter.SetChanged();
    }
  }

  vertices_.resize(std::remove_if(vertices_.begin(), vertices_.end(),
                                  [](UnitigGraphVertex &a) {
                                    return SudoVertexAdapter(a).GetFlag() &
                                           kDeleted;
                                  }) -
                   vertices_.begin());

  size_type num_changed = 0;
  //vertices_sort();
#pragma omp parallel for reduction(+ : num_changed)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    auto adapter = MakeSudoAdapter(i);
    assert(adapter.IsStandalone() || adapter.GetFlag());
    adapter.SetFlag(0);
    id_map_.at(adapter.b()) = i;
    id_map_.at(adapter.rb()) = i;
    num_changed += adapter.IsChanged();
  }
}

std::string UnitigGraph::VertexToDNAString(VertexAdapter v) {
  v.ToUniqueFormat();
  std::string label;
  label.reserve(k() + v.GetLength());
  uint64_t cur_edge = v.e();

  for (unsigned i = 1; i < v.GetLength(); ++i) {
    int8_t cur_char = sdbg_->GetW(cur_edge);
    label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

    cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
    if (cur_edge == SDBG::kNullID) {
      xfatal("{}, {}, {}, {}, ({}, {}), {}, {}\n", v.b(), v.e(), v.rb(), v.re(),
             sdbg_->EdgeReverseComplement(v.e()),
             sdbg_->EdgeReverseComplement(v.b()), v.GetLength(), i);
    }
  }

  int8_t cur_char = sdbg_->GetW(cur_edge);
  label.push_back("ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)]);

  if (cur_edge != v.b()) {
    xfatal("fwd: {}, {}, rev: {}, {}, ({}, {}) length: {}\n", v.b(), v.e(),
           v.rb(), v.re(), sdbg_->EdgeReverseComplement(v.e()),
           sdbg_->EdgeReverseComplement(v.b()), v.GetLength());
  }

  uint8_t seq[kMaxK];
  sdbg_->GetLabel(v.b(), seq);

  for (int i = sdbg_->k() - 1; i >= 0; --i) {
    assert(seq[i] >= 1 && seq[i] <= 4);
    label.append(1, "ACGT"[seq[i] - 1]);
  }

  std::reverse(label.begin(), label.end());
  return label;
}

void UnitigGraph_Allreduce_Vertices_Op(void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) {
  UnitigGraphVertex* in = static_cast<UnitigGraphVertex*>(invec);
  UnitigGraphVertex* inout = static_cast<UnitigGraphVertex*>(inoutvec);
  
  uint32_t length = static_cast<uint32_t>(*len);

  for (uint32_t i = 0; i < length; ++i) {
    uint8_t flag_in = in[i].flag.v.load(std::memory_order::memory_order_relaxed);
    auto old_val = inout[i].flag.v.fetch_or(flag_in, std::memory_order::memory_order_relaxed);
  }
}

void UnitigGraph::Mpi_Allreduce_vertices() {
  MPI_Datatype MPI_vertices;
  
  // 定义每个字段的长度
  const int kFieldCount = 7;
  int block_lengths[kFieldCount] = {4, 1, 1, 1, 1, 1, 1};

  // 定义每个字段的偏移量
  MPI_Aint displacements[kFieldCount];
  displacements[0] = offsetof(UnitigGraphVertex, strand_info);
  displacements[1] = offsetof(UnitigGraphVertex, total_depth);
  displacements[2] = offsetof(UnitigGraphVertex, length);
  displacements[3] = offsetof(UnitigGraphVertex, is_looped);
  displacements[4] = offsetof(UnitigGraphVertex, is_palindrome);
  displacements[5] = offsetof(UnitigGraphVertex, is_changed);
  displacements[6] = offsetof(UnitigGraphVertex, flag);

  // 定义每个字段的 MPI 数据类型
  MPI_Datatype types[kFieldCount] = {
    MPI_UINT64_T, // strand_info
    MPI_UINT64_T, // total_depth
    MPI_UINT32_T, // length
    MPI_C_BOOL, // is_looped
    MPI_C_BOOL, // is_palindrome
    MPI_C_BOOL, // is_changed
    MPI_UINT8_T, // flag
  };

  // 创建自定义 MPI 类型
  MPI_Type_create_struct(kFieldCount, block_lengths, displacements, types, &MPI_vertices);
  MPI_Type_commit(&MPI_vertices);
  MPI_Op Vertices_Ar_Op;
  MPI_Op_create(UnitigGraph_Allreduce_Vertices_Op, 1, &Vertices_Ar_Op);

  MPI_Allreduce(MPI_IN_PLACE, vertices_.data(), vertices_.size(), MPI_vertices, Vertices_Ar_Op, MPI_COMM_WORLD);

  MPI_Type_free(&MPI_vertices);
  MPI_Op_free(&Vertices_Ar_Op);
}

void UnitigGraph::Mpi_Bcast_vertices() {
  MPI_Datatype MPI_vertices;
  
  // 定义每个字段的长度
  const int kFieldCount = 7;
  int block_lengths[kFieldCount] = {4, 1, 1, 1, 1, 1, 1};

  // 定义每个字段的偏移量
  MPI_Aint displacements[kFieldCount];
  displacements[0] = offsetof(UnitigGraphVertex, strand_info);
  displacements[1] = offsetof(UnitigGraphVertex, total_depth);
  displacements[2] = offsetof(UnitigGraphVertex, length);
  displacements[3] = offsetof(UnitigGraphVertex, is_looped);
  displacements[4] = offsetof(UnitigGraphVertex, is_palindrome);
  displacements[5] = offsetof(UnitigGraphVertex, is_changed);
  displacements[6] = offsetof(UnitigGraphVertex, flag);

  // 定义每个字段的 MPI 数据类型
  MPI_Datatype types[kFieldCount] = {
    MPI_UINT64_T, // strand_info
    MPI_UINT64_T, // total_depth
    MPI_UINT32_T, // length
    MPI_C_BOOL, // is_looped
    MPI_C_BOOL, // is_palindrome
    MPI_C_BOOL, // is_changed
    MPI_UINT8_T, // flag
  };

  // 创建自定义 MPI 类型
  MPI_Type_create_struct(kFieldCount, block_lengths, displacements, types, &MPI_vertices);
  MPI_Type_commit(&MPI_vertices);

  MPI_Bcast(vertices_.data(), vertices_.size(), MPI_vertices, 0, MPI_COMM_WORLD);

  MPI_Type_free(&MPI_vertices);
}

void UnitigGraph::UniGather() {
  MPI_Datatype MPI_vertices;

  const int kFieldCount = 7;
  int block_lengths[kFieldCount] = {4, 1, 1, 1, 1, 1, 1};
  MPI_Aint displacements[kFieldCount];
  displacements[0] = offsetof(UnitigGraphVertex, strand_info);
  displacements[1] = offsetof(UnitigGraphVertex, total_depth);
  displacements[2] = offsetof(UnitigGraphVertex, length);
  displacements[3] = offsetof(UnitigGraphVertex, is_looped);
  displacements[4] = offsetof(UnitigGraphVertex, is_palindrome);
  displacements[5] = offsetof(UnitigGraphVertex, is_changed);
  displacements[6] = offsetof(UnitigGraphVertex, flag);

  MPI_Datatype types[kFieldCount] = {
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT32_T,
    MPI_C_BOOL,
    MPI_C_BOOL,
    MPI_C_BOOL,
    MPI_UINT8_T,
  };

  MPI_Type_create_struct(kFieldCount, block_lengths, displacements, types, &MPI_vertices);
  MPI_Type_commit(&MPI_vertices);

  // Step 1: gather local sizes as MPI_Count
  MPI_Count local_count = static_cast<MPI_Count>(vertices_.size());
  std::vector<MPI_Count> recv_counts(mpienv_.nprocs);
  MPI_Allgather_c(&local_count, 1, MPI_COUNT, recv_counts.data(), 1, MPI_COUNT, MPI_COMM_WORLD);

  // Step 2: compute displacements
  std::vector<MPI_Count> displs(mpienv_.nprocs, 0);
  MPI_Count total_count = 0;
  for (int i = 0; i < mpienv_.nprocs; ++i) {
    displs[i] = total_count;
    total_count += recv_counts[i];
  }

  // Step 3: resize vertices_ to hold all gathered data
  MPI_Count old_local_count = local_count;
  vertices_.resize(total_count);

  // Step 4: move local data to correct offset if using IN_PLACE
  if (displs[mpienv_.rank] != 0) {
    std::memmove(
      vertices_.data() + displs[mpienv_.rank],
      vertices_.data(),
      old_local_count * sizeof(UnitigGraphVertex)
    );
  }

  // Step 5: use IN_PLACE gather
  MPI_Allgatherv_c(
    MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,         // sendbuf ignored
    vertices_.data(), recv_counts.data(), displs.data(),
    MPI_vertices, MPI_COMM_WORLD
  );

  MPI_Type_free(&MPI_vertices);
}

void UnitigGraph::show_info(int rank) {
  std::string filename = std::to_string(rank) + "_output.txt"; // 文件名
    
  // 1. 打开文件（"w" 表示写入，覆盖原有内容）
  FILE* file = fopen(filename.c_str(), "w");
  if (file == NULL) {
      perror("Failed to open file");
      return;
  }

  for (size_t i = 0; i < vertices_.size(); i++)
  {
    fprintf(file, "beginend[%llu]: %llu %llu %llu %llu\n", i, vertices_[i].strand_info[0].begin, vertices_[i].strand_info[0].end, vertices_[i].strand_info[1].begin, vertices_[i].strand_info[1].end);
    fprintf(file, "length[%llu]: %u\n", i, vertices_[i].length);
    fprintf(file, "flag[%llu]: %u\n", i, (unsigned int)vertices_[i].flag.v);
    fprintf(file, "is_changed[%llu]: %u\n", i, (unsigned int)vertices_[i].is_changed);
    fprintf(file, "is_looped[%llu]: %u\n", i, (unsigned int)vertices_[i].is_looped);
    fprintf(file, "is_palindrome[%llu]: %u\n", i, (unsigned int)vertices_[i].is_palindrome);
  }

  // 3. 关闭文件
  fclose(file);
}

void UnitigGraph::vertices_resize(size_t size) {
  vertices_.resize(size);
}

//按照strandinfo的begin升序排序
void UnitigGraph::vertices_sort() {
    tbb::parallel_sort(vertices_.begin(), vertices_.end(),
        [](const UnitigGraphVertex& a, const UnitigGraphVertex& b) { return a.strand_info[0].begin < b.strand_info[0].begin; }
    );
}

size_t UnitigGraph::vertices_size() {
    return vertices_.size();
}

uint32_t UnitigGraph::VerticesIndexWithSdbgId(uint64_t sdbg_id) {
    auto it = std::lower_bound(
        vertices_.begin(), vertices_.end(), sdbg_id,
        [](const UnitigGraphVertex& node, uint64_t value) {
            return node.strand_info[0].begin < value;
        });

    if (it != vertices_.end() && it->strand_info[0].begin == sdbg_id) {
        return static_cast<int>(std::distance(vertices_.begin(), it));
    }

    // 试图查找反向互补
    uint64_t prev_edge;
    uint64_t rc_sdbg_id = sdbg_->EdgeReverseComplement(sdbg_id);
    while ((prev_edge = sdbg_->PrevSimplePathEdge(rc_sdbg_id)) !=
            SDBG::kNullID) {
        rc_sdbg_id = prev_edge;
    }
    //printf("rc_sdbg_id: {%ld};sdbg_id: {%d}\n", rc_sdbg_id, sdbg_id);
    it = std::lower_bound(
        vertices_.begin(), vertices_.end(), rc_sdbg_id,
        [](const UnitigGraphVertex& node, uint64_t value) {
            return node.strand_info[0].begin < value;
        });

    if (it != vertices_.end() && it->strand_info[0].begin == rc_sdbg_id) {
        return static_cast<int>(std::distance(vertices_.begin(), it));
    }

    uint64_t last_sdbg_id = sdbg_->EdgeReverseComplement(rc_sdbg_id);
    while ((prev_edge = sdbg_->PrevSimplePathEdge(last_sdbg_id)) !=
            SDBG::kNullID) {
        last_sdbg_id = prev_edge;
    }
    return last_sdbg_id;
    // 一定不能到这
    assert(false && "BUG: sdbg_id not found in VerticesIndexWithSdbgId()");
    return -1;  // 只为了编译器警告消除
}