//
// Created by dinghua.li on 6/30/19.
//

#ifndef MEGAHIT_EDGE_COUNTER_H
#define MEGAHIT_EDGE_COUNTER_H

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include "sdbg/sdbg_def.h"

class EdgeMultiplicityRecorder {
 public:
  EdgeMultiplicityRecorder() = default;

  void SetNumThreads(unsigned n) {
    counters_.resize(n);
    for (auto &counter : counters_) {
      std::fill(counter.begin(), counter.end(), 0);
    }
    std::fill(local_counter_sum_.begin(), local_counter_sum_.end(), 0);
    std::fill(local_counter_sum_reduce_.begin(), local_counter_sum_reduce_.end(), 0);
  }

  size_t size_in_byte() const {
    return (kMaxMul + 1) * (counters_.size() + 2) * sizeof(int64_t);
  }

  template <typename T>
  void Add(T multiplicity, unsigned thread_id) {
    ++counters_[thread_id][std::min(static_cast<T>(kMaxMul), multiplicity)];
  }

  int64_t GetNumSolidEdges(int solid_threshold) const {
    int64_t sum = 0;
    for (const auto &counter : counters_) {
      for (int i = solid_threshold; i <= kMaxMul; ++i) {
        sum += counter[i];
      }
    }
    return sum;
  }

  int64_t GetNumSolidEdges_mpi(int solid_threshold) const {
    int64_t sum = 0;
    for (int i = solid_threshold; i <= kMaxMul; ++i) {
      sum += local_counter_sum_reduce_[i];
    }
    return sum;
  }

  void DumpStat(std::ostream &os) const {
    for (int i = 1; i <= kMaxMul; ++i) {
      os << i << ' ' << local_counter_sum_reduce_[i] << '\n';
    }
  }

  void addlocal() {
    for (int i = 1; i <= kMaxMul; i++)
    {
      for (const auto &counter : counters_) {
        local_counter_sum_[i] += counter[i];
      }
    }
  }

 private:
  std::vector<std::array<int64_t, kMaxMul + 1>> counters_;

 public:
  std::array<int64_t, kMaxMul + 1> local_counter_sum_;
  std::array<int64_t, kMaxMul + 1> local_counter_sum_reduce_;
};

#endif  // MEGAHIT_EDGE_COUNTER_H
