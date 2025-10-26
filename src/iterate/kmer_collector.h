//
// Created by vout on 11/24/18.
//

#ifndef MEGAHIT_SPANNING_KMER_COLLECTOR_H
#define MEGAHIT_SPANNING_KMER_COLLECTOR_H

#include "parallel_hashmap/phmap.h"
#include "sdbg/sdbg_def.h"
#include "sequence/io/edge/edge_writer.h"
#include "sequence/kmer_plus.h"
#include "utils/mutex.h"

template <class KmerType>
class KmerCollector {
 public:
  using kmer_type = KmerType;
  //using kmer_plus = KmerPlus<KmerType, mul_t>;
  using kmer_plus = KmerPlus<KmerType, bool>;
  // using hash_set = phmap::parallel_flat_hash_set<
  //     kmer_plus, KmerHash,
  //     phmap::container_internal::hash_default_eq<kmer_plus>,
  //     phmap::container_internal::Allocator<kmer_plus>, 12, SpinLock>;
  using hash_set = phmap::parallel_flat_hash_map<
      KmerType, mul_t, KmerHash,
      phmap::container_internal::hash_default_eq<KmerType>,
      phmap::container_internal::Allocator<phmap::container_internal::Pair<KmerType, bool>>, 12, SpinLock>;

  KmerCollector(unsigned k, const std::string &out_prefix, MPIEnviroment &mpienv)
      : k_(k), output_prefix_(out_prefix), mpienv_(mpienv) {
    last_shift_ = k_ % 16;
    last_shift_ = (last_shift_ == 0 ? 0 : 16 - last_shift_) * 2;
    words_per_kmer_ = DivCeiling(k_ * 2 + kBitsPerMul, 32);
    buffer_.resize(words_per_kmer_);

    // writer_.SetFilePrefix(out_prefix);
    // writer_.SetUnordered();
    // writer_.SetKmerSize(k_ - 1);
    // writer_.InitFiles(mpienv_);
  }

  // void Insert(const KmerType &kmer, mul_t mul) {
    //collection_.insert({kmer, mul});
    //collection_.emplace(kmer, mul);
  // }

  void Insert(const KmerType &kmer, mul_t mul) {
    mul_t multi = mul >> 1;
    auto res = collection_.emplace(kmer, multi);
    if (!res.second) {
      // res.first->second = true;
      res.first->second |= 0x8000;
    }
  }

  const hash_set &collection() const { return collection_; }
  void FlushToFile(MPIEdgeWriter<KmerType> &mpi_edgewiriter, int64_t &num_iterative_edges) {
    for (const auto &item : collection_) {
      if (item.second & 0x8000) {
        mul_t mul = item.second & 0x7fff;
        mpi_edgewiriter.WriteToBuf(item.first, mul);
        num_iterative_edges++;
        if (mpi_edgewiriter.check_buf()) {
          mpi_edgewiriter.MPIFileWrite();
        }
      }
    }
  }

 private:
  // void WriteToFile(const KmerType &kmer, mul_t mul) {
  //   std::fill(buffer_.begin(), buffer_.end(), 0);

  //   auto ptr = buffer_.begin();
  //   uint32_t w = 0;

  //   for (unsigned j = 0; j < k_; ++j) {
  //     w = (w << 2) | kmer.GetBase(k_ - 1 - j);
  //     if (j % 16 == 15) {
  //       *ptr = w;
  //       w = 0;
  //       ++ptr;
  //     }
  //   }

  //   assert(ptr - buffer_.begin() < words_per_kmer_);
  //   *ptr = (w << last_shift_);
  //   assert((buffer_.back() & kMaxMul) == 0);
  //   buffer_.back() |= mul;
  //   //writer_.WriteUnordered(buffer_.data());
  // }

 private:
  unsigned k_;
  std::string output_prefix_;
  hash_set collection_;
  //EdgeWriter writer_;
  unsigned last_shift_;
  unsigned words_per_kmer_;
  std::vector<uint32_t> buffer_;
  MPIEnviroment mpienv_;
};

#endif  // MEGAHIT_SPANNING_KMER_COLLECTOR_H
