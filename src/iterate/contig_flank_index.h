//
// Created by vout on 11/23/18.
//

#ifndef MEGAHIT_JUNCION_INDEX_H
#define MEGAHIT_JUNCION_INDEX_H

#include <sdbg/sdbg_def.h>
#include <utils/utils.h>
#include "parallel_hashmap/phmap.h"
#include "sequence/kmer_plus.h"
#include "sequence/sequence_package.h"
#include "utils/mutex.h"
#include "ankerl/unordered_dense.h"
#include "sparsepp/spp.h"
#include "libbloom/bloom.h"
#include "sequence/io/edge/edge_writer.h"
#include <rocksdb/db.h>
#include <rocksdb/iterator.h>
#include <rocksdb/write_batch.h>

template <class KmerType>
class ContigFlankIndex {
 public:
  struct FlankInfo {
    uint64_t ext_seq : 58;
    unsigned ext_len : 6;
    float mul;
  } __attribute__((packed));
  using Flank = KmerPlus<KmerType, FlankInfo>;
  using HashSet = phmap::flat_hash_set<Flank, KmerHash>;
  // using HashSet = spp::sparse_hash_set<Flank, KmerHash>;

 public:
  ContigFlankIndex(unsigned k, unsigned step) : k_(k), step_(step) {}
  size_t size() const { return hash_index_.size(); }

  void FeedBatchContigs(SeqPackage &seq_pkg, const std::vector<float> &mul) {
    SpinLock lock;
#pragma omp parallel for
    for (size_t i = 0; i < seq_pkg.seq_count(); ++i) {
      auto seq_view = seq_pkg.GetSeqView(i);
      size_t seq_len = seq_view.length();
      if (seq_len < k_ + 1) {
        continue;
      }
      for (int strand = 0; strand < 2; ++strand) {
        auto get_jth_char = [&seq_view, strand,
                             seq_len](unsigned j) -> uint8_t {
          uint8_t c = seq_view.base_at(strand == 0 ? j : (seq_len - 1 - j));
          return strand == 0 ? c : 3u ^ c;
        };

        KmerType kmer;
        for (unsigned j = 0; j < k_ + 1; ++j) {
          kmer.ShiftAppend(get_jth_char(j), k_ + 1);
        }
        if (kmer.IsPalindrome(k_ + 1)) {
          continue;
        }

        unsigned ext_len =
            std::min(static_cast<size_t>(step_ - 1), seq_len - (k_ + 1));
        uint64_t ext_seq = 0;
        for (unsigned j = 0; j < ext_len; ++j) {
          ext_seq |= uint64_t(get_jth_char(k_ + 1 + j)) << j * 2;
        }

        {
          std::lock_guard<SpinLock> lk(lock);
          auto res = hash_index_.emplace(kmer, FlankInfo{ext_seq, ext_len});
          if (!res.second) {
            auto old_len = res.first->aux.ext_len;
            auto old_seq = res.first->aux.ext_seq;
            if (old_len < ext_len ||
                (old_len == ext_len && old_seq < ext_seq)) {
              hash_index_.erase(res.first);
              res = hash_index_.emplace(kmer, FlankInfo{ext_seq, ext_len});
              assert(res.second);
            }
          }
        }
        if (seq_len == k_ + 1) {
          break;
        }
      }
    }
  }

  template <class CollectorType, class WriterType>
  size_t FindNextKmersFromReads(const SeqPackage &seq_pkg,
                                CollectorType *out, int rank, int nprocs, WriterType *mpiwiriter, int64_t *num_edges, Bloom *blf) const {
    std::vector<bool> kmer_exist;
    std::vector<float> kmer_mul;
    size_t num_aligned_reads = 0;
    size_t num_iter_edges = 0;

    KmerHash hasher;

  #pragma omp parallel
  {
    //rocksdb::WriteBatch thread_batch;

    #pragma omp for reduction(+ : num_aligned_reads, num_iter_edges) private(kmer_exist, kmer_mul)
    for (unsigned seq_id = 0; seq_id < seq_pkg.seq_count(); ++seq_id) {
      auto seq_view = seq_pkg.GetSeqView(seq_id);
      size_t length = seq_view.length();
      if (length < k_ + step_ + 1) {
        continue;
      }

      bool success = false;
      kmer_exist.clear();
      kmer_exist.resize(length, false);
      kmer_mul.clear();
      kmer_mul.resize(length, 0);

      Flank flank, rflank;
      auto &kmer = flank.kmer;
      auto &rkmer = rflank.kmer;

      for (unsigned j = 0; j < k_ + 1; ++j) {
        kmer.ShiftAppend(seq_view.base_at(j), k_ + 1);
      }
      rkmer = kmer;
      rkmer.ReverseComplement(k_ + 1);

      unsigned cur_pos = 0;
      while (cur_pos + k_ + 1 <= length) {
        unsigned next_pos = cur_pos + 1;

        if (!kmer_exist[cur_pos]) {
          auto iter = hash_index_.find(flank);
          if (iter != hash_index_.end()) {
            kmer_exist[cur_pos] = true;
            uint64_t ext_seq = iter->aux.ext_seq;
            unsigned ext_len = iter->aux.ext_len;
            float mul = iter->aux.mul;
            kmer_mul[cur_pos] = mul;

            for (unsigned j = 0; j < ext_len && cur_pos + k_ + 1 + j < length;
                 ++j, ++next_pos) {
              if (seq_view.base_at(cur_pos + k_ + 1 + j) ==
                  ((ext_seq >> j * 2u) & 3u)) {
                kmer_exist[cur_pos + j + 1] = true;
                kmer_mul[cur_pos + j + 1] = mul;
              } else {
                break;
              }
            }
          }
          if ((iter = hash_index_.find(rflank)) != hash_index_.end()) {
            uint64_t ext_seq = iter->aux.ext_seq;
            unsigned ext_len = iter->aux.ext_len;
            float mul = iter->aux.mul;
            kmer_mul[cur_pos] =
                kmer_exist[cur_pos] ? (kmer_mul[cur_pos] + mul) / 2 : mul;
            kmer_exist[cur_pos] = true;

            for (unsigned j = 0; j < ext_len && cur_pos >= j + 1; ++j) {
              if ((3u ^ seq_view.base_at(cur_pos - 1 - j)) ==
                  ((ext_seq >> j * 2u) & 3u)) {
                kmer_mul[cur_pos - 1 - j] =
                    kmer_exist[cur_pos - 1 - j]
                        ? (kmer_mul[cur_pos - 1 - j] + mul) / 2
                        : mul;
                kmer_exist[cur_pos - 1 - j] = true;
              } else {
                break;
              }
            }
          }
        }

        if (next_pos + k_ + 1 <= length) {
          while (cur_pos < next_pos) {
            ++cur_pos;
            uint8_t c = seq_view.base_at(cur_pos + k_);
            kmer.ShiftAppend(c, k_ + 1);
            rkmer.ShiftPreappend(3u ^ c, k_ + 1);
          }
        } else {
          break;
        }
      }

      for (unsigned j = 1; j + k_ + 1 <= length; ++j) {
        kmer_mul[j] += kmer_mul[j - 1];
      }

      typename CollectorType::kmer_type new_kmer, new_rkmer, final_kmer;

      for (unsigned accumulated_len = 0, j = 0, end_pos = 0; j + k_ < length;
           ++j) {
        accumulated_len = kmer_exist[j] ? accumulated_len + 1 : 0;
        if (accumulated_len >= step_ + 1) {
          unsigned target_end = j + k_ + 1;
          if (end_pos + 8 < target_end) {
            while (end_pos < target_end) {
              auto c = seq_view.base_at(end_pos);
              new_kmer.ShiftAppend(c, k_ + step_ + 1);
              new_rkmer.ShiftPreappend(3u ^ c, k_ + step_ + 1);
              end_pos++;
            }
          } else {
            if (end_pos + k_ + step_ + 1 < target_end) {
              end_pos = target_end - (k_ + step_ + 1);
            }
            while (end_pos < target_end) {
              auto c = seq_view.base_at(end_pos);
              new_kmer.ShiftAppend(c, k_ + step_ + 1);
              end_pos++;
            }
            new_rkmer = new_kmer;
            new_rkmer.ReverseComplement(k_ + step_ + 1);
          }
          float mul =
              (kmer_mul[j] - (j >= step_ + 1 ? kmer_mul[j - (step_ + 1)] : 0)) /
              (step_ + 1);
          assert(mul <= kMaxMul + 1);

          final_kmer = new_kmer < new_rkmer ? new_kmer : new_rkmer;
          auto hash_result = hasher(final_kmer);

          if (hash_result % nprocs == rank) {
            // out->Insert(final_kmer,
            //             static_cast<mul_t>(
            //                 std::min(kMaxMul, static_cast<int>(mul + 0.5))));
            mul_t value_to_write = std::min(kMaxMul, static_cast<int>(mul + 0.5));
            if (!blf->bloom_check_add(final_kmer)) {
              mpiwiriter->WriteToBuf(final_kmer, value_to_write);
              num_iter_edges++;
            } else {
              // rocksdb::Slice key_slice(reinterpret_cast<const char*>(final_kmer.data()), sizeof(final_kmer));
              //mul_t value_to_write = std::min(kMaxMul, static_cast<int>(mul + 0.5));
              // rocksdb::Slice value_slice(reinterpret_cast<const char*>(&value_to_write), sizeof(mul_t));
  
              //rocksdb::Status write_status = db->Put(write_options, key_slice, value_slice);
              // thread_batch.Put(key_slice, value_slice);
              out->Insert(final_kmer);
            }
            

            // 检查批次大小，如果达到阈值则写入
            // if (thread_batch.Count() > 4194304) { // 例如，4M次操作
            //     // 提交线程的批次
            //     rocksdb::Status write_status = db->Write(write_options, &thread_batch);
                
            //     // 提交后清空，准备新的批次
            //     thread_batch.Clear();
            // }

            // if (!blf->bloom_check_add(final_kmer)) {
            //   mpiwiriter->WriteToBuf(final_kmer, static_cast<mul_t>(std::min(kMaxMul, static_cast<int>(mul + 0.5))));
            //   num_iter_edges++;
            // }
          }
          success = true;
        }
      }
      num_aligned_reads += success;

      //仅主线程检查并刷新 buffer
      if (omp_get_thread_num() == 0 && mpiwiriter->check_buf()) {
        mpiwiriter->MPIFileWrite();
      }
    }

    // 每个线程检查自己的批次是否为空，然后提交
    // if (thread_batch.Count() != 0) {
    //     db->Write(write_options, &thread_batch);
    //     thread_batch.Clear();
    // }
  }//#pragma omp parallel

    // xinfo("finish db construct and start writing file\n");
    // rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
    // for (it->SeekToFirst(); it->Valid(); it->Next()) {
    //   typename CollectorType::kmer_type tmp_kmer(it->key().data());
    //   mul_t mul = *reinterpret_cast<const mul_t*>(it->value().data());
    //   mpiwiriter->WriteToBuf(tmp_kmer, mul);
    //   num_iter_edges++;
    //   if (mpiwiriter->check_buf()) {
    //     mpiwiriter->MPIFileWrite();
    //   }
    // }

    // delete db;

    // rocksdb::Options cleanup_options;
    // status = rocksdb::DestroyDB(db_path, cleanup_options);
    // if (!status.ok()) {
    //    std::cerr << "Failed to clean up database: " << status.ToString() << std::endl;
    // }

    *num_edges += num_iter_edges;
    return num_aligned_reads;
  }
 private:
  HashSet hash_index_;
  unsigned k_{};
  unsigned step_{};
};

#endif  // MEGAHIT_JUNCION_INDEX_H
