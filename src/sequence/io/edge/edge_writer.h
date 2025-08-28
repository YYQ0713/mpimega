#ifndef MEGAHIT_EDGE_WRITER_H
#define MEGAHIT_EDGE_WRITER_H

#include <cassert>
#include <cstdint>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "definitions.h"
#include "edge_io_meta.h"
#include "mpienv/mpienv.hpp"
#include "utils/buffered_reader.h"
#include "utils/utils.h"
#include "sdbg/sdbg_def.h"

class EdgeWriter {
 private:
  std::string file_prefix_;
  std::vector<std::unique_ptr<std::ofstream>> files_;
  std::vector<int64_t> n_edges_at_thread_;

  bool is_opened_;

 public:
  EdgeIoMetadata metadata_{};

  class Snapshot {
   private:
    EdgeIoBucketInfo bucket_info;
    int bucket_id{-1};
    friend class EdgeWriter;
  };

  EdgeWriter() : is_opened_(false){};
  ~EdgeWriter() { Finalize(); }

  void SetKmerSize(uint32_t k) {
    metadata_.kmer_size = k;
    metadata_.words_per_edge = DivCeiling((k + 1) * 2 + 16, 32); // (k + 1)-mer and counts(16bits)
  }

  void SetNumThreads(int32_t num_threads) { metadata_.num_files = num_threads; }
  void SetFilePrefix(const std::string &prefix) { file_prefix_ = prefix; }
  void SetNumBuckets(int num_buckets) {
    metadata_.buckets.clear();
    //metadata_.buckets_reduce.clear();
    metadata_.buckets.resize(num_buckets);
    //metadata_.buckets_reduce.resize(num_buckets);
  }
  void SetUnordered() {
    metadata_.buckets.clear();
    metadata_.is_sorted = false;
    metadata_.num_files = 1;
    metadata_.num_edges = 0;
  }

  void InitFiles(MPIEnviroment& mpienv) {
    assert(!is_opened_);
    n_edges_at_thread_.resize(metadata_.num_files, 0);

    for (unsigned i = 0; i < metadata_.num_files; ++i) {
      files_.emplace_back(new std::ofstream(
          (file_prefix_ + ".rank." + std::to_string(mpienv.rank) + ".edges." + std::to_string(i)).c_str(),
          std::ofstream::binary | std::ofstream::out));
    }

    is_opened_ = true;
  }

  // TO DO : 需要映射一下fileid和不同进程下各个线程之间的对应，不然tid对应不到是哪个进程的tid
  void Write(MPIEnviroment& mpienv, int n_threads, uint32_t *edge_ptr, int32_t bucket, int tid,
             Snapshot *snapshot) const {
    assert(metadata_.is_sorted);
    if (bucket != snapshot->bucket_id) {
      if (snapshot->bucket_id != -1) {
        xinfo("bucket: {} | snapshot->bucket_id {} | snapshot->bucket_info.file_id: | {}\n", bucket, snapshot->bucket_id, snapshot->bucket_info.file_id);
      }
      assert(snapshot->bucket_id == -1);
      assert(snapshot->bucket_info.file_id == -1);
      snapshot->bucket_id = bucket;
      snapshot->bucket_info.file_id = mpienv.rank * n_threads + tid;
      snapshot->bucket_info.file_offset = n_edges_at_thread_[tid];
    }
    assert(snapshot->bucket_id == bucket);
    assert(snapshot->bucket_info.file_id == mpienv.rank * n_threads + tid);

    files_[tid]->write(reinterpret_cast<const char *>(edge_ptr),
                       sizeof(uint32_t) * metadata_.words_per_edge);
    ++snapshot->bucket_info.total_number;
  }


  // TO DO : 需要映射一下fileid和不同进程下各个线程之间的对应，不然tid对应不到是哪个进程的tid
  void SaveSnapshot(const Snapshot &snapshot, int tid) {
    if (snapshot.bucket_id != -1) {
      metadata_.buckets[snapshot.bucket_id] = snapshot.bucket_info;
      n_edges_at_thread_[tid] =
          snapshot.bucket_info.total_number + snapshot.bucket_info.file_offset;
    }
  }

  void WriteUnordered(uint32_t *edge_ptr) {
    assert(!metadata_.is_sorted);
    files_[0]->write(reinterpret_cast<const char *>(edge_ptr),
                     sizeof(uint32_t) * metadata_.words_per_edge);
    ++metadata_.num_edges;
  }

  void Finalize(MPIEnviroment &mpienv) {
    if (is_opened_) {
      for (auto &file : files_) {
        file->close();
      }

      if (mpienv.rank == 0)
      {
        std::ofstream info_file(file_prefix_ + ".edges.info");
        metadata_.Serialize(info_file, mpienv);
        info_file.close();
      }
      is_opened_ = false;
    }
  }

  void Finalize() {
    if (is_opened_) {
      for (auto &file : files_) {
        file->close();
      }
      is_opened_ = false;
    }
  }
};

template <class KmerType>
class MPIEdgeWriter {
 private:
  std::mutex buffer_mutex;
  std::string file_prefix_;
  std::vector<std::unique_ptr<std::ofstream>> files_;
  std::vector<int64_t> n_edges_at_thread_;
  unsigned k_;
  unsigned last_shift_;
  unsigned num_threads_;
  std::vector<uint32_t> buffer_;
  MPIEnviroment mpienv_;
  MPI_File mpi_file;                       // MPI 文件句柄
  static const unsigned BUFSIZE = 268435456; // 256 MB
  bool is_opened_;
  
  public:
  EdgeIoMetadata metadata_{};
  unsigned words_per_kmer_;

  MPIEdgeWriter(unsigned k, const std::string &out_prefix, MPIEnviroment &mpienv)
   : is_opened_(false), k_(k), mpienv_(mpienv) {
    last_shift_ = k_ % 16;
    last_shift_ = (last_shift_ == 0 ? 0 : 16 - last_shift_) * 2;
    words_per_kmer_ = DivCeiling(k_ * 2 + kBitsPerMul, 32);
    buffer_.reserve(BUFSIZE);

    SetFilePrefix(out_prefix);
    SetUnordered();
    SetKmerSize(k_ - 1);
    InitFiles(mpienv_);

    int rc = MPI_File_open(
    MPI_COMM_WORLD,             // 通信域（单进程可用 MPI_COMM_SELF）
    (file_prefix_ + ".rank.0.edges.0").c_str(),                 // 文件名
    MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
    MPI_INFO_NULL,
    &mpi_file
    );

    if (rc != MPI_SUCCESS) {
        std::cerr << "Error: Failed to open MPI file!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
   };
  ~MPIEdgeWriter() { Finalize(); }

  void SetKmerSize(uint32_t k) {
    metadata_.kmer_size = k;
    metadata_.words_per_edge = DivCeiling((k + 1) * 2 + 16, 32); // (k + 1)-mer and counts(16bits)
  }

  void SetNumThreads(int32_t num_threads) { metadata_.num_files = num_threads; }
  void SetFilePrefix(const std::string &prefix) { file_prefix_ = prefix; }
  void SetNumBuckets(int num_buckets) {
    metadata_.buckets.clear();
    //metadata_.buckets_reduce.clear();
    metadata_.buckets.resize(num_buckets);
    //metadata_.buckets_reduce.resize(num_buckets);
  }
  void SetUnordered() {
    metadata_.buckets.clear();
    metadata_.is_sorted = false;
    metadata_.num_files = 1;
    metadata_.num_edges = 0;
  }

  void InitFiles(MPIEnviroment& mpienv) {
    assert(!is_opened_);
    n_edges_at_thread_.resize(metadata_.num_files, 0);

    // for (unsigned i = 0; i < metadata_.num_files; ++i) {
    //   // files_.emplace_back(new std::ofstream(
    //   //     (file_prefix_ + ".rank." + std::to_string(mpienv.rank) + ".edges." + std::to_string(i)).c_str(),
    //   //     std::ofstream::binary | std::ofstream::out));
    // }

    is_opened_ = true;
  }

  void WriteUnordered(uint32_t *edge_ptr) {
    assert(!metadata_.is_sorted);
    files_[0]->write(reinterpret_cast<const char *>(edge_ptr),
                     sizeof(uint32_t) * metadata_.words_per_edge);
    ++metadata_.num_edges;
  }

  void Finalize(MPIEnviroment &mpienv) {
    if (is_opened_) {
      for (auto &file : files_) {
        file->close();
      }

      if (mpienv.rank == 0)
      {
        std::ofstream info_file(file_prefix_ + ".edges.info");
        metadata_.UnorderSerialize(info_file, mpienv);
        info_file.close();
      }
      is_opened_ = false;
    }
    MPI_File_close(&mpi_file);
  }

  void Finalize() {
    if (is_opened_) {
      for (auto &file : files_) {
        file->close();
      }
      is_opened_ = false;
    }
  }

  void WriteToBuf(const KmerType &kmer, mul_t mul) {
    //std::fill(buffer.begin(), buffer.end(), 0);
    std::vector<uint32_t> buffer(words_per_kmer_, 0);
    auto ptr = buffer.begin();
    uint32_t w = 0;

    for (unsigned j = 0; j < k_; ++j) {
      w = (w << 2) | kmer.GetBase(k_ - 1 - j);
      if (j % 16 == 15) {
        *ptr = w;
        w = 0;
        ++ptr;
      }
    }

    assert(ptr - buffer.begin() < words_per_kmer_);
    *ptr = (w << last_shift_);
    assert((buffer.back() & kMaxMul) == 0);
    buffer.back() |= mul;

    {
      std::lock_guard<std::mutex> lock(buffer_mutex);
      buffer_.insert(buffer_.end(), buffer.begin(), buffer.end());
      ++metadata_.num_edges;
    }
    //writer_.WriteUnordered(buffer.data());
  }

  bool check_buf() {
    return buffer_.size() >= 268435456;
  }

  void MPIFileWrite() {
    {
      std::lock_guard<std::mutex> lock(buffer_mutex);
      if (buffer_.empty()) {
        return;
      }
      // 原子写入（自动追加到文件末尾）
      MPI_File_write_shared(
        mpi_file,
        buffer_.data(),
        buffer_.size(),
        MPI_INT32_T,
        MPI_STATUS_IGNORE
      );
      buffer_.clear();  // 清空缓冲区
    }
  }

  void allreduce() {
    MPI_Allreduce(MPI_IN_PLACE, &metadata_.num_edges, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE, &n_bases_, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  }
};

#endif  // MEGAHIT_EDGE_WRITER_H