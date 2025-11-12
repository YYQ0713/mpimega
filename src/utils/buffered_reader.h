//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_BUFFERED_READER_H
#define MEGAHIT_BUFFERED_READER_H

#include <cstring>
#include <fstream>
#include "mpienv/mpienv.hpp"

/**
 * A buffered wrapper to speed up ifstream::read
 */
class BufferedReader {
 public:
  static constexpr size_t kMaxBufferSize = 65536;
  static constexpr size_t kMaxBigBufferSize = 256 * 1024 * 1024;
  explicit BufferedReader() = default;
  void reset(std::istream *is, size_t buffer_size = kMaxBufferSize) {
    is_ = is;
    head_ = tail_ = 0;
    buffer_size_ = std::min(buffer_size, kMaxBufferSize * 1);
  }

  void reset_mpi(std::istream *is, MPI_File mpi_file, size_t buffer_size = kMaxBufferSize) {
    is_ = is;
    head_ = tail_ = 0;
    mpi_file_ = mpi_file;
    buffer_size_ = std::min(buffer_size, kMaxBufferSize * 1);
    big_buffer_.resize(kMaxBigBufferSize);
  }

  template <typename T>
  size_t read(T *dst, size_t size = 1) {
    if (is_ == nullptr) {
      return 0;
    }
    size_t wanted = sizeof(T) * size;
    size_t remained = wanted;
    auto dst_ptr = reinterpret_cast<char *>(dst);
    while (remained > 0) {
      if (remained <= tail_ - head_) {
        memcpy(dst_ptr, buffer_ + head_, remained);
        head_ += remained;
        return wanted;
      } else {
        if (tail_ > head_) {
          memcpy(dst_ptr, buffer_ + head_, tail_ - head_);
          remained -= tail_ - head_;
          dst_ptr += tail_ - head_;
        }
        if (refill() == 0) {
          return wanted - remained;
        }
      }
    }
    return 0;
  }

  template <typename T>
  size_t read_mpi(T *dst, size_t size = 1) {
    if (is_ == nullptr) {
      return 0;
    }
    size_t wanted = sizeof(T) * size;
    size_t remained = wanted;
    auto dst_ptr = reinterpret_cast<char *>(dst);
    while (remained > 0) {
      if (remained <= tail_ - head_) {
        memcpy(dst_ptr, big_buffer_.data() + head_, remained);
        head_ += remained;
        return wanted;
      } else {
        if (tail_ > head_) {
          memcpy(dst_ptr, big_buffer_.data() + head_, tail_ - head_);
          remained -= tail_ - head_;
          dst_ptr += tail_ - head_;
        }
        if (refill_mpi() == 0) {
          return wanted - remained;
        }
      }
    }
    return 0;
  }

 private:
  size_t refill() {
    head_ = 0;
    is_->read(buffer_, buffer_size_);
    tail_ = *is_ ? buffer_size_ : is_->gcount();
    return tail_;
  }

  size_t refill_mpi() {
    head_ = 0;
    MPI_Status status;
    MPI_File_read_all(mpi_file_, big_buffer_.data(), big_buffer_size_, MPI_BYTE, &status);
    int count;
    MPI_Get_count(&status, MPI_BYTE, &count);
    tail_ = (count == big_buffer_size_) ? big_buffer_size_ : static_cast<size_t>(count);
    
    return tail_;
  }
 private:
  std::istream *is_{};
  MPI_File mpi_file_;
  char buffer_[kMaxBufferSize]{};
  std::vector<char> big_buffer_;
  size_t buffer_size_{kMaxBufferSize};
  size_t big_buffer_size_{kMaxBigBufferSize};
  size_t head_{0};
  size_t tail_{0};
};

#endif  // MEGAHIT_BUFFERED_READER_H
