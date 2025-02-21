//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_BINARY_READER_H
#define MEGAHIT_BINARY_READER_H

#include <fstream>
#include "base_reader.h"

#include "utils/buffered_reader.h"

class BinaryReader : public BaseSequenceReader {
 public:
  explicit BinaryReader(const std::string &filename)
      : is_(filename), buf_(120) {
    if (is_.bad()) {
      throw std::invalid_argument("Failed to open file " + filename);
    }
    reader_.reset(&is_);
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse) override {
    max_num = (max_num + 1) / 2 * 2;
    int64_t num_bases = 0;
    uint32_t read_len;

    for (int64_t i = 0; i < max_num; ++i) {
      if (reader_.read(&read_len) == 0) {
        return i;
      }
      auto num_words = DivCeiling(read_len, SeqPackage::kBasesPerWord);
      if (buf_.size() < num_words) {
        buf_.resize(num_words);
      }
      auto bytes_read = reader_.read(buf_.data(), num_words);
      assert(bytes_read == num_words * sizeof(buf_[0]));
      (void)(bytes_read);

      if (!reverse) {
        pkg->AppendCompactSequence(buf_.data(), read_len);
      } else {
        pkg->AppendReversedCompactSequence(buf_.data(), read_len);
      }

      num_bases += read_len;
      if (read_len >= max_num_bases && i % 2 == 1) {
        return i + 1;
      }
    }
    return max_num;
  }

  int64_t ReadAllGA(SeqPackageGA *pkg, MPIEnviroment& mpienv, bool reverse) {
    return ReadGA(pkg, kMaxNumSeq, kMaxNumBases, mpienv, reverse);
  }

  int64_t ReadGA(SeqPackageGA *pkg, int64_t max_num, int64_t max_num_bases, MPIEnviroment& mpienv,
               bool reverse = false) {
    max_num = (max_num + 1) / 2 * 2;
    int64_t num_bases = 0;
    uint32_t read_len;
    int64_t ga_flag = 0;
    bool put_flag = false;
    int64_t lo[2], hi[2];
    int64_t local_size;
    NGA_Distribution64(pkg->seq_pkg_ga_, ga_flag, lo, hi);
    local_size = hi[1] - lo[1] + 1;
    pkg->sequences_.reserve_word(local_size);

    for (int64_t i = 0; i < max_num; ++i) {
      if (reader_.read(&read_len) == 0) {
        //todo
        if (mpienv.ga_rank == ga_flag && !put_flag) {
          NGA_Put64(pkg->seq_pkg_ga_, lo, hi, pkg->sequences_.data(), &pkg->seq_pkg_ga_ld_);
          pkg->sequences_.clear();
          put_flag = true;
        }
        return i;
      }
      auto num_words = DivCeiling(read_len, SeqPackage::kBasesPerWord);
      if (buf_.size() < num_words) {
        buf_.resize(num_words);
      }
      auto bytes_read = reader_.read(buf_.data(), num_words);
      assert(bytes_read == num_words * sizeof(buf_[0]));
      (void)(bytes_read);

      if (!put_flag) {
        if (!reverse) {
          pkg->AppendCompactSequence(buf_.data(), read_len);
        } else {
          pkg->AppendReversedCompactSequence(buf_.data(), read_len);
        }
      } else {
        if (!reverse) {
          pkg->AppendCompactSequenceV(buf_.data(), read_len);
        } else {
          pkg->AppendReversedCompactSequenceV(buf_.data(), read_len);
        }
      }

      num_bases += read_len;
      if (read_len >= max_num_bases && i % 2 == 1) {
        if (mpienv.ga_rank == ga_flag && !put_flag) {
          NGA_Put64(pkg->seq_pkg_ga_, lo, hi, pkg->sequences_.data(), &pkg->seq_pkg_ga_ld_);
          pkg->sequences_.clear();
          put_flag = true;
        }
        return i + 1;
      }

      if (!put_flag && pkg->size_in_word_of_sequences() > local_size && ga_flag < mpienv.ga_nprocs) {
        if (mpienv.ga_rank == ga_flag) {
          NGA_Put64(pkg->seq_pkg_ga_, lo, hi, pkg->sequences_.data(), &pkg->seq_pkg_ga_ld_);
          pkg->sequences_.clear();
          put_flag = true;
        }

        ga_flag++;
        if (!put_flag) {
          pkg->sequences_.update_ga(local_size, pkg->size_in_word_of_sequences());
          //update new range
          NGA_Distribution64(pkg->seq_pkg_ga_, ga_flag, lo, hi);
          local_size = hi[1] - lo[1] + 1;
          pkg->sequences_.reserve_word(local_size);
        }
      }
    }

    if (mpienv.ga_rank == ga_flag && !put_flag) {
      NGA_Put64(pkg->seq_pkg_ga_, lo, hi, pkg->sequences_.data(), &pkg->seq_pkg_ga_ld_);
      pkg->sequences_.clear();
      put_flag = true;
    }
    return max_num;
  }

 private:
  std::ifstream is_;
  BufferedReader reader_;
  std::vector<SeqPackage::TWord> buf_;
};

#endif  // MEGAHIT_BINARY_READER_H
