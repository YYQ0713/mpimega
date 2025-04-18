/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#ifndef MEGAHIT_SEQUENCE_PACKAGE_GA_H
#define MEGAHIT_SEQUENCE_PACKAGE_GA_H

#include <cassert>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <memory>
#include <vector>
#include "kmlib/kmbit.h"
#include "kmlib/kmcompactvector.h"
#include "utils/mutex.h"
#include "utils/utils.h"
#include "ga.h"

/**
 * @brief hold a set of sequences
 */
template <class WordType = unsigned long>
class SequencePackageGA {
 public:
  /**
   * The sequence view of a sequence
   */
  class SeqView {
   public:
    SeqView(const SequencePackageGA *pkg, size_t seq_id)
        : package_(pkg), seq_id_(seq_id) {
      if (seq_id >= pkg->seq_count())
      {
        xinfo("seqid: {} pkg->seqcount: {}\n", seq_id, pkg->seq_count());
      }
      assert(seq_id < pkg->seq_count());
    }

    unsigned length() const { return package_->GetSeqLength(seq_id_); }

    std::pair<const WordType *, unsigned> raw_address(
        unsigned offset = 0) const {
      assert(offset < length());
      return package_->GetRawAddress(seq_id_, offset);
    }

    void raw_address_ga(std::pair<std::vector<uint32_t>, unsigned> &rhs,
        unsigned offset = 0) const {
      assert(offset < length());
      return package_->GetRawAddressGA(seq_id_, rhs, length(), offset);
    }

    size_t base_at(unsigned index) const {
      assert(index < length());
      return package_->GetBase(seq_id_, index);
    }

    size_t base_at_ga(std::pair<std::vector<uint32_t>, unsigned>& start_ptr_and_offset, unsigned index) const {
      assert(index < length());
      return package_->GetBaseGA(start_ptr_and_offset, seq_id_, index);
    }

    size_t id() const { return seq_id_; }

    size_t full_offset_in_pkg() const { return package_->StartPos(seq_id_); }

   private:
    const SequencePackageGA *package_;
    size_t seq_id_;
  };

  using TWord = WordType;
  using TVector = kmlib::CompactVector<2, TWord, kmlib::kBigEndian>;
  using TAddress = std::pair<const TWord *, unsigned>;
  using TAddressGA = std::pair<std::vector<uint32_t>, unsigned>;
  const static unsigned kBasesPerWord = TVector::kBasesPerWord;

 public:
  SequencePackageGA() {
    Clear();
    for (int i = 0; i < 10; ++i) {
      dna_map_[static_cast<int>("ACGTNacgtn"[i])] = "0123201232"[i] - '0';
    }
  }

  void Clear() {
    sequences_.clear();
    start_pos_.clear();
    start_pos_.push_back(0);
    pos_to_id_.clear();
    max_len_ = 0;
    fixed_len_ = 0;
    num_fixed_len_ = 0;
  }

  void ReserveBases(size_t num_bases) { sequences_.reserve(num_bases); }

  void ReserveSequences(size_t num_seq) { start_pos_.reserve(num_seq + 1); }

  void SeqGaReserve(size_t num_bases, size_t num_seq) {
    int64_t ndim, dims[2], chunk[2];
    int64_t lo[2], hi[2];

    // seq_pkg_ga_ create set
    ndim = 2;
    dims[0] = 1;
    dims[1] = DivCeiling(num_bases, kBasesPerWord);
    chunk[0] = -1;
    chunk[1] = -1;
    seq_pkg_ga_ld_ = dims[1];
    seq_pkg_ga_ = NGA_Create64(C_INT, ndim, dims, "seq_pkg_ga_", chunk);
    assert(seq_pkg_ga_);
    NGA_Distribution64(seq_pkg_ga_, 0, lo, hi);
    local_size_ga = hi[1] - lo[1] + 1;
  
  }

  size_t seq_count() const { return num_fixed_len_ + start_pos_.size() - 1; }

  size_t base_count() const { return start_pos_.back(); }

  size_t size_in_byte() const {
    return sizeof(TWord) * local_size_ga +
           sizeof(uint64_t) * start_pos_.capacity() +
           sizeof(uint64_t) * pos_to_id_.capacity();
  }

  size_t size_in_word_of_sequences() const {return sequences_.word_count(); }

  unsigned max_length() const { return max_len_; }

  SeqView GetSeqView(size_t seq_id) const { return SeqView(this, seq_id); }

  SeqView GetSeqViewByOffset(size_t offset) const {
    return SeqView(this, GetSeqID(offset));
  }

 private:
  unsigned GetSeqLength(size_t seq_id) const {
    if (seq_id < num_fixed_len_) {
      return fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_ + 1] -
             start_pos_[seq_id - num_fixed_len_];
    }
  }

  uint8_t GetBase(size_t seq_id, unsigned offset) const {
    return sequences_[StartPos(seq_id) + offset];
  }

  //修改适用于GA的GetBase函数
  uint8_t GetBaseGA(std::pair<std::vector<uint32_t>, unsigned>& start_ptr_and_offset, size_t seq_id, unsigned offset) const {
    size_t totoffset = start_ptr_and_offset.second + offset;
    size_t which_word = totoffset / 16;
    size_t locoffset = totoffset % 16;
    return (start_ptr_and_offset.first[which_word] >> ((15 - locoffset) * 2)) & 0x03;
  }

  uint64_t StartPos(size_t seq_id) const {
    if (seq_id < num_fixed_len_) {
      return seq_id * fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_];
    }
  }

  TAddress GetRawAddress(size_t seq_id, unsigned offset = 0) const {
    size_t index = StartPos(seq_id) + offset;
    return {sequences_.data() + index / kBasesPerWord, index % kBasesPerWord};
  }

  void GetRawAddressGA(size_t seq_id, std::pair<std::vector<uint32_t>, unsigned> &rhs, size_t length, unsigned offset = 0) const {
    size_t index = StartPos(seq_id) + offset;
    size_t which_word = index / kBasesPerWord;
    size_t rem_word = index % kBasesPerWord;
    size_t seq_len = ((length + kBasesPerWord - 1) / kBasesPerWord);
    //rhs.first = (uint32_t*)malloc((seq_len + 1) * sizeof(uint32_t));
    rhs.first = std::vector<uint32_t>(seq_len + 1);
    int64_t lo[2], hi[2], ld;
    lo[0] = hi[0] = 0;
    lo[1] = which_word;
    hi[1] = (which_word + seq_len) < seq_pkg_ga_ld_ ? (which_word + seq_len) : (seq_pkg_ga_ld_ - 1);
    ld = seq_pkg_ga_ld_;
    //assert(hi[1] < ld);
    #pragma omp critical
    {
      NGA_Get64(seq_pkg_ga_, lo, hi, rhs.first.data(), &ld);
    }
    rhs.second = rem_word;
  }

  uint64_t GetSeqID(size_t full_offset) const {
    if (full_offset < num_fixed_len_ * fixed_len_) {
      return full_offset / fixed_len_;
    } else {
      size_t look_up_entry =
          (full_offset - num_fixed_len_ * fixed_len_) / kLookupStep;
      size_t l = pos_to_id_[look_up_entry], r = pos_to_id_[look_up_entry + 1];

      while (l < r) {
        size_t mid = (l + r) / 2;
        if (start_pos_[mid - num_fixed_len_] > full_offset) {
          r = mid - 1;
        } else if (start_pos_[mid - num_fixed_len_ + 1] <= full_offset) {
          l = mid + 1;
        } else {
          return mid;
        }
      }
      return l;
    }
  }

 public:
  void AppendStringSequence(const char *s, unsigned len) {
    AppendStringSequence(s, s + len, len);
  }

  void AppendReversedStringSequence(const char *s, unsigned len) {
    AppendStringSequence(s + len - 1, s - 1, len);
  }

  void AppendCompactSequence(const TWord *s, unsigned len) {
    AppendCompactSequence(s, len, false);
  }

  void AppendReversedCompactSequence(const TWord *s, unsigned len) {
    AppendCompactSequence(s, len, true);
  }

  void AppendCompactSequenceV(const TWord *s, unsigned len) {
    AppendCompactSequenceV(s, len, false);
  }

  void AppendReversedCompactSequenceV(const TWord *s, unsigned len) {
    AppendCompactSequenceV(s, len, true);
  }

  void FetchSequence(size_t seq_id, std::vector<TWord> *s) const {
    TVector cvec(s);
    auto ptr_and_offset = GetRawAddress(seq_id);
    auto ptr = ptr_and_offset.first;
    auto offset = ptr_and_offset.second;
    auto len = GetSeqLength(seq_id);
    if (offset != 0) {
      unsigned remaining_len = std::min(len, kBasesPerWord - offset);
      cvec.push_word(*ptr, offset, remaining_len);
      len -= remaining_len;
      ++ptr;
    }
    unsigned n_full = len / kBasesPerWord;
    for (unsigned i = 0; i < n_full; ++i) {
      cvec.push_word(ptr[i]);
    }
    if (len % kBasesPerWord > 0) {
      cvec.push_word(ptr[n_full], 0, len % kBasesPerWord);
    }
  }

  void BuildIndex() {
    pos_to_id_.clear();
    pos_to_id_.reserve(start_pos_.back() / kLookupStep + 4);
    size_t abs_offset = num_fixed_len_ * fixed_len_;
    size_t cur_id = num_fixed_len_;

    while (abs_offset <= start_pos_.back()) {
      while (cur_id < seq_count() &&
             start_pos_[cur_id - num_fixed_len_ + 1] <= abs_offset) {
        ++cur_id;
      }

      pos_to_id_.push_back(cur_id);
      abs_offset += kLookupStep;
    }

    pos_to_id_.push_back(seq_count());
    pos_to_id_.push_back(seq_count());
  }

  void WriteSequences(std::ostream &os, int64_t from = 0,
                      int64_t to = -1) const {
    if (to == -1) {
      to = seq_count() - 1;
    }

    uint32_t len;
    std::vector<uint32_t> s;

    for (int64_t i = from; i <= to; ++i) {
      len = GetSeqView(i).length();
      FetchSequence(i, &s);
      os.write(reinterpret_cast<const char *>(&len), sizeof(uint32_t));
      os.write(reinterpret_cast<const char *>(s.data()),
               sizeof(uint32_t) * s.size());
    }
  }

 private:
  bool IsFixedLength() const { return start_pos_.size() == 1; }

  void UpdateLength(unsigned len) {
    if (num_fixed_len_ == 0) {
      num_fixed_len_ = 1;
      fixed_len_ = len;
      start_pos_.back() += len;
    } else if (IsFixedLength() && len == fixed_len_) {
      num_fixed_len_++;
      start_pos_.back() += len;
    } else {
      start_pos_.push_back(start_pos_.back() + len);
    }
    if (len > max_len_) {
      max_len_ = len;
    }
  }

  void AppendStringSequence(const char *from, const char *to, unsigned len) {
    if (len == 0) {
      // Fake a sequence whose length is 1, as we need all sequences' length > 0
      // to make `GetSeqID` working
      auto fake_sequence = "A";
      return AppendStringSequence(fake_sequence, fake_sequence + 1, 1);
    }
    UpdateLength(len);
    std::ptrdiff_t step = from < to ? 1 : -1;
    for (auto ptr = from; ptr != to; ptr += step) {
      sequences_.push_back(dna_map_[static_cast<int>(*ptr)]);
    }
  }

  void AppendCompactSequence(const TWord *ptr, unsigned len, bool rev) {
    if (len == 0) {
      // Fake a sequence whose length is 1, as we need all sequences' length > 0
      // to make `GetSeqID` working
      TWord fake_sequence = 0;
      return AppendCompactSequence(&fake_sequence, 1, false);
    }
    UpdateLength(len);

    if (rev) {
      auto rptr = ptr + DivCeiling(len, kBasesPerWord) - 1;
      unsigned bases_in_last_word = len % kBasesPerWord;
      if (bases_in_last_word > 0) {
        auto val = kmlib::bit::Reverse<2>(*rptr);
        sequences_.push_word(val, kBasesPerWord - bases_in_last_word,
                             bases_in_last_word);
        --rptr;
      }
      for (auto p = rptr; p >= ptr; --p) {
        sequences_.push_word(kmlib::bit::Reverse<2>(*p));
      }
    } else {
      while (len >= kBasesPerWord) {
        sequences_.push_word(*ptr);
        len -= kBasesPerWord;
        ++ptr;
      }
      if (len > 0) {
        sequences_.push_word(*ptr, 0, len);
      }
    }
  }

  void AppendCompactSequenceV(const TWord *ptr, unsigned len, bool rev) {
    if (len == 0) {
      // Fake a sequence whose length is 1, as we need all sequences' length > 0
      // to make `GetSeqID` working
      TWord fake_sequence = 0;
      return AppendCompactSequenceV(&fake_sequence, 1, false);
    }
    UpdateLength(len);
  }

 private:
  unsigned fixed_len_{0};
  size_t num_fixed_len_{0};
  char dna_map_[256]{};
  unsigned max_len_{0};

  // for looking up the seq_id of a full offset
  std::vector<uint64_t> pos_to_id_;
  const static unsigned kLookupStep = 1024;

 public:
  std::vector<uint64_t>
      start_pos_;  // the index of the starting position of a sequence
  TVector sequences_;
  int seq_pkg_ga_;
  int64_t seq_pkg_ga_ld_;
  int64_t local_size_ga;
};

using SeqPackageGA = SequencePackageGA<uint32_t>;

#endif
