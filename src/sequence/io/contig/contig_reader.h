//
// Created by vout on 4/28/19.
//

#ifndef MEGAHIT_CONTIG_READER_H
#define MEGAHIT_CONTIG_READER_H

#include <fstream>
#include "sequence/io/fastx_reader.h"

class ContigReader : public FastxReader {
 public:
  explicit ContigReader(const std::string &file_name)
      : FastxReader(file_name), file_name_(file_name) {}
  ContigReader *SetMinLen(unsigned min_len) {
    min_len_ = min_len;
    return this;
  }
  ContigReader *SetExtendLoop(unsigned k_from, unsigned k_to) {
    k_from_ = k_from;
    k_to_ = k_to;
    return this;
  }
  ContigReader *SetDiscardFlag(unsigned flag) {
    discard_flag_ = flag;
    return this;
  }

  std::pair<int64_t, int64_t> GetNumContigsAndBases() const {
    std::ifstream info_fs(file_name_ + ".info");
    int64_t num_contigs, num_bases;
    info_fs >> num_contigs >> num_bases;
    if (!info_fs) {
      xfatal("Invalid format of contig info file: {s}.info",
             file_name_.c_str());
    }
    return {num_contigs, num_bases};
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse) override {
    return ReadWithMultiplicity<float>(pkg, nullptr, max_num, max_num_bases,
                                       reverse);
  }

  int64_t ReadDistbt(SeqPackage *pkg, uint64_t ctg_sz, MPIEnviroment &mpienv, bool reverse) {
    // 每个进程负责的区间计算
    uint64_t base = ctg_sz / mpienv.nprocs;
    uint64_t rem  = ctg_sz % mpienv.nprocs;
    uint64_t start_ctg = mpienv.rank < rem ? (mpienv.rank * (base + 1)) : (rem * (base + 1) + (mpienv.rank - rem) * base);
    uint64_t end_ctg   = start_ctg + (mpienv.rank < rem ? (base + 1) : base);
    
    int64_t max_num = 1LL << 60;
    int64_t max_num_bases = 1LL << 60;

    // 先跳过自己的起始 contig 之前的部分
    uint64_t cur_ctg_idx = 0;
    while (cur_ctg_idx < start_ctg) {
        auto record = ReadNext();
        if (!record) return 0; // 文件提前结束
        ++cur_ctg_idx;
    }

    // 原本的逻辑保留，只是加上区间结束判断
    bool extend_loop = k_from_ < k_to_ && !(discard_flag_ & contig_flag::kLoop);
    int64_t num_bases = 0;
    int64_t ri;

    for (ri = 0; cur_ctg_idx < end_ctg; ++ri, ++cur_ctg_idx) {
      auto record = ReadNext();
      if (!record) {
          return ri;
      }

      if (record->seq.l < min_len_) {
          --ri;
          continue;
      }

      unsigned flag = record->comment.s[5] - '0';
      if (discard_flag_ & flag) {
          --ri;
          continue;
      }

      if (extend_loop && (flag & contig_flag::kLoop)) {
          if (record->seq.l < k_to_ + 1U) {
              continue;
          }
          std::string ss(record->seq.s);
          for (unsigned i = k_from_; i < k_to_; ++i) {
              ss.push_back(ss[i]);
          }
          if (reverse) {
              pkg->AppendReversedStringSequence(ss.c_str(), ss.length());
          } else {
              pkg->AppendStringSequence(ss.c_str(), ss.length());
          }
      } else {
          if (reverse) {
              pkg->AppendReversedStringSequence(record->seq.s, record->seq.l);
          } else {
              pkg->AppendStringSequence(record->seq.s, record->seq.l);
          }
      }

      num_bases += record->seq.l;
      if (num_bases >= max_num_bases) {
          return ri + 1;
      }
    }

    return ri;
  }

  template <typename TMul>
  int64_t ReadAllWithMultiplicity(SeqPackage *pkg, std::vector<TMul> *mul,
                                  bool reverse) {
    return ReadWithMultiplicity(pkg, mul, kMaxNumSeq, kMaxNumBases, reverse);
  }

  template <typename TMul>
  int64_t ReadWithMultiplicity(SeqPackage *pkg, std::vector<TMul> *mul,
                               int64_t max_num, int64_t max_num_bases,
                               bool reverse) {
    bool extend_loop = k_from_ < k_to_ && !(discard_flag_ & contig_flag::kLoop);

    int64_t num_bases = 0;
    for (int64_t ri = 0; ri < max_num; ++ri) {
      auto record = ReadNext();
      if (record) {
        if (record->seq.l < min_len_) {
          --ri;
          continue;
        }
        // comment = "flag=x multi=xx.xxxx"
        unsigned flag = record->comment.s[5] - '0';
        if (discard_flag_ & flag) {
          --ri;
          continue;
        }

        if (extend_loop && (flag & contig_flag::kLoop)) {
          if (record->seq.l < k_to_ + 1U) {
            continue;
          }
          std::string ss(record->seq.s);
          for (unsigned i = k_from_; i < k_to_; ++i) {
            ss.push_back(ss[i]);
          }

          if (reverse) {
            pkg->AppendReversedStringSequence(ss.c_str(), ss.length());
          } else {
            pkg->AppendStringSequence(ss.c_str(), ss.length());
          }
        } else {
          if (reverse) {
            pkg->AppendReversedStringSequence(record->seq.s, record->seq.l);
          } else {
            pkg->AppendStringSequence(record->seq.s, record->seq.l);
          }
        }

        if (mul) {
          mul->push_back(GetMultiplicity<TMul>(record->comment.s));
        }

        num_bases += record->seq.l;
        if (num_bases >= max_num_bases) {
          return ri + 1;
        }
      } else {
        return ri;
      }
    }
    return max_num;
  }

 private:
  template <typename TMul>
  static TMul GetMultiplicity(const char *fastx_comment) {
    auto m = atof(fastx_comment + 13);
    if (std::is_integral<TMul>::value) {
      return m + .5;
    } else {
      return m;
    }
  }

 private:
  unsigned min_len_{0};
  unsigned k_from_{0}, k_to_{0};
  unsigned discard_flag_{0};
  std::string file_name_;
};

#endif  // MEGAHIT_CONTIG_READER_H
