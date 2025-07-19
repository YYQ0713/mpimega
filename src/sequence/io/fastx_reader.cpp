//
// Created by vout on 4/28/19.
//

#include "fastx_reader.h"
#include <cassert>
#include <stdexcept>
#include <emmintrin.h>
#include <cstdint>

FastxReader::FastxReader(const std::string &file_name) {
  fp_ = file_name == "-" ? gzdopen(fileno(stdin), "r")
                         : gzopen(file_name.c_str(), "r");
  if (fp_ == nullptr) {
    throw std::invalid_argument("Cannot open file " + file_name);
  }
  kseq_reader_ = kseq_init(fp_);
  assert(kseq_reader_ != nullptr);
}

FastxReader::~FastxReader() {
  if (kseq_reader_) {
    kseq_destroy(kseq_reader_);
  }
  if (fp_) {
    gzclose(fp_);
  }
}

int64_t FastxReader::Read(SeqPackage *pkg, int64_t max_num,
                          int64_t max_num_bases, bool reverse) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; ++i) {
    auto record = ReadNext();
    if (record) {
      int b = 0, e = record->seq.l;
      if (trim_n_) {
        TrimN(record->seq.s, record->seq.l, &b, &e);
      }

      if (reverse) {
        pkg->AppendReversedStringSequence(record->seq.s + b, e - b);
      } else {
        pkg->AppendStringSequence(record->seq.s + b, e - b);
      }

      num_bases += e - b;
      if (num_bases >= max_num_bases && i % 2 == 1) {
        return i + 1;
      }
    } else {
      return i;
    }
  }
  return max_num;
}

void FastxReader::TrimN(const char *s, int len, int *out_bpos, int *out_epos) {
  *out_bpos = *out_epos = len;
  int i;
  for (i = 0; i < len; ++i) {
    if (s[i] == 'N' || s[i] == 'n') {
      if (*out_bpos < len) {
        break;
      }
    } else {
      if (*out_bpos == len) {
        *out_bpos = i;
      }
    }
  }
  *out_epos = i;
}

// bool FastxReader::TrimQC(const char *s, int len, int *out_epos, int window = 16, int threshold = 20) {
//     const int Q_THRESHOLD_ASCII = 33 + 15;  // ASCII < 48 → Q < 15
//     const int max_low_q = static_cast<int>(len * 0.2);
//     int low_q_count = 0;

//     int i = len;
//     const __m128i threshold = _mm_set1_epi8(Q_THRESHOLD_ASCII);

//     // Process 16 characters at a time, from high index to low
//     while (i >= 16) {
//         i -= 16;
//         __m128i data = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s + i));
//         __m128i cmp = _mm_cmplt_epi8(data, threshold);  // compare < 48
//         int mask = _mm_movemask_epi8(cmp);
//         low_q_count += __builtin_popcount(mask);
//         if (low_q_count > max_low_q)
//             return false;
//     }

//     return true;
// }

bool FastxReader::TrimQC(const char *s, int len, int *out_epos, int window, int threshold) {
    const int Q_BASE = 33;
    const int AVG_Q_CUTOFF = threshold * window;  // e.g., Q15×4 = 60
    const int MAX_CHECK_LEN = 30;
    *out_epos = len;  // 默认不截断
    
    if (len < window) return false;

    int start = std::max(0, len - MAX_CHECK_LEN);
    int end = len;

    // 初始化第一个窗口的总和
    int sum = 0;
    for (int i = end - window; i < end; ++i) {
        sum += static_cast<unsigned char>(s[i]) - Q_BASE;
    }
    if (sum < AVG_Q_CUTOFF) {
        *out_epos = end - window;
        return true;
    }

    // 向前滑动窗口（从末尾向前）
    for (int i = end - 1; i >= start + window; --i) {
        int out_val = static_cast<unsigned char>(s[i]) - Q_BASE;
        int in_val  = static_cast<unsigned char>(s[i - window]) - Q_BASE;
        sum = sum - out_val + in_val;

        if (sum < AVG_Q_CUTOFF) {
            *out_epos = i - window;
            return true;
        }
    }

    return false;
}
