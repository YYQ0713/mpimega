
#ifndef PRL_BLOOM_FILTER
#define PRL_BLOOM_FILTER


#include "kmlib/kmbitvector.h"
#include "xxhash/xxh3.h"
#include <cmath>
#include <iostream>

class Bloom {
  public:
    Bloom(uint64_t elements, double error) : elements_(elements), error_(error) { bloom_init(); }
    
    template <unsigned NumWords, class WordType>
    bool bloom_check_add(const Kmer<NumWords, WordType> &kmer) {
      int hits = 0;
      unsigned int a1 = XXH3_64bits_withSeed(static_cast<const void *>(kmer.data()), sizeof(WordType) * NumWords, 0x9747b28c);
      unsigned int a2 = XXH3_64bits_withSeed(static_cast<const void *>(kmer.data()), sizeof(WordType) * NumWords, a1);
      unsigned int b1 = XXH3_64bits_withSeed(static_cast<const void *>(kmer.data()), sizeof(WordType) * NumWords, a2);
      unsigned int b2 = XXH3_64bits_withSeed(static_cast<const void *>(kmer.data()), sizeof(WordType) * NumWords, b1);
      uint64_t a = (((uint64_t) a1) << 32) | ((uint64_t) a2);
      uint64_t b = (((uint64_t) b1) << 32) | ((uint64_t) b2);
      uint64_t x;
      uint64_t byte;
      unsigned int mask;
      unsigned int i;
      unsigned char c;

      for (i = 0; i < hashes_; i++) {
        x = (a + i * b) % bits_;
        hits += !bloom_.try_lock(x);
      }

      if (hits == hashes_) {
        return true;                   // 1 == element already in (or collision)
      }

      return false;
    }

  private:
    int bloom_init() {
      if (elements_ < 1 || error_ <= 0.0 || error_ >= 1.0) {
        return 1;
      }

      double num = std::log(error_);
      double denom = 0.480453013918201; // ln(2)^2
      bpe_ = -(num / denom);

      double dentries = (double)elements_;
      bits_ = (int64_t)(dentries * bpe_);

      hashes_ = (int)std::ceil(0.693147180559945 * bpe_);  // ln(2)
      bloom_.reset(bits_);

      return 0;
    }

  private:
    AtomicBitVector bloom_;
    uint64_t elements_;
    uint64_t bits_;
    double bpe_;
    double error_;
    int hashes_;

};


#endif