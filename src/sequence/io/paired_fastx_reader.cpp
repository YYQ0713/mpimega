//
// Created by vout on 5/11/19.
//

#include "paired_fastx_reader.h"

int64_t PairedFastxReader::Read(SeqPackage *pkg, int64_t max_num,
                                int64_t max_num_bases, bool reverse) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; i += 2) {
    auto r0 = readers_[0]->ReadNext();
    auto r1 = readers_[1]->ReadNext();
    if (r0 && r1) {
      int b0 = 0, e0 = r0->seq.l;
      int b1 = 0, e1 = r1->seq.l;

      // if (trim_qc_) {
      //   FastxReader::TrimQC(r0->qual.s, r0->qual.l, &e0);
      //   FastxReader::TrimQC(r1->qual.s, r1->qual.l, &e1);
      // }

      if (trim_n_) {
        FastxReader::TrimN(r0->seq.s, e0, &b0, &e0);
        FastxReader::TrimN(r1->seq.s, e1, &b1, &e1);
      }

      if (reverse) {
        pkg->AppendReversedStringSequence(r0->seq.s + b0, e0 - b0);
        pkg->AppendReversedStringSequence(r1->seq.s + b1, e1 - b1);
      } else {
        pkg->AppendStringSequence(r0->seq.s + b0, e0 - b0);
        pkg->AppendStringSequence(r1->seq.s + b1, e1 - b1);
      }

      num_bases += e0 - b0 + e1 - b1;

      if (num_bases >= max_num_bases) {
        return i + 2;
      }
    } else {
      return i;
    }
  }
  return max_num;
}

// int64_t PairedFastxReader::Read(SeqPackage *pkg, int64_t max_num,
//                                 int64_t max_num_bases, bool reverse) {
//   int64_t num_bases = 0;
//   int64_t num_reads = 0;
//   const int BATCH_SIZE = 1u << 20u;
//   std::vector<std::string> batch_seqs0;
//   std::vector<std::string> batch_seqs1;
//   std::vector<std::string> batch_quals0;
//   std::vector<std::string> batch_quals1;
//   batch_seqs0.reserve(BATCH_SIZE);
//   batch_seqs1.reserve(BATCH_SIZE);
//   batch_quals0.reserve(BATCH_SIZE);
//   batch_quals1.reserve(BATCH_SIZE);
//   std::vector<int> trim_begins0, trim_begins1;
//   std::vector<int> trim_ends0, trim_ends1;
//   trim_begins0.resize(BATCH_SIZE);
//   trim_begins1.resize(BATCH_SIZE);
//   trim_ends0.resize(BATCH_SIZE);
//   trim_ends1.resize(BATCH_SIZE);

//   for (int64_t i = 0; i < max_num; i += BATCH_SIZE) {
//     batch_seqs0.clear(); batch_seqs1.clear();
//     batch_quals0.clear(); batch_quals1.clear();

//     for (int i = 0; i < BATCH_SIZE; ++i) {
//       auto r0 = readers_[0]->ReadNext();
//       auto r1 = readers_[1]->ReadNext();
//       if (r0 && r1) {
//         // 这里构造string，会拷贝内容一次
//         batch_seqs0.emplace_back(r0->seq.s, r0->seq.l);
//         batch_seqs1.emplace_back(r1->seq.s, r1->seq.l);
//         batch_quals0.emplace_back(r0->qual.s, r0->qual.l);
//         batch_quals1.emplace_back(r1->qual.s, r1->qual.l);
//       }
//     }

//     #pragma omp parallel for
//     for (size_t j = 0; j < batch_quals0.size(); j++) {
//       if (trim_qc_) {
//         FastxReader::TrimQC(batch_quals0[j].c_str(), batch_quals0[j].size(), &trim_ends0[j]);
//         FastxReader::TrimQC(batch_quals1[j].c_str(), batch_quals1[j].size(), &trim_ends1[j]);
//       }
//     }

//     #pragma omp parallel for
//     for (size_t j = 0; j < batch_quals0.size(); j++) {
//       if (trim_n_) {
//         FastxReader::TrimN(batch_quals0[j].c_str(), trim_ends0[j], &trim_begins0[j], &trim_ends0[j]);
//         FastxReader::TrimN(batch_quals1[j].c_str(), trim_ends1[j], &trim_begins1[j], &trim_ends1[j]);
//       }
//     }
    
//     if (reverse) {
//       for (size_t k = 0; k < batch_seqs0.size(); k++) {
//         auto l0 = trim_ends0[k] - trim_begins0[k];
//         auto l1 = trim_ends1[k] - trim_begins1[k];
//         pkg->AppendReversedStringSequence(batch_seqs0[k].c_str() + trim_begins0[k], l0);
//         pkg->AppendReversedStringSequence(batch_seqs1[k].c_str() + trim_begins1[k], l1);
//         num_reads += 2;
//         num_bases += l0 + l1;
//       }
//     } else {
//       for (size_t k = 0; k < batch_seqs0.size(); k++) {
//         auto l0 = trim_ends0[k] - trim_begins0[k];
//         auto l1 = trim_ends1[k] - trim_begins1[k];
//         pkg->AppendStringSequence(batch_seqs0[k].c_str() + trim_begins0[k], l0);
//         pkg->AppendStringSequence(batch_seqs1[k].c_str() + trim_begins1[k], l1);
//         num_reads += 2;
//         num_bases += l0 + l1;
//       }
//     }

//     if (num_bases >= max_num_bases) {
//       return num_reads;
//     }
    
//   }

//   return num_reads;
// }