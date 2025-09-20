//
// Created by vout on 11/21/18.
//

#include "contig_output.h"
#include <cassert>
#include "definitions.h"
#include "unitig_graph.h"
#include <omp.h>

namespace {

inline char Complement(char c) {
  if (c >= 0 && c < 4) {
    return 3 - c;
  }
  switch (c) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    default:
      assert(false);
  }
  return 0;
}

inline void ReverseComplement(std::string &s) {
  int i, j;
  for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
    std::swap(s[i], s[j]);
    s[i] = Complement(s[i]);
    s[j] = Complement(s[j]);
  }
  if (i == j) {
    s[i] = Complement(s[i]);
  }
}

void FoldPalindrome(std::string &s, unsigned kmer_k, bool is_loop) {
  if (is_loop) {
    for (unsigned i = 1; i + kmer_k <= s.length(); ++i) {
      std::string rc = s.substr(i, kmer_k);
      ReverseComplement(rc);
      if (rc == s.substr(i - 1, kmer_k)) {
        assert(i <= s.length() / 2);
        s = s.substr(i, s.length() / 2);
        break;
      }
    }
  } else {
    int num_edges = s.length() - kmer_k;
    assert(num_edges % 2 == 1);
    s.resize((num_edges - 1) / 2 + kmer_k + 1);
  }
}

}  // namespace

void OutputContigs(UnitigGraph &graph, ContigWriter *contig_writer,
                   ContigWriter *final_contig_writer, bool change_only,
                   uint32_t min_standalone) {
  assert(!(change_only && final_contig_writer != nullptr));  // if output
                                                             // changed contigs,
                                                             // must not output
                                                             // final contigs

#pragma omp parallel for
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    double multi = change_only ? 1
                               : std::min(static_cast<double>(kMaxMul),
                                          adapter.GetAvgDepth());
    std::string ascii_contig = graph.VertexToDNAString(adapter);
    if (change_only && !adapter.IsChanged()) {
      continue;
    }

    if (adapter.IsLoop()) {
      int flag = contig_flag::kLoop | contig_flag::kStandalone;
      auto writer = contig_writer;

      if (adapter.IsPalindrome()) {
        FoldPalindrome(ascii_contig, graph.k(), adapter.IsLoop());
        flag = contig_flag::kStandalone;
      }

      if (final_contig_writer != nullptr) {
        if (ascii_contig.length() < min_standalone) {
          continue;
        } else {
          writer = final_contig_writer;
        }
      }
      writer->WriteContig(ascii_contig, graph.k(), i, flag, multi);
    } else {
      auto out_file = contig_writer;
      int flag = 0;

      if (adapter.IsStandalone() ||
          (graph.InDegree(adapter) == 0 && graph.OutDegree(adapter) == 0)) {
        if (adapter.IsPalindrome()) {
          FoldPalindrome(ascii_contig, graph.k(), adapter.IsLoop());
        }
        flag = contig_flag::kStandalone;
        if (final_contig_writer != nullptr) {
          if (ascii_contig.length() < min_standalone) {
            continue;
          } else {
            out_file = final_contig_writer;
          }
        }
      }
      out_file->WriteContig(ascii_contig, graph.k(), i, flag, multi);
    }
  }
}

void MPIOutputContigs(UnitigGraph &graph, MPIContigWriter *contig_writer,
                   MPIContigWriter *final_contig_writer, bool change_only,
                   uint32_t min_standalone, MPIEnviroment &mpienv) {
  assert(!(change_only && final_contig_writer != nullptr));  // if output
                                                             // changed contigs,
                                                             // must not output
                                                             // final contigs


  int64_t num_edges_mean = graph.size() / mpienv.nprocs;
  int64_t remain = graph.size() % mpienv.nprocs;
  int64_t start_index = mpienv.rank * num_edges_mean + (mpienv.rank < remain ? mpienv.rank : remain);
  int64_t end_index = start_index + num_edges_mean + (mpienv.rank < remain ? 1 : 0);
  graph.OpenReadOnly_db();
#pragma omp parallel for
  for (UnitigGraph::size_type i = start_index; i < end_index; ++i) {
    if (graph.is_del(i)) {
      continue;
    }
    auto adapter = graph.MakeVertexAdapter(i);
    double multi = change_only ? 1
                               : std::min(static_cast<double>(kMaxMul),
                                          adapter.GetAvgDepth());
    std::string ascii_contig = graph.VertexToDNAString(adapter);
    if (change_only && !adapter.IsChanged()) {
      continue;
    }

    if (adapter.IsLoop()) {
      int flag = contig_flag::kLoop | contig_flag::kStandalone;
      auto writer = contig_writer;

      if (adapter.IsPalindrome()) {
        FoldPalindrome(ascii_contig, graph.k(), adapter.IsLoop());
        flag = contig_flag::kStandalone;
      }

      if (final_contig_writer != nullptr) {
        if (ascii_contig.length() < min_standalone) {
          continue;
        } else {
          writer = final_contig_writer;
        }
      }
      writer->WriteContig(ascii_contig, graph.k(), i, flag, multi);
    } else {
      auto out_file = contig_writer;
      int flag = 0;

      if (adapter.IsStandalone() ||
          (graph.InDegree(adapter) == 0 && graph.OutDegree(adapter) == 0)) {
        if (adapter.IsPalindrome()) {
          FoldPalindrome(ascii_contig, graph.k(), adapter.IsLoop());
        }
        flag = contig_flag::kStandalone;
        if (final_contig_writer != nullptr) {
          if (ascii_contig.length() < min_standalone) {
            continue;
          } else {
            out_file = final_contig_writer;
          }
        }
      }
      out_file->WriteContig(ascii_contig, graph.k(), i, flag, multi);
    }

    // 仅主线程检查并刷新 buffer
    if (omp_get_thread_num() == 0 && contig_writer->check_buf()) {
      if (final_contig_writer != nullptr) {
        final_contig_writer->MPIFileWrite();
      } else {
        contig_writer->MPIFileWrite();
      }
    }
  }
  
  if (final_contig_writer != nullptr) {
    final_contig_writer->MPIFileWrite();
    final_contig_writer->allreduce();
  } else {
    contig_writer->MPIFileWrite();
    contig_writer->allreduce();
  }
  graph.Delete_db();
}