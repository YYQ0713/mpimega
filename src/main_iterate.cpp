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

#include <assert.h>
#include <omp.h>
#include <stdio.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "definitions.h"
#include "iterate/contig_flank_index.h"
#include "iterate/kmer_collector.h"
#include "sequence/io/async_sequence_reader.h"
#include "utils/options_description.h"
#include "libbloom/bloom.h"
#include <rocksdb/db.h>
#include <rocksdb/iterator.h>
#include <rocksdb/options.h>
#include <rocksdb/table.h>
#include <rocksdb/filter_policy.h>

using std::string;
using std::vector;

namespace {

struct IterOption {
  string contig_file;
  string bubble_file;
  string read_file;
  int num_cpu_threads{0};
  int kmer_k{0};
  int step{0};
  string output_prefix;
} opt;

// static void ParseIterOptions(int argc, char *argv[]) {
//   OptionsDescription desc;

//   desc.AddOption("contig_file", "c", opt.contig_file,
//                  "(*) contigs file, fasta/fastq format, output by assembler");
//   desc.AddOption("bubble_file", "b", opt.bubble_file,
//                  "(*) bubble file, fasta/fastq format, output by assembler");
//   desc.AddOption("read_file", "r", opt.read_file,
//                  "(*) reads to be aligned. \"-\" for stdin. Can be gzip'ed.");
//   desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads,
//                  "number of cpu threads, at least 2. 0 for auto detect.");
//   desc.AddOption("kmer_k", "k", opt.kmer_k, "(*) current kmer size.");
//   desc.AddOption("step", "s", opt.step,
//                  "(*) step for iteration (<= 28). i.e. this iteration is from "
//                  "kmer_k to (kmer_k + step)");
//   desc.AddOption("output_prefix", "o", opt.output_prefix,
//                  "(*) output_prefix.edges.0 will be created.");

//   try {
//     desc.Parse(argc, argv);

//     if (opt.step + opt.kmer_k >=
//         static_cast<int>(
//             std::max(Kmer<4>::max_size(), GenericKmer::max_size()))) {
//       std::ostringstream os;
//       os << "kmer_k + step must less than "
//          << std::max(Kmer<4>::max_size(), GenericKmer::max_size());
//       throw std::logic_error(os.str());
//     } else if (opt.contig_file == "") {
//       throw std::logic_error("No contig file!");
//     } else if (opt.bubble_file == "") {
//       throw std::logic_error("No bubble file!");
//     } else if (opt.read_file == "") {
//       throw std::logic_error("No reads file!");
//     } else if (opt.kmer_k <= 0) {
//       throw std::logic_error("Invalid kmer size!");
//     } else if (opt.step <= 0 || opt.step > 28 || opt.step % 2 == 1) {
//       throw std::logic_error("Invalid step size!");
//     } else if (opt.output_prefix == "") {
//       throw std::logic_error("No output prefix!");
//     }
//     if (opt.num_cpu_threads == 0) {
//       opt.num_cpu_threads = omp_get_max_threads();
//     }
//     if (opt.num_cpu_threads > 1) {
//       omp_set_num_threads(opt.num_cpu_threads - 1);
//     } else {
//       omp_set_num_threads(1);
//     }
//   } catch (std::exception &e) {
//     std::cerr << e.what() << std::endl;
//     std::cerr << "Usage: " << argv[0] << " [opt]" << std::endl;
//     std::cerr << "opt with (*) are must" << std::endl;
//     std::cerr << "opt:" << std::endl;
//     std::cerr << desc << std::endl;
//     exit(1);
//   }
// }

}  // namespace

void ParseItOptions(int argc, char **argv, IterOption &opt) {
  OptionsDescription desc;

  desc.AddOption("contig_file", "c", opt.contig_file,
                 "(*) contigs file, fasta/fastq format, output by assembler");
  desc.AddOption("bubble_file", "b", opt.bubble_file,
                 "(*) bubble file, fasta/fastq format, output by assembler");
  desc.AddOption("read_file", "r", opt.read_file,
                 "(*) reads to be aligned. \"-\" for stdin. Can be gzip'ed.");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads,
                 "number of cpu threads, at least 2. 0 for auto detect.");
  desc.AddOption("kmer_k", "k", opt.kmer_k, "(*) current kmer size.");
  desc.AddOption("step", "s", opt.step,
                 "(*) step for iteration (<= 28). i.e. this iteration is from "
                 "kmer_k to (kmer_k + step)");
  desc.AddOption("output_prefix", "o", opt.output_prefix,
                 "(*) output_prefix.edges.0 will be created.");

  for (int i = 1; i < argc; ++i) {
    std::string option = argv[i];
    if (option == "-c" || option == "--contig_file") {
      if (i + 1 <= argc) {
        opt.contig_file = argv[++i];
      }
    }
    else if (option == "-b" || option == "--bubble_file") {
      if (i + 1 <= argc) {
        opt.bubble_file = argv[++i];
      }
    }
    else if (option == "-r" || option == "--read_file") {
      if (i + 1 <= argc) {
        opt.read_file = argv[++i];
      }
    }
    else if (option == "-k" || option == "--kmer_k") {
      if (i + 1 <= argc) {
        opt.kmer_k = std::stoi(argv[++i]);
      }
    }
    else if (option == "-s" || option == "--step") {
      if (i + 1 <= argc) {
        opt.step = std::stoi(argv[++i]);
      }
    }
    else if (option == "-t" || option == "--num_cpu_threads") {
      if (i + 1 <= argc) {
        opt.num_cpu_threads = std::stoi(argv[++i]);
      }
    }
    else if (option == "-o" || option == "--output_prefix") {
      if (i + 1 <= argc) {
        opt.output_prefix = argv[++i];
      }
    }
  }

  try {
    if (opt.step + opt.kmer_k >=
        static_cast<int>(
            std::max(Kmer<4>::max_size(), GenericKmer::max_size()))) {
      std::ostringstream os;
      os << "kmer_k + step must less than "
         << std::max(Kmer<4>::max_size(), GenericKmer::max_size());
      throw std::logic_error(os.str());
    } else if (opt.contig_file == "") {
      throw std::logic_error("No contig file!");
    } else if (opt.bubble_file == "") {
      throw std::logic_error("No bubble file!");
    } else if (opt.read_file == "") {
      throw std::logic_error("No reads file!");
    } else if (opt.kmer_k <= 0) {
      throw std::logic_error("Invalid kmer size!");
    } else if (opt.step <= 0 || opt.step > 28 || opt.step % 2 == 1) {
      throw std::logic_error("Invalid step size!");
    } else if (opt.output_prefix == "") {
      throw std::logic_error("No output prefix!");
    }
    if (opt.num_cpu_threads == 0) {
      opt.num_cpu_threads = omp_get_max_threads();
    }
    if (opt.num_cpu_threads > 1) {
      omp_set_num_threads(opt.num_cpu_threads - 1);
    } else {
      omp_set_num_threads(1);
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " [opt]" << std::endl;
    std::cerr << "opt with (*) are must" << std::endl;
    std::cerr << "opt:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
};

template <class KmerType, class IndexType>
static bool ReadReadsAndProcessKernel(const IterOption &opt,
                                      const IndexType &index,
                                      MPIEnviroment &mpienv) {
  if (KmerType::max_size() < static_cast<unsigned>(opt.kmer_k + opt.step + 1)) {
    return false;
  }
  xinfo("Selected kmer type size for next k: {}\n", sizeof(KmerType));
  BinaryReader binary_reader(opt.read_file);
  AsyncSequenceReader reader(&binary_reader);
  KmerCollector<KmerType> collector(opt.kmer_k + opt.step + 1,
                                    opt.output_prefix,
                                    mpienv);
  int64_t num_aligned_reads = 0;
  int64_t num_total_reads = 0;
  int64_t num_iterative_edges = 0;
  MPIEdgeWriter<KmerType> mpi_edgewiriter(opt.kmer_k, opt.output_prefix, mpienv);

  rocksdb::DB* db;
  rocksdb::Options options;
  options.create_if_missing = true;

  options.IncreaseParallelism(opt.num_cpu_threads / 2);
  options.max_background_jobs = opt.num_cpu_threads / 2; 

  // Write Buffer & Memtable 设置
  options.write_buffer_size = 1ULL * 1024ULL * 1024ULL * 1024ULL; // 每个 memtable 1GB
  options.max_write_buffer_number = 16;          // 最多 16 个 memtable
  options.min_write_buffer_number_to_merge = 2; // 减少 flush 频率

  // Level-Compaction 优化
  options.OptimizeUniversalStyleCompaction();   // Universal compaction 更适合写多
  options.compaction_style = rocksdb::kCompactionStyleUniversal;
  options.compaction_options_universal.size_ratio = 20;
  options.compaction_options_universal.min_merge_width = 2;
  options.compaction_options_universal.max_size_amplification_percent = 300;

  // 文件大小 & compaction 触发阈值
  options.target_file_size_base = 1024ULL * 1024ULL * 1024ULL; // 每个 SST 文件更大，减少文件数量
  options.level0_file_num_compaction_trigger = 8;    // L0 文件超过 8 就触发 compaction
  options.level0_slowdown_writes_trigger = 20;       // 超过 20 慢写
  options.level0_stop_writes_trigger = 40;           // 超过 40 停止写入，保护系统

  // 压缩算法选择
  options.compression = rocksdb::kNoCompression;     // 写入密集，直接不压缩
  options.bottommost_compression = rocksdb::kLZ4Compression; // 只有最后一层压缩，节省空间

  // I/O 优化
  options.use_direct_io_for_flush_and_compaction = true;
  options.use_direct_reads = true;
  options.bytes_per_sync = 16 * 1024 * 1024;  // 每写 2MB 刷一次，平滑写入
  options.wal_bytes_per_sync = 1 * 1024 * 1024;

  // Table Format & Cache
  rocksdb::BlockBasedTableOptions table_options;
  table_options.block_cache = rocksdb::NewLRUCache(1024ULL * 1024ULL * 1024ULL); // 缓存 1GB block
  table_options.cache_index_and_filter_blocks = true;
  table_options.pin_l0_filter_and_index_blocks_in_cache = true;
  table_options.filter_policy.reset(rocksdb::NewBloomFilterPolicy(10, false)); // 可选
  options.table_factory.reset(NewBlockBasedTableFactory(table_options));

  // WriteOptions (批量写入场景可以关 WAL)
  rocksdb::WriteOptions write_options;
  write_options.disableWAL = true; 
  write_options.sync = false;      // 不强制同步 WAL
  // options.OptimizeLevelStyleCompaction();

  std::string db_path = mpi_edgewiriter.RetFilePrefix() + ".rank." + std::to_string(mpienv.rank);
  rocksdb::Status status = rocksdb::DB::Open(options, db_path, &db);
  if (!status.ok()) {
      std::cerr << "Unable to open database: " << status.ToString() << std::endl;
      return 1;
  }

  //Bloom bloom(1000000000, 0.000001);

  while (true) {
    const auto &read_pkg = reader.Next();
    if (read_pkg.seq_count() == 0) {
      break;
    }
    num_aligned_reads += index.FindNextKmersFromReads(read_pkg, &collector, mpienv.rank, mpienv.nprocs, &mpi_edgewiriter, &num_iterative_edges, db);
    num_total_reads += read_pkg.seq_count();
    // xinfo("Processed: {}, aligned: {}. Iterative edges: {}\n", num_total_reads,
    //        num_aligned_reads, num_iterative_edges);
    xinfo("Processed: {}, aligned: {}.\n", num_total_reads, num_aligned_reads);
  }
  
  
  xinfo("finish db construct and start writing file\n");
  rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    KmerType tmp_kmer(it->key().data());
    mul_t mul = *reinterpret_cast<const mul_t*>(it->value().data());
    mpi_edgewiriter.WriteToBuf(tmp_kmer, mul);
    num_iterative_edges++;
    if (mpi_edgewiriter.check_buf()) {
      mpi_edgewiriter.MPIFileWrite();
    }
  }
  xinfo("Total: {}, aligned: {}. Iterative edges: {}\n", num_total_reads, num_aligned_reads, num_iterative_edges);
  //collector.FlushToFile();
  mpi_edgewiriter.MPIFileWrite();
  mpi_edgewiriter.allreduce();
  mpi_edgewiriter.Finalize(mpienv);

  delete db;

  rocksdb::Options cleanup_options;
  status = rocksdb::DestroyDB(db_path, cleanup_options);
  if (!status.ok()) {
      std::cerr << "Failed to clean up database: " << status.ToString() << std::endl;
  }
  
  return true;
}

template <class IndexType>
static void ReadReadsAndProcess(const IterOption &opt,
                                const IndexType &index,
                                MPIEnviroment &mpienv) {
  if (ReadReadsAndProcessKernel<Kmer<1, uint64_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<3, uint32_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<2, uint64_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<5, uint32_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<3, uint64_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<7, uint32_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<4, uint64_t>>(opt, index, mpienv)) return;
  if (ReadReadsAndProcessKernel<Kmer<kUint32PerKmerMaxK, uint32_t>>(opt, index, mpienv))
    return;
  xfatal("k is too large!\n");
}

template <class IndexType>
static void ReadContigsAndBuildIndex(const IterOption &opt,
                                     const std::string &file_name,
                                     IndexType *index,
                                     MPIEnviroment &mpienv) {
  AsyncContigReader reader(file_name);
  while (true) {
    auto &pkg = reader.Next();
    auto &contig_pkg = pkg.first;
    auto &mul = pkg.second;
    if (contig_pkg.seq_count() == 0) {
      break;
    }
    xinfo("Read {} contigs\n", contig_pkg.seq_count());
    index->FeedBatchContigs(contig_pkg, mul);
    xinfo("Number of flank kmers: {}\n", index->size());
  }
}

struct BaseRunner {
  virtual ~BaseRunner() = default;
  virtual void Run(const IterOption &opt, MPIEnviroment &mpienv) = 0;
  virtual uint32_t max_k() const = 0;
};

template <class KmerType>
struct Runner : public BaseRunner {
  ~Runner() override = default;
  void Run(const IterOption &opt, MPIEnviroment &mpienv) override {
    xinfo("Selected kmer type size for k: {}\n", sizeof(KmerType));
    ContigFlankIndex<KmerType> index(opt.kmer_k, opt.step);
    ReadContigsAndBuildIndex(opt, opt.contig_file, &index, mpienv);
    ReadContigsAndBuildIndex(opt, opt.bubble_file, &index, mpienv);
    ReadReadsAndProcess(opt, index, mpienv);
    size_t vmrss_kb = getCurrentRSS_kb();
    xinfo("End of iter currentRSS: {} KB\n", vmrss_kb);
  }
  uint32_t max_k() const override { return KmerType::max_size(); }
};

int main_iterate(int argc, char **argv, MPIEnviroment &mpienv) {
  size_t vmrss_kb = getCurrentRSS_kb();
  xinfo("start of iter currentRSS: {} KB\n", vmrss_kb);
  AutoMaxRssRecorder recorder;
  IterOption opt;
  ParseItOptions(argc, argv, opt);

  std::list<std::unique_ptr<BaseRunner>> runners;
  runners.emplace_back(new Runner<Kmer<1, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<3, uint32_t>>());
  runners.emplace_back(new Runner<Kmer<2, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<5, uint32_t>>());
  runners.emplace_back(new Runner<Kmer<3, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<7, uint32_t>>());
  runners.emplace_back(new Runner<Kmer<4, uint64_t>>());
  runners.emplace_back(new Runner<Kmer<kUint32PerKmerMaxK, uint32_t>>());

  for (auto &runner : runners) {
    if (runner->max_k() >= static_cast<uint32_t>(opt.kmer_k + 1)) {
      runner->Run(opt, mpienv);
      return 0;
    }
  }

  xfatal("k is too large!\n");
  return 1;
}