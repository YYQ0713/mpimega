//
// Created by vout on 6/24/19.
//

#ifndef MEGAHIT_CONTIG_WRITER_H
#define MEGAHIT_CONTIG_WRITER_H

#include <definitions.h>
#include <atomic>
#include <string>
#include "utils/utils.h"
#include <vector>
#include <mutex>
#include "mpienv/mpienv.hpp"
#define FMT_HEADER_ONLY 
#include "fmt/core.h"

class ContigWriter {
 public:
  explicit ContigWriter(const std::string file_name) : file_name_(file_name) {
    file_ = xfopen(file_name.c_str(), "w");
  }

  ~ContigWriter() {
    std::FILE *info_file = xfopen((file_name_ + ".info").c_str(), "w");
    pfprintf(info_file, "{} {}\n", n_contigs_.load(), n_bases_.load());
    fclose(info_file);
    fclose(file_);
  }

  void WriteContig(const std::string &ascii_contig, unsigned k_size,
                   long long id, int flag, double multi) {
    pfprintf(file_, ">k{}_{} flag={} multi={.4} len={}\n{s}\n", k_size, id,
             flag, multi, ascii_contig.length(), ascii_contig.c_str());
    n_contigs_.fetch_add(1, std::memory_order_relaxed);
    n_bases_.fetch_add(
        ascii_contig.length() + (flag & contig_flag::kLoop) ? 28 : 0,
        std::memory_order_relaxed);
  }

  void WriteLocalContig(const std::string &ascii_contig,
                        int64_t origin_contig_id, int strand,
                        int64_t contig_id) {
    pfprintf(file_, ">lc_{}_strand_{}_id_{} flag=0 multi=1\n{s}\n",
             origin_contig_id, strand, contig_id, ascii_contig.c_str());
    n_contigs_.fetch_add(1, std::memory_order_relaxed);
    n_bases_.fetch_add(ascii_contig.length(), std::memory_order_relaxed);
  }

 private:
  std::string file_name_;
  std::FILE *file_;
  std::atomic<int64_t> n_contigs_{0};
  std::atomic<int64_t> n_bases_{0};
};

class MPIContigWriter {
 public:
  explicit MPIContigWriter(const std::string file_name, int rank) : file_name_(file_name), rank_(rank) {
    write_buffer.reserve(BUFSIZE); // reserve 256 MB
    int rc = MPI_File_open(
        MPI_COMM_WORLD,             // 通信域（单进程可用 MPI_COMM_SELF）
        file_name_.c_str(),                 // 文件名
        MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
        MPI_INFO_NULL,
        &mpi_file
    );
    if (rc != MPI_SUCCESS) {
        std::cerr << "Error: Failed to open MPI file!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
  }

  ~MPIContigWriter() {
    if (rank_ == 0) {
      std::FILE *info_file = xfopen((file_name_ + ".info").c_str(), "w");
      pfprintf(info_file, "{} {}\n", n_contigs_.load(), n_bases_.load());
      fclose(info_file);
    }
    MPI_File_close(&mpi_file);
  }

  void WriteContig(const std::string &ascii_contig, unsigned k_size,
                   long long id, int flag, double multi) {
    // 格式化字符串（原 pfprintf 的逻辑）
    //std::string formatted_str = ">k" + std::to_string(k_size) + "_" + std::to_string(id) +
    //                            " flag=" + std::to_string(flag) + " multi=" + std::to_string() +
    //                            " len=" + std::to_string(ascii_contig.length()) + "\n" + ascii_contig + "\n";
    //std::string formatted_str = fmt::format(
    //    ">k{}_{} flag={} multi={:.4f} len={}\n{}\n",
    //    k_size, id, flag, multi, ascii_contig.length(), ascii_contig);
    // 加锁，将数据存入缓冲区
    {
      std::lock_guard<std::mutex> lock(buffer_mutex);
      // 直接格式化并追加（避免临时 string）
      fmt::format_to(
          std::back_inserter(write_buffer),  // 直接追加到 combined_data
          ">k{}_{} flag={} multi={:.4f} len={}\n{}\n",
          k_size, id, flag, multi, ascii_contig.length(), ascii_contig
      );
    }
    
    n_contigs_.fetch_add(1, std::memory_order_relaxed);
    n_bases_.fetch_add(
        ascii_contig.length() + (flag & contig_flag::kLoop) ? 28 : 0,
        std::memory_order_relaxed);
  }

  void WriteLocalContig(const std::string &ascii_contig,
                        int64_t origin_contig_id, int strand,
                        int64_t contig_id) {
    pfprintf(file_, ">lc_{}_strand_{}_id_{} flag=0 multi=1\n{s}\n",
             origin_contig_id, strand, contig_id, ascii_contig.c_str());
    n_contigs_.fetch_add(1, std::memory_order_relaxed);
    n_bases_.fetch_add(ascii_contig.length(), std::memory_order_relaxed);
  }

  void MPIFileWrite() {
    {
      std::lock_guard<std::mutex> lock(buffer_mutex);
      if (write_buffer.empty()) {
        return;
      }
      // 原子写入（自动追加到文件末尾）
      MPI_File_write_shared(
        mpi_file,
        write_buffer.data(),
        write_buffer.size(),
        MPI_CHAR,
        MPI_STATUS_IGNORE
      );
      write_buffer.clear();  // 清空缓冲区
    }
  }

  bool check_buf() {
    return write_buffer.size() >= 268435456;
  }

  void allreduce() {
    MPI_Allreduce(MPI_IN_PLACE, &n_contigs_, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &n_bases_, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  }

 private:
  std::string file_name_;
  std::FILE *file_;
  int rank_;
  std::atomic<int64_t> n_contigs_{0};
  std::atomic<int64_t> n_bases_{0};
  MPI_File mpi_file;                       // MPI 文件句柄
  std::string write_buffer;                // 缓冲区
  std::mutex buffer_mutex;                 // 保护缓冲区的互斥锁
  static const unsigned BUFSIZE = 268435456; // 256 MB
};

#endif  // MEGAHIT_CONTIG_WRITER_H
