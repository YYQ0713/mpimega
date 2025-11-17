//
// Created by vout on 11/5/18.
//

#include "sdbg_writer.h"

#include <utils/utils.h>
#include <algorithm>
#include <cassert>
#include "sdbg_item.h"

void SdbgWriter::InitFiles(MPIEnviroment &mpienv) {
  cur_thread_offset_.resize(num_threads_, 0);
  bucket_rec_.resize(num_buckets_);

  for (size_t i = 0; i < num_threads_; ++i) {
    files_.emplace_back(
        new std::ofstream((file_prefix_ + ".rank." + std::to_string(mpienv.rank) + ".sdbg." + std::to_string(i)).c_str(),
                          std::ofstream::binary | std::ofstream::out));
    assert(files_[i]->is_open());
  }
  is_opened_ = true;
}

void SdbgWriter::InitFiles() {
 cur_thread_offset_.resize(num_threads_, 0);
 bucket_rec_.resize(num_buckets_);

 for (size_t i = 0; i < num_threads_; ++i) {
   files_.emplace_back(
       new std::ofstream((file_prefix_ + ".sdbg." + std::to_string(i)).c_str(),
                         std::ofstream::binary | std::ofstream::out));
   assert(files_[i]->is_open());
 }
 is_opened_ = true;
}

/*
void SdbgWriter::Write(unsigned tid, uint32_t bucket_id, uint8_t w,
                       uint8_t last, uint8_t tip, mul_t multiplicity,
                       label_word_t *packed_tip_label, Snapshot *snapshot) {
  if (bucket_id != snapshot->bucket_record.bucket_id) {
    assert(snapshot->bucket_record.bucket_id == SdbgBucketRecord::kNullID);
    assert(snapshot->bucket_record.file_id == SdbgBucketRecord::kNullID);
    snapshot->bucket_record.file_id = tid;
    snapshot->bucket_record.bucket_id = bucket_id;
    snapshot->bucket_record.starting_offset = cur_thread_offset_[tid];
    snapshot->cur_thread_offset = cur_thread_offset_[tid];
  }
  assert(tid < num_threads_ && tid == snapshot->bucket_record.file_id);

  SdbgItem item(w, last, tip, std::min(multiplicity, mul_t{kSmallMulSentinel}));
  files_[tid]->write(reinterpret_cast<const char *>(&item), sizeof(item));
  ++snapshot->bucket_record.num_items;
  ++snapshot->bucket_record.num_w[w];
  snapshot->bucket_record.ones_in_last += last;
  snapshot->cur_thread_offset += sizeof(uint16_t);

  if (multiplicity > kMaxSmallMul) {
    files_[tid]->write(reinterpret_cast<const char *>(&multiplicity),
                       sizeof(multiplicity));
    ++snapshot->bucket_record.num_large_mul;
    snapshot->cur_thread_offset += sizeof(mul_t);
  }

  if (tip) {
    files_[tid]->write(reinterpret_cast<const char *>(packed_tip_label),
                       sizeof(label_word_t) * words_per_tip_label_);
    ++snapshot->bucket_record.num_tips;
    snapshot->cur_thread_offset += sizeof(label_word_t) * words_per_tip_label_;
  }
}
*/

void SdbgWriter::Write(unsigned rank, unsigned nthread, unsigned tid, uint32_t bucket_id, uint8_t w,
                       uint8_t last, uint8_t tip, mul_t multiplicity,
                       label_word_t *packed_tip_label, Snapshot *snapshot) {
  if (bucket_id != snapshot->bucket_record.bucket_id) {
    assert(snapshot->bucket_record.bucket_id == SdbgBucketRecord::kNullID);
    assert(snapshot->bucket_record.file_id == SdbgBucketRecord::kNullID);
    snapshot->bucket_record.file_id = rank * nthread + tid;
    snapshot->bucket_record.bucket_id = bucket_id;
    snapshot->bucket_record.starting_offset = cur_thread_offset_[tid];
    snapshot->cur_thread_offset = cur_thread_offset_[tid];
  }
  assert(tid < num_threads_ && ((rank * nthread + tid) == snapshot->bucket_record.file_id));

  SdbgItem item(w, last, tip, std::min(multiplicity, mul_t{kSmallMulSentinel}));
  files_[tid]->write(reinterpret_cast<const char *>(&item), sizeof(item));
  //MPI_Request req1;
  //MPI_File_write(files_[tid], reinterpret_cast<const char *>(&item), 1, MPI_UINT16_T, MPI_STATUS_IGNORE);
  //requests.push_back(req1);
  ++snapshot->bucket_record.num_items;
  ++snapshot->bucket_record.num_w[w];
  snapshot->bucket_record.ones_in_last += last;
  snapshot->cur_thread_offset += sizeof(uint16_t);

  if (multiplicity > kMaxSmallMul) {
    //MPI_Request req2;
    //MPI_File_write(files_[tid], reinterpret_cast<const char *>(&multiplicity), 1, MPI_UINT16_T, MPI_STATUS_IGNORE);
    //requests.push_back(req2);
    files_[tid]->write(reinterpret_cast<const char *>(&multiplicity),
                       sizeof(multiplicity));
    ++snapshot->bucket_record.num_large_mul;
    snapshot->cur_thread_offset += sizeof(mul_t);
  }

  if (tip) {
    //MPI_Request req3;
    //MPI_File_write(files_[tid], reinterpret_cast<const char *>(&packed_tip_label), words_per_tip_label_, MPI_UINT32_T, MPI_STATUS_IGNORE);
    //requests.push_back(req3);

    files_[tid]->write(reinterpret_cast<const char *>(packed_tip_label),
                       sizeof(label_word_t) * words_per_tip_label_);
    ++snapshot->bucket_record.num_tips;
    snapshot->cur_thread_offset += sizeof(label_word_t) * words_per_tip_label_;
  }
}

void SdbgWriter::SaveSnapshot(const SdbgWriter::Snapshot &snapshot, int nthread) {
  if (snapshot.bucket_record.file_id != SdbgBucketRecord::kNullID) {
    bucket_rec_[snapshot.bucket_record.bucket_id] = snapshot.bucket_record;
    cur_thread_offset_[snapshot.bucket_record.file_id % nthread] =
        snapshot.cur_thread_offset;
  }
}

void SdbgWriter::Finalize() {
  if (is_opened_) {
    for (size_t i = 0; i < num_threads_; ++i) {
      files_[i]->close();
    }

    //for (auto& fh : files_) {
    //  MPI_File_close(&fh);
    //}

    is_opened_ = false;
  }
}

void SdbgWriter::Finalize(MPIEnviroment &mpienv) {
  if (is_opened_) {
    for (size_t i = 0; i < num_threads_; ++i) {
      files_[i]->close();
    }

    //for (auto& fh : files_) {
    //  MPI_File_close(&fh);
    //}

    if (mpienv.rank == 0)
    {
      std::ofstream os((file_prefix_ + ".sdbg_info").c_str());
      final_meta_.FromBucketRecord(bucket_rec_, k_, words_per_tip_label_)
          .Serialize(os);
      os.close();
    }
    
    is_opened_ = false;
  }
}