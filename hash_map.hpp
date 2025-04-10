#ifndef HASH_MAP_HPP
#define HASH_MAP_HPP

#include <upcxx/upcxx.hpp>
#include <cstdint>
#include <cstring>
#include <vector>
#include <functional>
#include <thread>    // For std::this_thread::yield

// Repository files that define the k‑mer types.
#include "pkmer_t.hpp"
#include "kmer_t.hpp"   // defines kmer_pair

// ---------------------------------------------------------------------
// Distributed bucket structure for the hash table.
// ---------------------------------------------------------------------
struct bucket_t {
  pkmer_t key;        // The packed k‑mer key.
  char forwardExt;    // Forward extension.
  char backwardExt;   // Backward extension.
  bool valid;         // Indicates if this bucket holds valid data.

  bucket_t() : forwardExt('F'), backwardExt('F'), valid(false) {
    std::memset(&key, 0, sizeof(pkmer_t));
  }
};

// ---------------------------------------------------------------------
// Distributed hash table (hash_map) using open addressing with
// linear probing and a coarse‑grained global lock (via RPC on rank 0).
// ---------------------------------------------------------------------
class hash_map {
public:
  size_t table_size;
  // Distributed array (global pointer) to an array of buckets.
  upcxx::global_ptr<bucket_t> table;
  // Global lock pointer (one int) shared by all ranks.
  upcxx::global_ptr<int> global_lock;

  // Constructor: allocate and initialize the table and global lock.
  hash_map(size_t size) : table_size(size) {
    // Allocate the distributed table.
    table = upcxx::new_array<bucket_t>(table_size);
    // Rank 0 creates the global lock.
    if (upcxx::rank_me() == 0) {
      global_lock = upcxx::new_<int>(0);
    }
    // Broadcast the global_lock pointer so that all ranks receive it.
    global_lock = upcxx::broadcast(global_lock, 0).wait();

    // Initialize each bucket.
    for (size_t i = 0; i < table_size; i++) {
      bucket_t* bp = (table + i).local();
      bp->valid = false;
      bp->forwardExt = 'F';
      bp->backwardExt = 'F';
      std::memset(&(bp->key), 0, sizeof(pkmer_t));
    }
    upcxx::barrier();  // Ensure every rank sees the initialized table.
  }

  ~hash_map() {
    upcxx::delete_array(table);
    if (upcxx::rank_me() == 0) {
      upcxx::delete_(global_lock);
    }
  }

  // Simple hash function for a pkmer_t key.
  size_t hash(const pkmer_t &key) const {
    const uint64_t* p = reinterpret_cast<const uint64_t*>(&key);
    return std::hash<uint64_t>()(*p) % table_size;
  }

private:
  // ---------------------------------------------------------------------
  // Global lock acquisition implemented via an RPC on rank 0.
  // ---------------------------------------------------------------------
  void lock_global() {
    while (true) {
      bool acquired = upcxx::rpc(
        0,
        [] (upcxx::global_ptr<int> lock_ptr) -> bool {
          // Access the lock’s underlying storage on rank 0.
          int old = *(lock_ptr.local());
          if(old == 0) {
            *(lock_ptr.local()) = 1;
            return true;
          }
          return false;
        },
        global_lock
      ).wait();
      if(acquired)
        break;
      std::this_thread::yield();
    }
  }

  void unlock_global() {
    upcxx::rpc(
      0,
      [] (upcxx::global_ptr<int> lock_ptr) {
         *(lock_ptr.local()) = 0;
      },
      global_lock
    ).wait();
  }

public:
  // ---------------------------------------------------------------------
  // Insert the given bucket (bk) into the distributed hash table.
  // Returns true if insertion is successful; false if the key already
  // exists or the table is full.
  // ---------------------------------------------------------------------
  bool insert(const bucket_t &bk) {
    lock_global();
    size_t idx = hash(bk.key);
    size_t start_idx = idx;
    bucket_t curr;
    bool result = false;
    while (true) {
      curr = upcxx::rget(table + idx).wait();
      if (!curr.valid) {
        // Bucket is empty; write our bucket.
        upcxx::rput(bk, table + idx).wait();
        result = true;
        break;
      }
      // Compare keys using the member operator== on pkmer_t.
      if (curr.valid && (curr.key.operator==(bk.key))) {
        result = false;
        break;
      }
      idx = (idx + 1) % table_size;
      if (idx == start_idx) {
        result = false;
        break;
      }
    }
    unlock_global();
    return result;
  }

  // ---------------------------------------------------------------------
  // Find the bucket with the given key.
  // If found, copies its data into result and returns true;
  // otherwise returns false.
  // ---------------------------------------------------------------------
  bool find(const pkmer_t &key, bucket_t &result) {
    lock_global();
    size_t idx = hash(key);
    size_t start_idx = idx;
    bucket_t curr;
    bool found = false;
    while (true) {
      curr = upcxx::rget(table + idx).wait();
      if (!curr.valid) {
        found = false;
        break;
      }
      if (curr.valid && (curr.key.operator==(key))) {
        result = curr;
        found = true;
        break;
      }
      idx = (idx + 1) % table_size;
      if (idx == start_idx) {
        found = false;
        break;
      }
    }
    unlock_global();
    return found;
  }
};

#endif // HASH_MAP_HPP
