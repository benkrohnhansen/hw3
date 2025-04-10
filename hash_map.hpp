#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include <iostream>

struct HashMap {
    using dist_data upcxx::dist_object<std::vector<kmer_pair>>;    
    dist_data local_data;
    upcxx::dist_object<std::vector<int>> local_used;

    size_t my_size;

    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);

    // Dist functionality 
    int get_target_rank(const std::string &hash);
    upcxx::future<> dist_insert(const kmer_pair& kmer);
};

HashMap::HashMap(size_t size) {
    // self_ptr = upcxx::new_<HashMap>(*this);
    my_size = size;
    data.resize(size);
    used.resize(size, 0);
}

int HashMap::get_target_rank(const std::string &hash) {
    return hash % upcxx::rank_n();
}

upcxx::future<> HashMap::dist_insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    return upcxx::rpc(get_target_rank(hash),
        // lambda to insert the key-value pair
        [](dist_data &ldata, const kmer_pair &kmer) {
        // insert into the local map at the target
        ldata->insert({kmer});
        }, local_data, kmer);
}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        uint64_t slot = (hash + probe++) % size();
        success = request_slot(slot);
        if (success) {
            write_slot(slot, kmer);
        }
    } while (!success && probe < size());
    return success;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    // uint64_t hash = key_kmer.hash();
    // int owner_rank = hash % upcxx::rank_n();

    // if (owner_rank == upcxx::rank_me()) {
    //     // Local lookup — same logic as serial
    //     uint64_t probe = 0;
    //     bool success = false;
    //     do {
    //         uint64_t slot = (hash + probe++) % size();
    //         if (slot_used(slot)) {
    //             val_kmer = read_slot(slot);
    //             if (val_kmer.kmer == key_kmer) {
    //                 success = true;
    //             }
    //         }
    //     } while (!success && probe < size());
    //     return success;
    // } 
    return false;
}


bool HashMap::slot_used(uint64_t slot) { return used[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data[slot]; }

bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0) {
        return false;
    } else {
        used[slot] = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return my_size; }
