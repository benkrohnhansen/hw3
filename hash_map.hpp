#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    std::vector<kmer_pair> data;
    std::vector<int> used;

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
};

HashMap::HashMap(size_t size) {
    my_size = size;
    data.resize(size);
    used.resize(size, 0);
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
    uint64_t hash = key_kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        uint64_t slot = (hash + probe++) % size();
        if (slot_used(slot)) {
            val_kmer = read_slot(slot);
            if (val_kmer.kmer == key_kmer) {
                success = true;
            }
        }
    } while (!success && probe < size());
    return success;
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


// class DistributedHashMap {
//     private:
//         upcxx::dist_object<upcxx::global_ptr<HashMap>> local_map_g;
//         HashMap *local_map;
        
//         int get_target_rank(const std::string &key) {
//             return std::hash<std::string>{}(key) % upcxx::rank_n();
//         }
//     public:
//         DistributedHashMap(size_t local_size)
//             : local_map_g(HashMap(local_size)) {
//                 local_map = local_map_g->local();
//             }

//         bool insert(const kmer_pair& kmer) {
//             return true;
//         }

//         size_t size() {
//             return local_map->size();
//         }
// };

class DistributedHashMap {
    private:
        using dist_hash_map = upcxx::dist_object<upcxx::global_ptr<HashMap>>;
        dist_hash_map local_map_g;
        HashMap* local_map;
    
        int get_target_rank(const uint64_t &hash) {
            return hash % upcxx::rank_n();
        }
    
    public:
        DistributedHashMap(size_t local_size)
            : local_map_g(upcxx::new_<HashMap>(local_size)) {
            // Only valid on the local rank!
            local_map = local_map_g->local();
        }

        upcxx::future<> insert(const kmer_pair& kmer) {
        uint64_t hash = kmer.hash();
        return upcxx::rpc(get_target_rank(hash),
            // lambda to insert the key-value pair
            [](dist_hash_map &lmap, const kmer_pair &kmer) {
            // insert into the local map at the target
            HashMap* local = lmap->local();
            local->insert(kmer);
            }, local_map_g, kmer);
        }   
    
        size_t size() const {
            return local_map->size();
        }
};


// #pragma once

// #include "kmer_t.hpp"
// #include <upcxx/upcxx.hpp>
// #include <iostream>

// struct HashMap {
//     using dist_data = upcxx::dist_object<std::vector<kmer_pair>>;    
//     dist_data local_data;

//     using dist_used = upcxx::dist_object<std::vector<int>>;
//     dist_used local_used;

//     size_t my_size;

//     size_t size() const noexcept;

//     HashMap(size_t size);

//     // Most important functions: insert and retrieve
//     // k-mers from the hash table.
//     bool insert(const kmer_pair& kmer);
//     bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

//     // Helper functions

//     // Write and read to a logical data slot in the table.
//     void write_slot(uint64_t slot, const kmer_pair& kmer);
//     kmer_pair read_slot(uint64_t slot);

//     // Request a slot or check if it's already used.
//     bool request_slot(uint64_t slot);
//     bool slot_used(uint64_t slot);

//     // Dist functionality 
//     int get_target_rank(const std::string &hash);
//     upcxx::future<> dist_insert(const kmer_pair& kmer);
// };

// HashMap::HashMap(size_t size) {
//     // self_ptr = upcxx::new_<HashMap>(*this);
//     my_size = size;
//     local_data.resize(size);
//     local_used.resize(size, 0);
// }

// int HashMap::get_target_rank(const std::string &hash) {
//     return hash % upcxx::rank_n();
// }

// upcxx::future<> HashMap::dist_insert(const kmer_pair& kmer) {
//     uint64_t hash = kmer.hash();
//     return upcxx::rpc(get_target_rank(hash),
//         // lambda to insert the key-value pair
//         [](dist_data &ldata, const kmer_pair &kmer) {
//         // insert into the local map at the target
//         ldata->insert({kmer});
//         }, local_data, kmer);
// }

// bool HashMap::insert(const kmer_pair& kmer) {
//     uint64_t hash = kmer.hash();
//     uint64_t probe = 0;
//     bool success = false;
//     do {
//         uint64_t slot = (hash + probe++) % size();
//         success = request_slot(slot);
//         if (success) {
//             write_slot(slot, kmer);
//         }
//     } while (!success && probe < size());
//     return success;
// }

// bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
//     // uint64_t hash = key_kmer.hash();
//     // int owner_rank = hash % upcxx::rank_n();

//     // if (owner_rank == upcxx::rank_me()) {
//     //     // Local lookup — same logic as serial
//     //     uint64_t probe = 0;
//     //     bool success = false;
//     //     do {
//     //         uint64_t slot = (hash + probe++) % size();
//     //         if (slot_used(slot)) {
//     //             val_kmer = read_slot(slot);
//     //             if (val_kmer.kmer == key_kmer) {
//     //                 success = true;
//     //             }
//     //         }
//     //     } while (!success && probe < size());
//     //     return success;
//     // } 
//     return false;
// }


// bool HashMap::slot_used(uint64_t slot) { return used[slot] != 0; }

// void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data[slot] = kmer; }

// kmer_pair HashMap::read_slot(uint64_t slot) { return data[slot]; }

// bool HashMap::request_slot(uint64_t slot) {
//     if (used[slot] != 0) {
//         return false;
//     } else {
//         used[slot] = 1;
//         return true;
//     }
// }

// size_t HashMap::size() const noexcept { return my_size; }
