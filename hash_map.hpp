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
    // bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    kmer_pair find(const pkmer_t& key_kmer);
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
        // std::cout << "Rank " << upcxx::rank_me() << " insert success: " << success << std::endl;
        if (success) {
            write_slot(slot, kmer);
        }
    } while (!success && probe < size());
    return success;
}

// bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
kmer_pair HashMap::find(const pkmer_t& key_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t probe = 0;
    kmer_pair val_kmer;
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
    // return success;
    return val_kmer;
}

bool HashMap::slot_used(uint64_t slot) { return used[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data[slot]; }

bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0) {
        return false;
    } else {
        used[slot] = 1;
        // std::cout << "Rank " << upcxx::rank_me() << " used " << slot << " to 1\n"; 
        return true;
    }
}

size_t HashMap::size() const noexcept { return my_size; }


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
            // std::cout << "from " << upcxx::rank_me() << " to " << get_target_rank(hash) << std::endl;
            return upcxx::rpc(get_target_rank(hash),
                // lambda to insert the key-value pair
                [](dist_hash_map &lmap, const kmer_pair &kmer) {
                // insert into the local map at the target
                HashMap* local = lmap->local();
                local->insert(kmer);
                }, local_map_g, kmer);
        }   

        // upcxx::future<kmer_pair> find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
        upcxx::future<kmer_pair> find(const pkmer_t& key_kmer) {
            uint64_t hash = key_kmer.hash();
            return upcxx::rpc(get_target_rank(hash),
                // lambda to find the hash in the local map
                [](dist_hash_map &lmap, const pkmer_t& key_kmer) -> kmer_pair {
                HashMap* local = lmap->local();
                return local->find(key_kmer);
                // if (val_kmer == lmap->end()) return kmer_pair; // not found
                }, local_map_g, key_kmer);
        }
    
        size_t size() const {
            return local_map->size();
        }
};
