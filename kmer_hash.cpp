#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"

#include "butil.hpp"

#include <iostream>

#include <map>

using namespace std;

class DistrMap {
    private:
      // store the local unordered map in a distributed object to access from RPCs
      using dobj_map_t = upcxx::dist_object<std::unordered_map<std::string, std::string> >;
      dobj_map_t local_map;
      // map the key to a target process
      int get_target_rank(const std::string &key) {
        return std::hash<std::string>{}(key) % upcxx::rank_n();
      }
    public:
      // initialize the local map
      DistrMap() : local_map({}) {}
      // insert a key-value pair into the hash table
      upcxx::future<> insert(const std::string &key, const std::string &val) {
        // the RPC returns an empty upcxx::future by default
        return upcxx::rpc(get_target_rank(key),
                          // lambda to insert the key-value pair
                          [](dobj_map_t &lmap, const std::string &key, const std::string &val) {
                            // insert into the local map at the target
                            lmap->insert({key, val});
                          }, local_map, key, val);
      }
      // find a key and return associated value in a future
      upcxx::future<std::string> find(const std::string &key) {
        return upcxx::rpc(get_target_rank(key),
                          // lambda to find the key in the local map
                          [](dobj_map_t &lmap, const std::string &key) -> std::string {
                            auto elem = lmap->find(key);
                            if (elem == lmap->end()) return std::string(); // not found
                            else return elem->second; // key found: return value
                          }, local_map, key);
      }
    };




int main(int argc, char** argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = "";

    if (argc >= 3) {
        run_type = std::string(argv[2]);
    }

    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = std::string(argv[3]);
    }

    int ks = kmer_size(kmer_fname);

    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) +
                                 "-mers.  Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);

    // Load factor of 0.5
    size_t hash_table_size = n_kmers * (1.0 / 0.5);
    size_t local_size = (hash_table_size / upcxx::rank_n());
    if (hash_table_size % upcxx::rank_n() > upcxx::rank_me()) {
        local_size++;
    }
    // DistributedHashMap hashmap(local_size);
    DistributedHashMap hashmap(local_size);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
                     n_kmers);
    }

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }

    upcxx::barrier();

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;

    if (upcxx::rank_me() == 0) {
        std::cout << "Size " << hashmap.size() << std::endl;
    }

    size_t rank_start = upcxx::rank_me() * (hash_table_size / upcxx::rank_n()) 
                     + hash_table_size % upcxx::rank_n();
    size_t rank_end = rank_start + local_size;

    std::cout << "rank " << upcxx::rank_me() << "Rank start " << rank_start << " rank end " << rank_end << "\n";
    
    upcxx::barrier();
    std::vector<upcxx::future<>> futures;

    for (size_t i = rank_start; i < rank_end; i++) {
        auto kmer = kmers[i];
        futures.push_back(hashmap.insert(kmers[i]));

        // std::cout << "kmer " << kmer.kmer_str() << std::endl;

        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
            std::cout << "rank " << upcxx::rank_me() << " kmer " << kmer.kmer_str() << " " << kmer.backwardExt << std::endl;
        }
        // if (i > 20) break;
    }

    upcxx::when_all(futures.begin(), futures.end()).wait();
    // if (kmer.backwardExt() == 'F') {
    //     start_nodes.push_back(kmer);
    // }

    upcxx::barrier();
    
    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    for (auto& start_node : start_nodes) {
        std::cout << "rank " << upcxx::rank_me() << " start node " << start_node.kmer_str() << std::endl;
        std::cout << start_node.backwardExt() << std::endl; 
    }
    
    
    // if (upcxx::rank_me() == 1) {
    //         hashmap.print_local();
    //     }
    
    // upcxx::barrier();
        
    // if (upcxx::rank_me() == 0) {
    //     hashmap.print_local();
    // }

    // ==============================================
    
    upcxx::global_ptr<int> x = upcxx::new_<int>(42);

    // if (upcxx::rank_me() == 0){
    //     x = 69;
    // }

    // if (upcxx::rank_me() == 0){
    //     std::cout << "X0 " << x << std::endl;
    // }

    // if (upcxx::rank_me() == 1){
    //     std::cout << "X1 " << x << std::endl;
    // }

    upcxx::delete_(x);
    
    upcxx::global_ptr<int> x_arr = upcxx::new_array<int>(10);
    UPCXX_ASSERT(x_arr.is_local());   // a precondition of global_ptr<T>::local()
    int *local_ptr = x_arr.local();
    local_ptr[1] = 4;

    if (upcxx::rank_me() == 0){
        std::cout << "Xarr1 " << local_ptr[1] << std::endl;
    }

    // ==============================================
    upcxx::finalize();
    return 0;
}
