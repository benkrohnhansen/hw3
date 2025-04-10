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
    HashMap hashmap(hash_table_size);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
                     n_kmers);
    }

    // ==============================================
    
    upcxx::global_ptr<int> x = upcxx::new_<int>(42);

    // if (upcxx::rank_me() == 0){
    //     x = 69;
    // }

    if (upcxx::rank_me() == 0){
        std::cout << "X0 " << x << std::endl;
    }

    if (upcxx::rank_me() == 1){
        std::cout << "X1 " << x << std::endl;
    }

    upcxx::delete_(x);

    // upcxx::global_ptr<int> x_arr = upcxx::new_array<int>(10);

    // upcxx::delete_array(x_arr);

    // ==============================================
    upcxx::finalize();
    return 0;
}
