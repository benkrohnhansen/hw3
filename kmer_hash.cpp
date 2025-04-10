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
    
    int rank = upcxx::rank_me();
    int n_ranks = upcxx::rank_n();

    // Each rank holds an integer in a dist_object
    upcxx::dist_object<int> dist_val(rank);

    // Ensure all ranks initialized their dist_object
    upcxx::barrier();

    // Rank 0 will do the following:
    // - rget value from rank 1 (if it exists)
    // - rput value to rank 1 (if it exists)
    // - rpc to rank 1 (if it exists)
    if (rank == 0 && n_ranks > 1) {
        // rget from rank 1
        upcxx::future<int> fut_val = upcxx::rget(dist_val.fetch(1));
        int val = fut_val.wait();

        std::cout << "Rank 0: rget from rank 1 got value = " << val << std::endl;

        // rput to rank 1 (send value 42)
        upcxx::rput(42, dist_val.fetch(1)).wait();

        std::cout << "Rank 0: rput 42 to rank 1" << std::endl;

        // rpc to rank 1
        upcxx::rpc(1, []() {
            std::cout << "Rank 1 (via rpc): Hello from rank 1!" << std::endl;
        }).wait();
    }

    upcxx::barrier();

    // Final output from rank 0 to confirm execution
    if (rank == 0) {
        std::cout << "Rank 0: Finished UPC++ demo." << std::endl;
    }
    // ==============================================

    upcxx::finalize();
    return 0;
}
