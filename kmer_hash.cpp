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

using namespace std;


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
    
    // if (run_type == "verbose") {
    //     BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
    //         n_kmers);
    //     }
        
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());
    
    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }
    
    // Load factor of 0.5
    // size_t hash_table_size = n_kmers * (1.0 / 0.5);
    size_t hash_table_size = kmers.size() * (1.0 / 0.5);

    // std::cout << "Rank " << upcxx::rank_me() << " kemrs size " << kmers 

    // DistributedHashMap hashmap(local_size);
    DistributedHashMap hashmap(hash_table_size);

    upcxx::barrier();

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;

    if (upcxx::rank_me() == 0) {
        std::cout << "Size " << hashmap.size() << std::endl;
    }

    upcxx::barrier();
    std::vector<upcxx::future<>> futures;

    for (auto& kmer : kmers) {
        futures.push_back(hashmap.insert(kmer));

        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
            std::cout << "rank " << upcxx::rank_me() << " kmer " << kmer.kmer_str() << " " << kmer.backwardExt() << std::endl;
        }
    }

    upcxx::when_all(futures.begin(), futures.end()).wait();

    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf\n", insert_time);
    }
    upcxx::barrier();

    std::cout << "Rank " << upcxx::rank_me() << " hash map size " << hashmap.size() << std::endl;

    // ===================== READ ===================

    auto start_read = std::chrono::high_resolution_clock::now();

    std::list<std::list<kmer_pair>> contigs;
    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair kmer = hashmap.find(contig.back().next_kmer()).wait();
            // bool success = hashmap.find(contig.back().next_kmer(), kmer);
            // if (!success) {
            //     throw std::runtime_error("Error: k-mer not found in hashmap.");
            // }
            contig.push_back(kmer);
        }
        contigs.push_back(contig);
    }

    if (upcxx::rank_me() == 1) {
        const auto& first_contig = contigs.front();
        int count = 0;
        for (const auto& kmer : first_contig) {
            std::cout << "kmer " << count << ": " << kmer.kmer_str() << " " << kmer.backwardExt() << kmer.forwardExt() << std::endl;
            if (++count >= 10) break;
        }
    }

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    auto end = std::chrono::high_resolution_clock::now();

    upcxx::finalize();
    return 0;
}
