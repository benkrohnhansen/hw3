#include <upcxx/upcxx.hpp>
#include <iostream>
#include <vector>
#include <iterator>
#include "hash_map.hpp"
#include "read_kmers.hpp"    // Provides: std::vector<kmer_pair> read_kmers(const char* filename)
#include "kmer_t.hpp"        // Defines kmer_pair and pkmer_t
#include "packing.hpp"       // Defines KMER_LEN

#include <cstdint>
#include <string>

//------------------------------------------------------------------------------
// Helper: convert a pkmer_t to a std::string by decoding its bits.
// Assumes:
//   - pkmer_t is stored as an unsigned 64-bit integer (or convertible to uint64_t)
//   - Each nucleotide is encoded in 2 bits: 0->'A', 1->'C', 2->'G', 3->'T'.
//   - KMER_LEN (from packing.hpp) gives the length of the k-mer.
//------------------------------------------------------------------------------
std::string pkmer_to_string(const pkmer_t &pk) {
    int len = KMER_LEN; // KMER_LEN should be defined in packing.hpp.
    std::string s(len, 'A');
    // Assume pkmer_t is convertible to uint64_t.
    uint64_t value = static_cast<uint64_t>(pk);
    for (int i = len - 1; i >= 0; i--) {
        int bits = value & 0x3; // extract 2 bits
        char nucleotide;
        switch(bits) {
            case 0: nucleotide = 'A'; break;
            case 1: nucleotide = 'C'; break;
            case 2: nucleotide = 'G'; break;
            case 3: nucleotide = 'T'; break;
            default: nucleotide = 'N'; break;
        }
        s[i] = nucleotide;
        value >>= 2;
    }
    return s;
}

//------------------------------------------------------------------------------
// Overloaded operator<< for kmer_pair.
// We assume that kmer_pair contains a member 'kmer' of type pkmer_t and
// that its getter functions for the extensions are forwardExt() and backwardExt().
//------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const kmer_pair &kp) {
    // Convert pkmer_t to string.
    std::string kmer_str = pkmer_to_string(kp.kmer);
    os << kmer_str;
    os << " (F:" << kp.forwardExt() << " B:" << kp.backwardExt() << ")";
    return os;
}

//------------------------------------------------------------------------------
// Helper function: construct a kmer_pair from a bucket.
// The repository provides a constructor for kmer_pair of the form:
//    kmer_pair(const std::string&, const std::string&)
// where the first argument is the DNA string and the second is a string
// containing the two extension characters concatenated.
//------------------------------------------------------------------------------
kmer_pair bucket_to_kmer(const bucket_t &b) {
    std::string kmer_str = pkmer_to_string(b.key);
    std::string ext_str;
    ext_str.push_back(b.forwardExt);
    ext_str.push_back(b.backwardExt);
    return kmer_pair(kmer_str, ext_str);
}

int main(int argc, char** argv) {
    upcxx::init();
    int rank = upcxx::rank_me();
  
    if(argc < 2) {
        if(rank == 0)
            std::cerr << "Usage: " << argv[0] << " <kmer_dataset>" << std::endl;
        upcxx::finalize();
        return 1;
    }
  
    const char* dataset = argv[1];
  
    // Choose an appropriate table size (example: 1M buckets).
    size_t table_size = 1 << 20;
  
    // Create the distributed hash table.
    hash_map dht(table_size);
  
    // Read the k-mer dataset.
    std::vector<kmer_pair> kmers = read_kmers(dataset);
  
    // Insert each k-mer into the distributed hash table.
    for (const auto &k : kmers) {
        bucket_t bk;
        bk.key = k.kmer;                      // k.kmer is of type pkmer_t.
        bk.forwardExt = k.forwardExt();       // Call the getter.
        bk.backwardExt = k.backwardExt();
        bk.valid = true;
        bool inserted = dht.insert(bk);
        (void)inserted; // Optionally handle duplicates.
    }
  
    upcxx::barrier();  // Ensure all insertions are complete.
  
    //------------------------------------------------------------------------------
    // Contig Generation:
    // For each k-mer that starts a contig (i.e. k.backwardExt() == 'F'),
    // traverse the de Bruijn graph (using next_kmer()) until reaching a k-mer
    // with forwardExt() == 'F'.
    //------------------------------------------------------------------------------
    std::vector< std::vector<kmer_pair> > contigs;
    for (const auto &k : kmers) {
        if (k.backwardExt() == 'F') {  // Contig start.
            std::vector<kmer_pair> contig;
            contig.push_back(k);
            kmer_pair current = k;
      
            while (current.forwardExt() != 'F') {
                // Compute the next k-mer key using the provided member function.
                pkmer_t next_key = current.next_kmer();
        
                bucket_t found_bucket;
                bool found = dht.find(next_key, found_bucket);
                if (!found)
                    break;  // End traversal if not found.
        
                // Construct the next kmer_pair from the found bucket using our helper.
                kmer_pair next_k = bucket_to_kmer(found_bucket);
        
                contig.push_back(next_k);
                current = next_k;
            }
            contigs.push_back(contig);
        }
    }
  
    // For demonstration, print each contig.
    for (const auto &contig : contigs) {
        std::cout << "Rank " << rank << " contig: ";
        for (const auto &k : contig) {
            std::cout << k << " ";  // Uses the overloaded operator<< for kmer_pair.
        }
        std::cout << std::endl;
    }
  
    upcxx::finalize();
    return 0;
}
