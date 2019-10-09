#ifndef HTS_DICT_HPP
#define HTS_DICT_HPP

#include <cinttypes>
#include <cstdio>
#include <string>
#include <htslib/sam.h>
#include "flat_hash_map.hpp"

namespace hts_util {
    template <typename T>
    using Int2T = ska::flat_hash_map< int32_t, T >;

    template <typename T, template<typename> typename U>
    using Str2Map = ska::flat_hash_map< std::string, U<T>>;

    template <typename T>
    using Str2T = ska::flat_hash_map< std::string, T >;

    using StrSet = ska::flat_hash_set< std::string >;

    /* T must have constructor that takes bam1_t* as sole argument */
    template<typename T>
    Str2T<T> bam_to_dict(samFile* fp, bam_hdr_t* hdr) {
        Str2T<T> read_dict;
        bam1_t* rec = bam_init1();
        while (sam_read1(fp, hdr, rec) >= 0) {
            std::string read_name(bam_get_qname(rec));
            if (read_dict.find(read_name) != read_dict.end()) {
                fprintf(stderr, "read already in dict, skipping...\n");
            } else {
                read_dict.insert_or_assign(read_name, T(hdr, rec)); 
            }
        }
        return read_dict;
    }
};

#endif // HTS_DICT_HPP
