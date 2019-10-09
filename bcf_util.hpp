#ifndef VCF_UTIL_HPP
#define VCF_UTIL_HPP

#include <cstdio>
#include <cmath>
#include <array>
#include <htslib/vcf.h>
#include "flat_hash_map.hpp"

namespace hts_util {
    /* very basic variant struct. use as template and customize for your own
     * purposes */
    struct Var {
        public:
        Var(bcf_hdr_t* h, bcf1_t* r) {
            hdr = h;
            rec = bcf_dup(r);
        }
        bcf_hdr_t* hdr;
        bcf1_t* rec;
    }

    /*** GENOTYPING FUNCTIONS ***/
    int32_t* get_genotype(bcf_hdr_t* hdr, bcf1_t* rec, int* ngt) {
        int32_t *gt_arr = NULL, ngt_arr = 0;
        *ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if ( ngt<=0 ) return NULL; // GT not present
        else return gt_arr;
    }

    /* input: ref allele read count, alt allele read count, uniform error rate
     * output: genotype likelihood for ref/ref, ref/alt, alt/alt genotypes, resp.
     * assumes uniform base qualities. Makes no distinction between SNPS and indels
     */
    inline std::array<int, 3> pl_naive(int rc, int ac, float e) {
        std::array<int, 3> pls;
        pls[0] = static_cast<int>( -10 * ((rc*std::log10(1-e)) + (ac*std::log(e))) );
        pls[1] = static_cast<int>( 10 * (rc + ac) * std::log10(2) );
        pls[2] = static_cast<int>( -10 * ((ac*std::log10(1-e)) + (rc*std::log(e))) );
        return pls;
    }

    inline std::array<int32_t, 3> pl_naive_normalized(int32_t rc, int32_t ac, float e) {
        std::array<int32_t, 3> pls;
        pls[0] = static_cast<int32_t>( -10 * ((rc*std::log10(1-e)) + (ac*std::log(e))) );
        pls[1] = static_cast<int32_t>( 10 * (rc + ac) * std::log10(2) );
        pls[2] = static_cast<int32_t>( -10 * ((ac*std::log10(1-e)) + (rc*std::log(e))) );
        int32_t min = pls[0];
        for (int i = 1; i < 3; ++i) {
            if (pls[i] < min) min = pls[i];
        }
        for (int i = 0; i < 3; ++i) {
            pls[i] -= min;
        }
        return pls;
    }

    /* assuming diploid gts, so len(gt_arr) = nsmpl * 2
    * returns: {RR count, RA count, AA count}
    */
    inline std::array<int32_t, 3> gt_freq(int32_t* gt_arr, int32_t nsmpl) {
        std::array<int32_t, 3> gf = {0, 0, 0};
        for (int i = 0; i < nsmpl; ++i) {
            int gt1 = 2 * i, gt2 = (2*i)+1;
            gf[0] += !gt1 & (gt1 == gt2); 
            gf[1] += (gt1 != gt2);
            gf[2] += gt1 & (gt1 == gt2);
        }
        return gf;
    }
};

#endif // VCF_UTIL_HPP
