#ifndef BAM_UTIL_HPP
#define BAM_UTIL_HPP

namespace hts_util {
    /* very basic alignment struct. use as template and customize for your own
     * purposes
     */
    struct Aln {
        Aln(bam_hdr_t* h, bam1_t* r) {
            hdr = h;
            rec = bam_dup1(r);
        }
        bam_hdr_t* hdr;
        bam1_t* rec;
    };
};

#endif
