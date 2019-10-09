#include <cstdio>
#include "hts_dict.hpp"
#include "bam_util.hpp"

void die(const char* s) {
    fprintf(stderr, s);
    exit(1);
}

/* returns "better" alignment. 1 if l >= r, 0 otherwise.
 * looks by AS:i field first
 * then by mapq score
 * if mapqs AND alignment scores are equal, return 1
 *
 * usage: ./merge_sams sam1.[bs]am sam2.[bsam]
 */
bool compare_bam1_t(bam1_t* l, bam1_t* r) {
    // l->core.flag & 4 -> 0
    // r->core.flag & 4 -> 1
    // l->core.flag & 4 && r->core.flag & 4 -> 1
    if ((l->core.flag & 4) || (r->core.flag & 4)) {
        return r->core.flag & 4;
    } else {
        uint8_t* laux = bam_aux_get(l, "AS");
        uint8_t* raux = bam_aux_get(r, "AS");
        int las = laux ? bam_aux2i(laux) : -1;
        int ras = raux ? bam_aux2i(raux) : -1;
        return las == ras ? l->core.qual >= r->core.qual : las > ras;
    }
}

/* merges two SAM/BAM files. If two records for the same reads exist, then report
 * the "better".
 * here, we define "better" as the alignment that has:
 * 1) is aligned
 * 2) OR has the better alignment score
 * 3) OR has the better MAPQ score
 * in the future, we will add options for customizing how to define "better"
 * */
int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "please specify truth sam file AND test sam file\n");
        exit(1);
    }
    samFile* fps[2];
    bam_hdr_t* hdrs[2];
    for (int i = 0; i < 2; ++i) {
        fps[i] = sam_open(argv[i+1], "r");
        hdrs[i] = sam_hdr_read(fps[i]);
    }
    // open output sam file
    samFile* outfp = sam_open("-", "w");
    /* TODO: 
     * merge the headers properly
     * add PG field
     */ 
    bam_hdr_t* outhdr = sam_hdr_dup(hdrs[0]); 
    if (sam_hdr_write(outfp, outhdr) != 0) {
        die("error writing header\n");
    }

    auto sam1_d = hts_util::bam_to_dict<hts_util::Aln>(fps[0], hdrs[0]);
    hts_util::StrSet sam2_s;
    bam1_t* b = bam_init1();
    int i = 0;
    while (sam_read1(fps[1], hdrs[1], b) >= 0) {
        const char* qname(bam_get_qname(b));
        sam2_s.emplace(qname);
        const char* rname = hdrs[1]->target_name[b->core.tid];
        int32_t pos = b->core.pos;
        auto pair = sam1_d.find(std::string(qname));
        if (pair != sam1_d.end()) {
            auto a = pair->second.rec;
            int j = compare_bam1_t(a,b);
            if (sam_write1(outfp, outhdr, compare_bam1_t(a,b) ? a : b) < 0) {
                die("error writing to sam file\n");
            }
        } else {
            if (sam_write1(outfp,outhdr, b) < 0) {
                die("error writing to sam file\n");
            }
        }
        ++i;
    }
    /* find all recs in sam1 that were not present in sam2 */
    for (auto a: sam1_d) {
        if (sam2_s.find(a.first) == sam2_s.end()) {
            if (sam_write1(outfp,outhdr,a.second.rec) < 0) { 
                die("error writing to sam file\n");
            }
        }
    }

    for (int i = 0; i < 2; ++i) {
        sam_hdr_destroy(hdrs[i]);
        sam_close(fps[i]);
    }
    sam_close(outfp);
    sam_hdr_destroy(outhdr);
    bam_destroy1(b);
}
