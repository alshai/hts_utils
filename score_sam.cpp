#include <cstdio>
#include "hts_dict.hpp"
#include "bam_util.hpp"

struct Aln {
    Aln(bam_hdr_t* hdr, bam1_t* rec) :
       r(hdr->target_name[rec->core.tid]),
       p(rec->core.pos) {}
    std::string r;
    int32_t p;
};

/* given a set of "true" alignments and a set of test alignments,
 * for each record in the test alignments, output:
 * READ_NAME TRUE_REF TRUE_POS TEST_REF TEST_POS CORRECT
 * where "CORRECT" is 0 if incorrect, 1 if correct
 * By correct, we mean that TEST_REF==TRUE_REF and TEST_POS==TRUE_POS
 *
 * usage: ./score_sam truth.[bs]am test.[bs]am
 */
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
    auto truth_set = hts_util::bam_to_dict<Aln>(fps[0], hdrs[0]);
    bam1_t* rec = bam_init1();
    // fprintf(stderr, "read truth\n");
    int i = 0;
    const char* fmt_str = "%s\t%s\t%d\t%s\t%d\t%d\n";
    while (sam_read1(fps[1], hdrs[1], rec) >= 0) {
        const char* qname(bam_get_qname(rec));
        // fprintf(stderr, "%s starting\n", qname);
        bool aligned = !(rec->core.flag & 4);
        const char* rname = aligned ? hdrs[1]->target_name[rec->core.tid] : "NA";
        int32_t pos = aligned ? rec->core.pos : -1;
        auto pair = truth_set.find(std::string(qname));
        if (pair != truth_set.end()) {
            auto a = pair->second;
            fprintf(stdout,fmt_str,qname,a.r.data(),a.p,rname,pos,a.p==pos);
        } else {
            fprintf(stdout,fmt_str,qname,"*",0,rname,pos,0);
        }
        // fprintf(stderr, "%s done\n", qname);
        ++i;
    }
}
