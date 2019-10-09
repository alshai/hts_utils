#include <cstdio>
#include "hts_dict.hpp"
#include "bam_util.hpp"

inline void report_read(FILE* fp, const char* qname, const char* r1, int p1, const char* r2, int p2, int c) {
    fprintf(fp,"%s\t%s\t%d\t%s\t%d\t%d\n",qname,r1,p1,r2,p2,c);
}

struct Aln {
    Aln(bam_hdr_t* hdr, bam1_t* rec) :
       r(hdr->target_name[rec->core.tid]),
       p(rec->core.pos) {}
    std::string r;
    int32_t p;
};

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "please specify truth sam file AND test sam file\n");
        exit(1);
    }
    samFile* fps[2];
    bam_hdr_t* hdrs[2];
    //fprintf(stderr, "opening files\n");
    for (int i = 0; i < 2; ++i) {
        fps[i] = sam_open(argv[i+1], "r");
        hdrs[i] = sam_hdr_read(fps[i]);
    }
    //fprintf(stderr, "storing truth in dict\n");
    auto truth_set = hts_util::bam_to_dict<Aln>(fps[0], hdrs[0]);
    for (auto x: truth_set) {
        fprintf(stdout, "%s\n", x.first.data());
    }
    bam1_t* rec = bam_init1();
    const char* fmt_str = "%s\t%s\t%d\t%s\t%d\t%d\n";
    while (sam_read1(fps[1], hdrs[1], rec) >= 0) {
        const char* qname(bam_get_qname(rec));
        const char* rname = hdrs[1]->target_name[rec->core.tid];
        int32_t pos = rec->core.pos;
        auto pair = truth_set.find(std::string(qname));
        if (pair != truth_set.end()) {
            auto a = pair->second;
            fprintf(stdout,fmt_str,qname,a.r.data(),a.p,rname,pos,a.p==pos);
        } else {
            fprintf(stdout,fmt_str,qname,"*",0,rname,pos,0);
        }
    }
}
