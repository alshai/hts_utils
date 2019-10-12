#include <cstdio>
#include <getopt.h>
#include "hts_dict.hpp"
#include "bam_util.hpp"

struct Args {
    int buffer = 0;
};

struct Aln {
    Aln(bam_hdr_t* hdr, bam1_t* rec) :
       r(hdr->target_name[rec->core.tid]),
       p(rec->core.pos) {}
    std::string r;
    int32_t p;
};

void print_help() {
    fprintf(stderr, "usage: ./score_sam.cpp truth.[sb]am test.[sb]am\n");
    exit(1);
}

/* given a set of "true" alignments and a set of test alignments,
 * for each record in the test alignments, output:
 * READ_NAME TRUE_REF TRUE_POS TEST_REF TEST_POS CORRECT
 * where "CORRECT" is 0 if incorrect, 1 if correct
 * By correct, we mean that TEST_REF==TRUE_REF and TEST_POS==TRUE_POS
 *
 * usage: ./score_sam truth.[bs]am test.[bs]am
 */
int main(int argc, char** argv) {
    Args args;
    static struct option long_options[] {
        {"buffer", required_argument, 0, 'b'},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "hb:", long_options, NULL)) != -1 ) {
        switch(ch) {
            case 'b':
                args.buffer = std::atoi(optarg);
                break;
            case 'h':
            default:
                print_help();
                exit(1);
        }
    }

    samFile* fps[2];
    bam_hdr_t* hdrs[2];
    if (argc - optind < 2) {
        print_help();
    }
    for (int i = 0; i < 2; ++i) {
        fps[i] = sam_open(argv[optind++], "r");
        hdrs[i] = sam_hdr_read(fps[i]);
    }
    auto truth_set = hts_util::bam_to_dict<Aln>(fps[0], hdrs[0]);
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
