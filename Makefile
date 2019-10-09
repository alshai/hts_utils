CXX=g++
CXX_FLAGS=-std=c++11

all: score_sam merge_sams

score_sam: score_sam.cpp hts_dict.hpp bam_util.hpp
	$(CXX) -o $@ $< -lhts

merge_sams: merge_sams.cpp hts_dict.hpp bam_util.hpp
	$(CXX) -o $@ $< -lhts
