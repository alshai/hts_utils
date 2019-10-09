CXX=g++
CXX_FLAGS=-std=c++11

all: score_sams

score_sams: score_sam.cpp hts_dict.hpp bam_util.hpp
	$(CXX) -o $@ $< -lhts

