CXX=g++
CXX_FLAGS=-std=c++11

all: score_sam

score_sam: score_sam.cpp hts_dict.hpp bam_util.hpp
	$(CXX) -o $@ $< -lhts

