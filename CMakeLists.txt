cmake_minimum_required(VERSION 3.10)

project(varcount CXX C)

set(PROJECT_URL "https://github.com/alshai/varcount")
set(PROJECT_DESCRIPTION "Outputs alignment coverage over a predefined set of variants")
set(CMAKE_CXX_STANDARD 11)

# htslib
find_package(PkgConfig REQUIRED)
pkg_search_module(PKG_HTSLIB REQUIRED htslib)

add_executable(score_sam score_sam.cpp bam_util.hpp)
target_link_libraries(score_sam ${PKG_HTSLIB_LDFLAGS})

add_executable(merge_sams merge_sams.cpp hts_dict.hpp bam_util.hpp)
target_link_libraries(merge_sams ${PKG_HTSLIB_LDFLAGS})

add_executable(gt_dist gt_dist.cpp hts_dict.hpp bcf_util.hpp)
target_link_libraries(gt_dist ${PKG_HTSLIB_LDFLAGS})

# INSTALL(TARGETS score_sam merge_sams gt_dist DESTINATION bin)
