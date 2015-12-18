cmake_minimum_required(VERSION 3.1)
project(hdr_src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    COPYING
    fix_chromatic_aberration.cc
    fix_chromatic_aberration.hh
    fix_noise.cc
    fix_noise.hh
    hdr_adjust.cc
    hdr_convert.cc
    hdr_create.cc
    hdr_denoise.cc
    hdr_fix_ca.cc
    hdr_squish.cc
    headers.hh
    helpers.hh
    image.cc
    image.hh
    Makefile
    shared_array.hh
    solve_poisson.cc
    solve_poisson.hh
    squish_luminance.cc
    squish_luminance.hh)

add_executable(hdr_src ${SOURCE_FILES})