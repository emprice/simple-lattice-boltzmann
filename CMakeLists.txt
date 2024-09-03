cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0079 NEW)
cmake_minimum_required(VERSION 3.28)

project(simplelb VERSION 1.0 LANGUAGES C CXX)

find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

add_executable(simplelb simplelb.cc)
target_link_libraries(simplelb PUBLIC hdf5::hdf5 OpenMP::OpenMP_CXX)

# vim: set ft=cmake:
