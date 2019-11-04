cmake_minimum_required(VERSION 2.8.12)
project(ConstraintForce)

find_package(drake CONFIG REQUIRED)
find_package(pybind11 REQUIRED)

include_directories(include)

pybind11_add_module(ConstraintForce ConstraintForce_py.cpp ConstraintForce.cpp)
#add_library(ConstraintForce ConstraintForce.cpp)
target_link_libraries(ConstraintForce
  drake::drake drake::drake-common-text-logging-gflags gflags_shared
)
