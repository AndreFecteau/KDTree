cmake_minimum_required(VERSION 3.4)

# set the project name
project(KDTree)

set(CMAKE_CXX_STANDARD 17)

add_executable(tests tests/test_KDTree.cpp)
target_include_directories(tests PUBLIC "lib/catch/include")