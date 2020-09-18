cmake_minimum_required(VERSION 3.17)
project(c3pOSEM)

set(CMAKE_CXX_STANDARD 14)

file(GLOB OSEM_SRC "*.h" "*.cpp")

add_library(c3pOSEM SHARED ${OSEM_SRC})