cmake_minimum_required(VERSION 3.16)
project(c3pOSEM)

set(CMAKE_CXX_STANDARD 14)

ADD_COMPILE_OPTIONS(-Wl,-Bsymbolic -fvisibility=hidden)

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

file(GLOB OSEM_SRC "*.h" "*.cpp")

add_library(c3pOSEM SHARED ${OSEM_SRC})