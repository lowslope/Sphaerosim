cmake_minimum_required(VERSION 3.10)

if (ScoreP)
    set(CMAKE_C_COMPILER "scorep-$ENV{CC}")
    set(CMAKE_CXX_COMPILER "scorep-$ENV{CXX}")
    set(LIBRARIES "${LIBRARIES}" "-L/work/xu077052/lib/")
else (ScoreP)
    set(CMAKE_C_COMPILER "C:/msys64/mingw64/bin/cc.exe")#"$ENV{CC}")
    set(CMAKE_CXX_COMPILER "C:/msys64/mingw64/bin/c++.exe")#"$ENV{CXX}")
endif (ScoreP)

set(CMAKE_CXX_STANDARD 17)
project(sphaerosim)

set(CMAKE_CXX_FLAGS "$ENV{FLAGS_OPENMP} -fopenmp -static -I/opt/intel/vtune/XE2019-u02/vtune_amplifier_xe/include")
#set(CMAKE_CXX_FLAGS_DEBUG "-O1 -ip -xSSE4.2 -axCORE-AVX2,AVX -fp-model fast=2 -g3 -debug")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -inline-level=0")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{FLAGS_FAST} $ENV{TEST_FLAGS_COMPILER}")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    if (THREAD_SANITIZER)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w0")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -inline-forceinline -xHost") #-inline-forceinline
    set(LIBRARIES "${LIBRARIES}" "-L/rwthfs/rz/SW/gcc/CENTOS-7.3/8.2.0/lib64")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif ()

#set(LIBRARIES "-L/opt/intel/vtune/XE2019-u02/vtune_amplifier_xe/lib64/" "-littnotify")


if (LIKWID_PERFMON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLIKWID_PERFMON -I/opt/likwid/include/")
    set(LIBRARIES "${LIBRARIES}" "-L/opt/likwid/lib -llikwid")
endif ()

if (INSPECTOR)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINSPECTOR")
endif ()

if (Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif ()

add_subdirectory(src)
include_directories(${INCLUDE})
add_executable(Sphaerosim "${SOURCES}")
target_link_libraries(Sphaerosim "${LIBRARIES}")