cmake_minimum_required(VERSION 3.14)
project(NTT_GPU_ACCEL LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenFHE REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenFHE_CXX_FLAGS}")

include_directories(include)
include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
include_directories(${OpenFHE_INCLUDE}/binfhe)
include_directories(${OpenFHE_INCLUDE}/third-party/include)

link_directories(${OpenFHE_LIBDIR})

set(OPENFHE_OMP_SOURCES
    src/ntt_cpu.cpp
    src/ntt_openfhe.cpp
    bench/main.cpp
)

add_executable(bench_ntt ${OPENFHE_OMP_SOURCES})
target_link_libraries(bench_ntt PRIVATE ${OpenFHE_SHARED_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_options(bench_ntt PRIVATE -O3)
target_compile_definitions(bench_ntt PRIVATE BENCH_NTT)


option(USE_BARRETT "Enable Barrett reduction" OFF)
if(USE_BARRETT)
    # target_compile_definitions(cpu-baseline PRIVATE USE_BARRETT)
    target_compile_definitions(bench_ntt PRIVATE USE_BARRETT)
endif()
