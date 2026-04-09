cmake_minimum_required(VERSION 3.14)
project(NTT_GPU_ACCEL LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_ARCHITECTURES "70;80")

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
    src/ntt_gpu.cu
    bench/main.cpp
    src/bench_utils.cpp    
    src/arg_parser.cpp
)

add_executable(bench_ntt ${OPENFHE_OMP_SOURCES})
target_link_libraries(bench_ntt PRIVATE ${OpenFHE_SHARED_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_options(bench_ntt PRIVATE -O3)
target_compile_definitions(bench_ntt PRIVATE BENCH_NTT)



