cmake_minimum_required(VERSION 3.14)
project(NTT_GPU_ACCEL LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenFHE REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenFHE_CXX_FLAGS}")

include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
include_directories(${OpenFHE_INCLUDE}/binfhe)
include_directories(${OpenFHE_INCLUDE}/third-party/include)

link_directories(${OpenFHE_LIBDIR})

add_executable(cpu-baseline cpu-baseline.cpp)

target_link_libraries(cpu-naive PRIVATE ${OpenFHE_SHARED_LIBRARIES})
target_link_libraries(cpu-rns PRIVATE ${OpenFHE_SHARED_LIBRARIES})

target_link_libraries(cpu-baseline PRIVATE ${OpenFHE_SHARED_LIBRARIES})
target_compile_options(cpu-baseline PRIVATE -O3)

option(USE_BARRETT "Enable Barrett reduction for modular arithmetic" OFF)
if(USE_BARRETT)
    target_compile_definitions(cpu-baseline PRIVATE USE_BARRETT)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cpu-baseline PRIVATE OpenMP::OpenMP_CXX)
endif()