set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 3.10)
project(DBITSET)

message(STATUS "Using C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "Using CXX compiler: ${CMAKE_CXX_COMPILER}")

set(CMAKE_CXX_STANDARD 23)

# add_compile_options(-O3 -march=native -funroll-loops)
# add_compile_options(-O0)
add_compile_options(-O3 -march=native)

# if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#     add_compile_options(-fopt-info-vec-all)
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     add_compile_options(-Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize)
# endif()

# Include project headers
include_directories(${CMAKE_SOURCE_DIR}/include)

# Include your headers if needed
include_directories(${CMAKE_SOURCE_DIR})

# Create ORAM library interface
add_library(dyna_bitset INTERFACE)
target_include_directories(dyna_bitset INTERFACE 
    ${CMAKE_SOURCE_DIR}/include
)

# Find Google Benchmark
find_package(benchmark REQUIRED)

# target_link_libraries(dyna_bitset INTERFACE
# )

# === Main Executable ===
add_executable(main
    src/main.cpp
)

target_link_libraries(main PRIVATE dyna_bitset benchmark::benchmark)

# === Test Executables ===
add_executable(shift_comp
    test/shift_comp.cpp
)

target_link_libraries(shift_comp PRIVATE dyna_bitset benchmark::benchmark)

add_executable(bitwise_comp
    test/bitwise_comp.cpp
)

target_link_libraries(bitwise_comp PRIVATE dyna_bitset benchmark::benchmark)

add_executable(mixed_comp
    test/mixed_comp.cpp
)

target_link_libraries(mixed_comp PRIVATE dyna_bitset benchmark::benchmark)