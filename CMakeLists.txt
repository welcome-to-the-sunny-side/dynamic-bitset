set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 3.10)
project(DBITSET)

set(CMAKE_CXX_STANDARD 17)

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