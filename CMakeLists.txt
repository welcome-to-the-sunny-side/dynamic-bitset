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

# target_link_libraries(dyna_bitset INTERFACE
# )

# === Test Executable ===
add_executable(test
    src/test.cpp
)

target_link_libraries(test PRIVATE dyna_bitset)