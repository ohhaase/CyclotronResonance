# This file is being included in ./CMakeLists.txt

# Create library for pcg (header only)
add_library(pcgLib INTERFACE)

# Add the pcg folder to the include directories
target_include_directories(pcgLib INTERFACE pcg)