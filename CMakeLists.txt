# Specify CMake version(s)
cmake_minimum_required(VERSION 3.15...4.3)

# Declare project
project(
    CyclotronResonance
    VERSION 1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Declare the final executable
add_executable(${PROJECT_NAME})

# Add source subdirectory
add_subdirectory(src)

# Add include subdirectory
add_subdirectory(include)

# Add Libraries from the include directory
target_link_libraries(${PROJECT_NAME} PRIVATE pcgLib)

# Add openMP
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
