# Specify CMake version(s)
cmake_minimum_required(VERSION 3.15...4.3)

# Declare project
project(
    CyclotronResonance
    VERSION 1.0
    LANGUAGES CXX
)

# Compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default build type (adds -O3)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Declare the executable
add_executable(${PROJECT_NAME})

# Add source files
target_sources(${PROJECT_NAME} PRIVATE
    # ./src/
    src/main.cpp
    src/global_vars.cpp

    # ./src/core/
    src/core/helper_functions.cpp
    src/core/parallel_functions.cpp
    src/core/sim_types.cpp

    # ./src/helper_objects/
    src/helper_objects/ElectronDistb.cpp
    src/helper_objects/Histogram.cpp
    src/helper_objects/Histogram2D.cpp
    src/helper_objects/NumericalDistb.cpp
    src/helper_objects/Timer.cpp

    # ./src/science_functions/
    src/science_functions/beta_calculations.cpp
    src/science_functions/cross_sections.cpp
    src/science_functions/sampling_functions.cpp
    src/science_functions/scattering_functions.cpp
)

# Add include subdirectories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

# Link external libraries
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
