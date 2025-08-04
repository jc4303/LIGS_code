cmake_minimum_required (VERSION 3.6)
project(hnsw_lib
    LANGUAGES CXX)
# Ensure CMake 3.12 or newer is used to avoid compatibility issues.
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    message(FATAL_ERROR "CMake 3.12 or higher is required.")
endif()

# Find Python 3 and its components
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

add_subdirectory(external/pybind11)
include_directories(external/pybind11/include)

find_package(MKL REQUIRED)
link_directories(${MKL_LIBRARY_DIRS})

pybind11_add_module(python_interface bin_methods_double.cpp)
target_link_libraries(python_interface PRIVATE ligslib)


add_library(ligslib INTERFACE)
target_include_directories(ligslib INTERFACE .) 

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 17)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      SET( CMAKE_CXX_FLAGS  "-Ofast -DNDEBUG -lstdc++fs -std=c++17 -DHAVE_CXX0X -openmp -march=native -fpic -ftree-vectorize")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      SET( CMAKE_CXX_FLAGS  "-Ofast -lrt -DNDEBUG -lstdc++fs -std=c++17 -DHAVE_CXX0X -march=native -fpic -w -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=0" )
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
      SET( CMAKE_CXX_FLAGS  "-Ofast -lrt -DNDEBUG -lstdc++fs -std=c++17 -DHAVE_CXX0X -openmp -march=native -fpic -w -fopenmp -ftree-vectorize" )
    endif()


endif()
