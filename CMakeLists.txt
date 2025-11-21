cmake_minimum_required(VERSION 3.12)
project(VectorMesonProduction)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------------------------------------------------------
# Locate GSL
# -----------------------------------------------------------------------------
find_package(GSL REQUIRED)

# -----------------------------------------------------------------------------
# Locate Cuba
# -----------------------------------------------------------------------------
find_library(CUBA_LIBRARY NAMES cuba)
find_path(CUBA_INCLUDE_DIR cuba.h)

if(NOT CUBA_LIBRARY OR NOT CUBA_INCLUDE_DIR)
    message(FATAL_ERROR "Cuba library not found. Please install it or set CUBA_LIBRARY and CUBA_INCLUDE_DIR manually.")
endif()

# -----------------------------------------------------------------------------
# Include ROOT
# -----------------------------------------------------------------------------
find_package(ROOT REQUIRED)

# Use ROOT's compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


# -----------------------------------------------------------------------------
# Include directories
# -----------------------------------------------------------------------------
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${GSL_INCLUDE_DIRS}
    ${CUBA_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

# -----------------------------------------------------------------------------
# Source files
# -----------------------------------------------------------------------------
file(GLOB SOURCES
    src/*.cpp
    src/*.cxx
    src/*.C
)

# -----------------------------------------------------------------------------
# Executable
# -----------------------------------------------------------------------------
add_executable(VMProduction ${SOURCES})

# -----------------------------------------------------------------------------
# Link libraries
# -----------------------------------------------------------------------------
target_link_libraries(VMProduction
    ${Python3_LIBRARIES}
    ${GSL_LIBRARIES}
    ${CUBA_LIBRARY}
    ${ROOT_LIBRARIES}  

    m
    dl
    crypt
)

# -----------------------------------------------------------------------------
# Compiler options
# -----------------------------------------------------------------------------
target_compile_options(VMProduction PRIVATE -Wall -Wextra -O2)
