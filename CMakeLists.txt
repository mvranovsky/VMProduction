cmake_minimum_required(VERSION 3.12)
project(MyPlotApp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------------------------------------------------------
# Locate Python and NumPy
# -----------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

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
# Include directories
# -----------------------------------------------------------------------------
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${CUBA_INCLUDE_DIR}
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
    m
    dl
    crypt
)

# -----------------------------------------------------------------------------
# Compiler options
# -----------------------------------------------------------------------------
target_compile_options(VMProduction PRIVATE -Wall -Wextra -O2)
