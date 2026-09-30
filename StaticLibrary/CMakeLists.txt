cmake_minimum_required(VERSION 3.16)

project(shasta2StaticLibrary)

# C++ dialect.
add_definitions(-std=c++26)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter -Wno-unused-function)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
else(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# This makes sure main is excluded.
add_definitions(-DSHASTA2_STATIC_LIBRARY)

# Sources files.
file(GLOB SOURCES ../src/*.cpp)

# Include directories.
include_directories(../src ~/.shasta2Build/abpoa/abPOA/include)
include_directories(~/.shasta2Build/theseus/theseus-lib/include)

# Define our library.
add_library(shasta2StaticLibrary STATIC ${SOURCES})

# Make sure the library is named shasta2.a.
set_target_properties(shasta2StaticLibrary PROPERTIES OUTPUT_NAME "shasta2")
set_target_properties(shasta2StaticLibrary PROPERTIES PREFIX "")

# Install the static library into the bin directory.
install(TARGETS shasta2StaticLibrary DESTINATION shasta2-install/bin)




