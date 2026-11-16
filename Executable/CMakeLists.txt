cmake_minimum_required(VERSION 3.16)
project(shasta2Executable)

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

# Source files
file(GLOB SOURCES ../src/*.cpp)

# Include directories.
include_directories(../src)
include_directories(~/.shasta2Build/abpoa/abPOA/include)
include_directories(~/.shasta2Build/theseus/theseus-lib/include)

# Link directories for abpoa, poasta, and theseus.
link_directories(~/.shasta2Build/abpoa/abPOA/lib)
link_directories(~/.shasta2Build/theseus/theseus-lib/staticBuild)

# Include and link directories for plegma, if requested.
# plegma is a static library, libplegma.a.
if(BUILD_WITH_PLEGMA)
    add_definitions(-DSHASTA2_PLEGMA)
    include_directories(~/.shasta2Build/plegma/include)
    link_directories(~/.shasta2Build/plegma/lib)
    set(SHASTA2_PLEGMA_LIBRARY plegma)
endif(BUILD_WITH_PLEGMA)

# Define our executable.
add_executable(shasta2Executable ${SOURCES})
set_target_properties(shasta2Executable PROPERTIES OUTPUT_NAME "shasta2")

# Request a static executable.
set_target_properties(shasta2Executable PROPERTIES LINK_FLAGS "-static" )

# Libraries to link with.
target_link_libraries(
    shasta2Executable
    atomic boost_program_options boost_chrono boost_serialization png z
    abpoa theseus ${SHASTA2_PLEGMA_LIBRARY} z m pthread)

  
# The static executable goes to the bin directory.
install(TARGETS shasta2Executable DESTINATION shasta2-install/bin)


