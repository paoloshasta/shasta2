cmake_minimum_required(VERSION 3.16)
project(shasta2Executable)

# C++ dialect.
add_definitions(-std=c++23)

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

# Include directories for abpoa and poasta..
include_directories(../src ~/.shasta2Build/abpoa/abPOA/include)
include_directories(~/.shasta2Build/poasta-c/poasta-c)

# Link directories for abpoa and poasta.
link_directories(~/.shasta2Build/abpoa/abPOA/lib)
link_directories(~/.shasta2Build/poasta-c/poasta-c/target/release)

# Define our executable.
add_executable(shasta2Executable ${SOURCES})
set_target_properties(shasta2Executable PROPERTIES OUTPUT_NAME "shasta2")

# Request a static executable.
set_target_properties(shasta2Executable PROPERTIES LINK_FLAGS "-static" )

# Libraries to link with.
target_link_libraries(
    shasta2Executable
    atomic boost_program_options boost_chrono boost_serialization png z
    abpoa poasta_c z m pthread)

  
# The static executable goes to the bin directory.
install(TARGETS shasta2Executable DESTINATION shasta2-install/bin)


