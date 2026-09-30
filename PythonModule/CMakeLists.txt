cmake_minimum_required(VERSION 3.16)

project(shasta2PythonModule)

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

# This makes sure the code in PythonModule.cpp is compiled.
add_definitions(-DSHASTA2_PYTHON_API)

# Sources files.
file(GLOB SOURCES ../src/*.cpp)

# Include directories.
include_directories(../src)
include_directories(~/.shasta2Build/abpoa/abPOA/include)
include_directories(~/.shasta2Build/theseus/theseus-lib/include)

# Link directories for abpoa and poasta.
link_directories(~/.shasta2Build/abpoa/abPOA/lib)
link_directories(~/.shasta2Build/theseus/theseus-lib/sharedBuild)

# Define our library.
add_library(shasta2PythonModule SHARED ${SOURCES})

# Make sure the library is named shasta.so.
set_target_properties(shasta2PythonModule PROPERTIES OUTPUT_NAME "shasta2")
set_target_properties(shasta2PythonModule PROPERTIES PREFIX "")
set_target_properties(shasta2PythonModule PROPERTIES DEFINE_SYMBOL "")
set_target_properties(shasta2PythonModule PROPERTIES INSTALL_RPATH "$ORIGIN")



# Include directories and libraries for Python3 and pybind11.
execute_process(COMMAND /usr/bin/python3-config --embed --libs OUTPUT_VARIABLE SHASTA2_PYTHON_LIBRARIES)
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE SHASTA2_PYTHON_INCLUDES)
add_definitions(${SHASTA2_PYTHON_INCLUDES})
string(STRIP "${SHASTA2_PYTHON_LIBRARIES}" SHASTA2_PYTHON_LIBRARIES)
SET(CMAKE_LINKER_FLAGS  "${CMAKE_LINKER_FLAGS} ${SHASTA2_PYTHON_LIBRARIES}")



# Libraries to link with.
target_link_libraries(
    shasta2PythonModule 
    atomic png boost_program_options boost_serialization boost_chrono pthread 
    abpoa theseus z m ${SHASTA2_PYTHON_LIBRARIES})

# Install the shared library and the abpoa, and theseus libraries it needs into the bin directory.
install(TARGETS shasta2PythonModule DESTINATION shasta2-install/bin)
install(FILES ~/.shasta2Build/abpoa/abPOA/lib/libabpoa.so DESTINATION shasta2-install/bin)
install(FILES ~/.shasta2Build/theseus/theseus-lib/sharedBuild/libtheseus.so DESTINATION shasta2-install/bin)




