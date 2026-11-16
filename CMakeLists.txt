# Use scripts/InstallPrerequisites.py to install required build dependencies.

cmake_minimum_required(VERSION 3.16)
project(shasta2)

# Summarize the environment.
message(STATUS "Processor architecture is " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "CMAKE_SYSTEM_NAME is " ${CMAKE_SYSTEM_NAME})

# Decide what we want to build.
option(BUILD_EXECUTABLE "Build the executable." ON)
option(BUILD_PYTHON_MODULE "Build the Python module." ON)
option(BUILD_STATIC_LIBRARY "Build the static library." OFF)

# Write out what we will actually build.
message(STATUS "BUILD_EXECUTABLE is " ${BUILD_EXECUTABLE})
message(STATUS "BUILD_PYTHON_MODULE is " ${BUILD_PYTHON_MODULE})
message(STATUS "BUILD_STATIC_LIBRARY is " ${BUILD_STATIC_LIBRARY})

# Option to request a debug build.
option(BUILD_DEBUG "Make a debuggable build." OFF)
message(STATUS "BUILD_DEBUG is " ${BUILD_DEBUG})

# Option to build with plegma, for local assembly method Plegma.
# Requires running scripts/BuildPlegma.py first.
# plegma is only available as a shared library,
# so with this option the executable is not static.
option(BUILD_WITH_PLEGMA "Build with plegma." OFF)
message(STATUS "BUILD_WITH_PLEGMA is " ${BUILD_WITH_PLEGMA})
if(BUILD_WITH_PLEGMA AND NOT EXISTS $ENV{HOME}/.shasta2Build/plegma/lib/libplegma.so)
    message(FATAL_ERROR "BUILD_WITH_PLEGMA requires ~/.shasta2Build/plegma/lib/libplegma.so. "
        "Run scripts/BuildPlegma.py to create it.")
endif(BUILD_WITH_PLEGMA AND NOT EXISTS $ENV{HOME}/.shasta2Build/plegma/lib/libplegma.so)

# Add the subdirectories we need.
if(BUILD_EXECUTABLE)
    add_subdirectory(Executable)
endif(BUILD_EXECUTABLE)

if(BUILD_PYTHON_MODULE)
    add_subdirectory(PythonModule)
endif(BUILD_PYTHON_MODULE)

if(BUILD_STATIC_LIBRARY)
    add_subdirectory(StaticLibrary)
endif(BUILD_STATIC_LIBRARY)



# Install to the shasta-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION shasta2-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION shasta2-install)





