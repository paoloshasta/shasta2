#!/usr/bin/python3

"""
This clones the albertjimenezbl/theseus-lib
Github repository and builds the theseus static and shared libraries. 

All output goes to ~/.shasta2Build/theseus
The static build goes to ~/.shasta2Build/theseus/theseus-lib/staticBuild
The shared build goes to ~/.shasta2Build/theseus/theseus-lib/sharedBuild
"""

import os

# Get the path to the root directory.
rootDirectoryWithTilda = "~/.shasta2Build"
rootDirectory = os.path.expanduser(rootDirectoryWithTilda)

# Create the root directory, if it does not already exists.
if os.path.exists(rootDirectory):
    if not os.path.isdir(rootDirectory):
        raise Exception(rootDirectory + " exists but is not a directory")
else:
    os.mkdir(rootDirectory)
    
# cd to the root directory.
os.chdir(rootDirectory)

# Create the theseus directory, if it does not already exist.
if os.path.exists("theseus"):
    raise Exception(rootDirectory + "/" + "theseus already exists.")
else:
    os.mkdir("theseus")

# cd to the theseus directory.
os.chdir("theseus")

# Clone the albertjimenezbl/theseus-lib repository.
os.system("git clone https://github.com/albertjimenezbl/theseus-lib.git")

# cd to the theseus-lib directory.
os.chdir("theseus-lib")

# Build the static library.
print("Building the theseus static library.")
os.mkdir("staticBuild")
os.chdir("staticBuild")
os.system("cmake .. -DCMAKE_BUILD_TYPE=Release")
os.system("make")
os.chdir("..")

# Build the shared library
print("Building the theseus shared library.")
os.mkdir("sharedBuild")
os.chdir("sharedBuild")
os.system("cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON")
os.system("make")

