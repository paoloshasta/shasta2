#!/usr/bin/python3

"""
This clones the plegma repository (https://github.com/kokyriakidis/plegma)
at a fixed commit and builds the plegma static library,
used by local assembly method Plegma.

plegma is written in Mojo. The build uses pixi (https://pixi.sh),
which downloads the Mojo compiler it needs. If pixi is not already
installed, this installs it into ~/.shasta2Build/plegma/pixi
without changing the shell startup files.

The static library includes the few Mojo runtime functions plegma uses,
so shasta2 needs no Mojo library at run time and the executable stays static.
On x86-64 it is built for x86-64-v3 (AVX2), so it runs on machines
other than the build machine.

Output goes to ~/.shasta2Build/plegma
Include files in ~/.shasta2Build/plegma/include
Library in ~/.shasta2Build/plegma/lib/libplegma.a
"""

import os
import shutil
import subprocess

# The plegma commit to build.
plegmaCommit = "4a6aa8f5dcabae0a0014552ee5edd4188edf4579"

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

# Create the plegma directory, if it does not already exist.
if os.path.exists("plegma"):
    raise Exception(rootDirectory + "/" + "plegma already exists.")
else:
    os.mkdir("plegma")

# cd to the plegma directory.
os.chdir("plegma")
plegmaDirectory = os.getcwd()

# Find pixi, or install it if necessary.
pixi = shutil.which("pixi")
if pixi is None and os.path.exists(os.path.expanduser("~/.pixi/bin/pixi")):
    pixi = os.path.expanduser("~/.pixi/bin/pixi")
if pixi is None:
    print("Installing pixi in " + plegmaDirectory + "/pixi")
    environment = dict(os.environ)
    environment["PIXI_HOME"] = plegmaDirectory + "/pixi"
    environment["PIXI_NO_PATH_UPDATE"] = "1"
    subprocess.run("curl -fsSL https://pixi.sh/install.sh | bash", shell=True, check=True, env=environment)
    pixi = plegmaDirectory + "/pixi/bin/pixi"
print("Using " + pixi)

# Clone the plegma repository and its submodules at the chosen commit.
subprocess.run(["git", "clone", "https://github.com/kokyriakidis/plegma.git"], check=True)
os.chdir("plegma")
subprocess.run(["git", "checkout", plegmaCommit], check=True)
subprocess.run(["git", "submodule", "update", "--init", "--recursive"], check=True)

# Build the static library, build/libplegma.a.
print("Building the plegma static library.")
subprocess.run([pixi, "run", "--locked", "build-static"], check=True)

# Copy the include files and the library.
os.mkdir("../include")
os.mkdir("../lib")
for name in ["plegma.h", "plegma.hpp"]:
    shutil.copy("capi/" + name, "../include/")
shutil.copy("build/libplegma.a", "../lib/")
print("plegma include files are in " + plegmaDirectory + "/include")
print("plegma library is in " + plegmaDirectory + "/lib")
