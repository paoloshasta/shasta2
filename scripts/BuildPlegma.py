#!/usr/bin/python3

"""
This clones the plegma repository (https://github.com/kokyriakidis/plegma)
at a fixed commit and builds the plegma shared library,
used by local assembly method Plegma.

plegma is written in Mojo. The build uses pixi (https://pixi.sh),
which downloads the Mojo compiler it needs. If pixi is not already
installed, this installs it into ~/.shasta2Build/plegma/pixi
without changing the shell startup files.

plegma is only available as a shared library, because the Mojo runtime
it uses is only available as shared libraries. So a shasta2 build
with plegma (cmake option -DBUILD_WITH_PLEGMA=ON) is not static.

Output goes to ~/.shasta2Build/plegma
Include files in ~/.shasta2Build/plegma/include
Libraries in ~/.shasta2Build/plegma/lib: libplegma.so and the Mojo
runtime libraries it needs, which it finds in its own directory.
"""

import os
import shutil
import subprocess

# The plegma commit to build.
plegmaCommit = "021b2c7f5c4a5b978f3e9c1b484f436528d2329e"

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

# Build the shared library. This puts libplegma.so and the
# Mojo runtime libraries it needs in the build directory.
print("Building the plegma shared library.")
subprocess.run([pixi, "run", "--locked", "build-capi"], check=True)

# Copy the include files and the libraries.
os.mkdir("../include")
os.mkdir("../lib")
for name in ["plegma.h", "plegma.hpp"]:
    shutil.copy("capi/" + name, "../include/")
for name in os.listdir("build"):
    if name.endswith(".so"):
        shutil.copy("build/" + name, "../lib/")
print("plegma include files are in " + plegmaDirectory + "/include")
print("plegma libraries are in " + plegmaDirectory + "/lib")
