# shasta2

**Development. Do not use.**

## Clone and build

* Requires gcc 15. Known to work on Ubuntu 26.04.
* Assumes `python3`, `git`, and `cmake` are available.

```
git clone https://github.com/paoloshasta/shasta2.git
shasta2/scripts/InstallBuildPrerequisites.py
mkdir shasta2-build
cd shasta2-build
cmake ../shasta2
make all
make install
```

* The call to `InstallBuildPrerequisites.py` requires `sudo` access and will ask for a password if needed.
* It installs a few required packages using `apt_get` and installs the `abpoa` and `theseus` libraries in `~/.shasta2Build`.

This creates:
* Static executable `shasta2-build/shasta2-install/bin/shasta2`. It is highly portable to other 64-bit `x86_64` Linux systems. It can be moved to other systems and immediately used without requiring installation.
* Python3 module `shasta2-build/shasta2-install/bin/shasta2.so`. This is not portable to other Linux distributions.

