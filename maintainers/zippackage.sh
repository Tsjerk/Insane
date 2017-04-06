#!/bin/bash
set -euo pipefail

# Build a zip package of insane
#
# In order to be easily used from the command line without installation, while
# being modular for developers and library users, we distribute insane an
# executable zip. That zip can be called directly from the command line (e.g.
# `./insane -h`) or through python (e.g. `python insane -h`).
#
# The distributed zip is a zipped version of the repository prepended with a
# shebang.

# Path from which to download simopt
SIMOPT_URL='https://raw.githubusercontent.com/Tsjerk/simopt/master/simopt.py'

# Zip is annoying when it comes to path inside the archive.
# We make a copy of the repository that we can modify without worrying to much.
mkdir -p build
rsync -ruv ../* build \
    --exclude='tests' \
    --exclude='maintainers' \
    --exclude='*.egg-info' \
    --exclude='setup.py'
cd build

# Download simopt
wget "${SIMOPT_URL}" -O simopt.py || curl -O "${SIMOPT_URL}" 

# __main__.py is the entry point of the script. We do not want it in the
# module as it makes no sense to have it there, so we write it here.
cat << EOF > __main__.py
import os
import sys
import insane.cli

# When running the zipped package, __file__ is set to insane.__main__.py
# Hide the internals to the user; he does not need to know what is inside the
# zip. If anything else than simopt uses __main__.__file__ it may be broken.
__file__ = os.path.dirname(__file__)

sys.exit(insane.cli.main(sys.argv))
EOF

# Actually make the zip file. Finally!
zip -r insane.zip *

# Make it executable
echo '#!/usr/bin/env python' | cat - insane.zip > package
chmod +x package

# Clean the mess
cd ../
mv build/package ./insane
rm -rf build

# Were are done. Let's say good bye.
echo 'The zip package is created. Try `./insane -h`.'
