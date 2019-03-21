#!/bin/bash

# Upgrade pip. We're using the Python 3 installation just made.
pip3 --disable-pip-version-check install --user -U pip

# Now install packages with updated pip3.
pip3 install --user -U networkx numpy scipy matplotlib pandas

echo "Versions installed:"
python3 -V
pip3 freeze | egrep 'networkx|numpy|scipy|matplotlib|pandas'
