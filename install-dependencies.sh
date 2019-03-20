#!/bin/bash

# Upgrade pip. We're using Python 2.7 since rlogin doesn't have pip3.
# We /could/ use Anaconda for Python 3, but that's probably asking too much.
pip --disable-pip-version-check install --user -U pip

# Now use the freshly installed version of pip.
~/.local/bin/pip --disable-pip-version-check \
                 install --user -U networkx numpy scipy matplotlib pandas

