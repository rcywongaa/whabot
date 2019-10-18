#!/bin/bash

sudo apt install python3-dill python3-sympy

# See
# https://drake.mit.edu/from_binary.html
# https://drake.mit.edu/python_bindings.html#python-bindings-binary
wget https://drake-packages.csail.mit.edu/drake/nightly/drake-latest-bionic.tar.gz
sudo tar -xvzf drake-latest-bionic.tar.gz -C /opt
sudo /opt/drake/share/drake/setup/install_prereqs

echo "export PYTHONPATH=/opt/drake/lib/python3.6/site-packages:${PYTHONPATH}" >> ~/.bashrc
