#!/bin/bash

# tell the script to exit on errors
set -e

# get the correct python environment
source /Users/christian/Library/Enthought/Canopy_64bit/User/bin/activate

# go to masai tests directory, and test
cd ~/Dropbox/PycharmProjects/masai/tests
py.test -x

# go to the tutorial directry and test ipynb   # deactivated as runipy seems not to work with my Plot1d routine. TODO
# cd ~/Dropbox/PycharmProjects/masai/doc/tutorial
# python test_tutorial.py_Deactivated

# return to masai directory
cd ~/Dropbox/PycharmProjects/masai
pwd

# commit
git add .
git commit -a -m "Test passed"
