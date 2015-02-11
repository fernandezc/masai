#!/bin/bash

# tell the script to exit on errors
set -e

# get the correct python environment
source /Users/christian/Library/Enthought/Canopy_64bit/User/bin/activate

# go to masai tests directory, and test
cd ~/Dropbox/PycharmProjects/masai/tests
py.test -x

# go to the tutorial directry and test ipynb
#cd ~/Dropbox/PycharmProjects/masai/doc/tutorial
#python test_tutorial.py

# return to masai directory
cd ~/Dropbox/PycharmProjects/masai
pwd

# commit
git add .
git commit -a -m "Test passed"
