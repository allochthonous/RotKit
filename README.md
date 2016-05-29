## What is RotKit?
RotKit is a python-based set of tools and routines for manipulating Euler rotations and using them for reconstructing plate motions and related tectonic parameters.

## Current Status
As of May 2016, RotKit is still in the very early stages of development; although many of the core routines are working, stability of code (or said code’s behaviour) is not guaranteed!

## Use of Fortran
Some functions are written in Fortran and have been linked to python using f2py. For these functions to work on a new machine, RotKit.f and RotKit-nonpython.f need to be compiled. The instructions listed in RotKit_compile.txt are currently very specific to my own set-up (compiled using gfortran on a Mac running OSX 10.10) but may provide some guidance.