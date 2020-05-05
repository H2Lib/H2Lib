H2Lib
=====

Welcome to H2Lib.

If you have questions, please feel free to contact the authors.
Most header files give the name of the author in charge, and you
can find the corresponding e-mail address in the file "AUTHORS".

Enjoy.

[![Build Status](https://travis-ci.org/d-meiser/H2Lib.png?branch=master)](https://travis-ci.org/d-meiser/H2Lib)
[![Coverage Status](https://coveralls.io/repos/d-meiser/H2Lib/badge.svg?branch=master)](https://coveralls.io/r/d-meiser/H2Lib?branch=master)


## Obtaining the public version of H2Lib

The public version of H2Lib can be obtained from the github repository
(https://github.com/H2Lib/H2Lib):

    git clone htts://github.com/H2Lib/H2Lib

This will download the H2Lib library into a directory called `H2Lib` in
the current working directory.


## Building 

To get started, you just need

  * a C build environment (compiler, linker, libraries) and
  * a version of the "make" tool.

We recommend that you also obtain optimized implementations of
BLAS and LAPACK, since these improve the performance of H2Lib
significantly.

The traditional build system consists of gnu make files that can be
customized.  Alternatively, there is a cmake based build system.
Configuration on Windows may be easier with cmake.


### gnu make

The build process is controlled by the files "system.inc" and
"options.inc".
The first file describes the properties of your system, e.g.,
which compiler to use and where certain libraries can be found.
The second file describes build options, e.g., whether you want
to use real or complex numbers.
The repository includes examples of both files for standard Linux
systems.

Once the configuration files are prepared, just enter

    make

to build the library in "Library" and the test programs in
"Tests".

If you have Doxygen installed, you can use

    make doc

to create the documentation of the library.


### cmake

To configure H2Lib with cmake do the following:

    mkdir build
    cd build
    cmake ..
    make

The test suite can be run using (still in the `build` directory)

    make test

The documentation can be build with

    make doc

Note that the documentation is currently still being built in the source
directory.
