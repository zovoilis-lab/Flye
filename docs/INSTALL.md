ABruijn Installation
====================

Availability
------------

ABuijn is available for Linux and MacOS platforms. Windows support is not guaranteed.


Build requirements
------------------

* C++ compiler with C++11 support (GCC 4.8+ / Clang 3.3+ / Apple Clang 5.0+)
* GNU make


Runtime Requirements
--------------------

* Python 2.7
* BLASR aligner [https://github.com/PacificBiosciences/blasr]


Installation
------------
ABruijn is written in both C++ and Python and requires Python 2.7 environment.

First, to build native C++ modules type

    make

ABruijn also requires the BLASR aligner to be installed and available through PATH environment variable.
For example if you have 'blasr' binary at '/aaa/bbb/blasr', you can type the following command
to make it visible for ABruijn:

    export $PATH=$PATH:/aaa/bbb/blasr

You can test if the binary is available by running:

    which blasr

in ABruijn directory - you should see '/aaa/bbb/blasr' as output.
