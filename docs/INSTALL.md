ABruijn installation
====================

Availability
------------

ABuijn is available for Linux and MacOS platforms.


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
ABruijn was written in both C++ and Python and requires Python 2.7 environment.

First, to build native C++ modules type

    make

ABruijn also requires the BLASR aligner to be installed (available through PATH environment variable).
For example if you have 'blasr' binary at '/aaa/bbb/blasr', you may type the following command
to make it available for ABruijn:

    export $PATH=$PATH:/aaa/bbb/blasr

You can test if the binary available by running:

    which blasr

in the ABruijn directory - it should output '/aaa/bbb/blasr'.
