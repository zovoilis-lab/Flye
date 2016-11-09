ABruijn Installation
====================

Availability
------------

ABruijn is available for Linux and MacOS platforms. Windows support is not guaranteed.


Build requirements
------------------

* C++ compiler with C++11 support (GCC 4.8+ / Clang 3.3+ / Apple Clang 5.0+)
* GNU make
* Python 2.7


Runtime Requirements
--------------------

* BLASR aligner [https://github.com/PacificBiosciences/blasr]


Installation
------------

ABruijn is written in both C++ and Python and requires a Python 2.7 environment.

First, to build ABruijn, run:

    python install.py build

ABruijn could be invoked with the following command:

    bin/abruijn

Additonally, you may install the package for the better OS integration:

    python setup.pu install

Alternatively, you can perform local user installation by adding '--user' or '--prefix'
options to the previous command.

BLASR installation
------------------

ABruijn also requires the BLASR aligner to be installed and available through the PATH environment variable.
For example if you have the 'blasr' binary at '/aaa/bbb/blasr', you can type the following command
to make it visible for ABruijn:

    export $PATH=$PATH:/aaa/bbb/blasr

You can test if the binary is available by running:

    which blasr

you should see '/aaa/bbb/blasr' as output.
