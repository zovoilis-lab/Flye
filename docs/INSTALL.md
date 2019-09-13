Flye Installation
=================

Availability
------------

Flye is available for Linux and MacOS platforms.

Bioconda Releases
-----------------

You can get the latest stable release through Bioconda:

    conda install flye

Alternatively, you can get a release version from the github releases page


Building Requirements
---------------------

* C++ compiler with C++11 support (GCC 4.8+ / Clang 3.3+ / Apple Clang 5.0+)
* GNU make
* Python 2.7 or 3.4+ with setuptools package installed
* Git
* Core OS development headers (zlib, etc)


Local building (without installation)
-------------------------------------

You may use the package locally without system installation.
To get and compile the latest git version, run:

    git clone https://github.com/fenderglass/Flye
    cd Flye
    make

Then, Flye will be available as:

    python3 bin/flye

This example is using Python 3, replace ```python3```
with ```python2``` to use the 2nd version. 


Installing from source
----------------------

To install the Flye package into your system, run:

    git clone https://github.com/fenderglass/Flye
	cd Flye
    python3 setup.py install

Depending on your OS, you might need to add
```--user``` or ```--prefix``` options to the 
install command for the local installation.

This example is using Python 3, replace ```python3```
with ```python2``` to use the 2nd version. 
After installation, Flye could be invoked via:

    flye
