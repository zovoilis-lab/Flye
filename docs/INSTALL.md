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
* Python 2.7
* Git
* Core OS development headers (zlib, etc)


Get the latest source version
-----------------------------

To get and compile the latest git version, run:

    git clone https://github.com/fenderglass/Flye
	cd Flye
    python setup.py build


After building, Flye could be invoked with the following command:

    bin/flye

Installation (optional)
-----------------------
You may install the package for the better OS integration:

    python setup.py install

Alternatively, you can perform local user installation by adding ```--user``` or ```--prefix```
options to the previous command.
