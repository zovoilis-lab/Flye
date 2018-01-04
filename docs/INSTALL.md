Flye Installation
=================

Availability
------------

Flye is available for Linux and MacOS platforms.

Requirements
------------

* C++ compiler with C++11 support (GCC 4.8+ / Clang 3.3+ / Apple Clang 5.0+)
* GNU make
* Python 2.7
* Git
* Basic OS development headers (zlib, etc.)

Get the latest version (recommended)
------------------------------------

To get and compile the latest git version, run:

    git clone https://github.com/fenderglass/Flye
	cd Flye
    python setup.py build

Alternatively, you can get a release verson from the "releases" page.

After building, Flye could be invoked with the following command:

    bin/flye

Installation (optional)
-----------------------
You may install the package for the better OS integration:

    python setup.py install

Alternatively, you can perform local user installation by adding ```--user``` or ```--prefix```
options to the previous command.
