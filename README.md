HClib
=============================================

HClib is a task-based parallel programming model that supports the finish-async,
parallel-for, and future-promise parallel programming patterns through both C
and C++ APIs. HClib explicitly
exposes hardware locality of the hardware, while allowing the programmer to fall
back on sane defaults. The HClib runtime is a lightweight, work-stealing, and
locality-aware runtime. HClib is not itself an exascale programming system, but
is intended to be the intra-node resource management and scheduling component
within an exascale programming system, integrating with inter-node communication
models such as MPI, UPC++, or OpenSHMEM.

[![Build Status](https://travis-ci.org/habanero-rice/hclib.svg?branch=master)](https://travis-ci.org/habanero-rice/hclib)

Installation
---------------------------------------------

HClib follows your standard bootstrap, configure, and make installation
procedure. At its simplest, the manual installation process consists of:

    ./bootstrap.sh
    ./configure --prefix=<installation-dir>
    make install

However, an install.sh script is also provided for your convenience that will
build and install HClib. Simply run the script to install:

    ./install.sh

In addition, the install script will auto-generate a file that can be used
to configure your environment (more on that later).

By default, HClib will be installed by the installation script to
`$PWD/hclib-install`. If you would like to use
a different installation location, you can override the `INSTALL_PREFIX`
environment variable:

    INSTALL_PREFIX=/opt/local ./install.sh

Likewise, if you would like to use different C/C++ compilers other than the
system defaults, then you can specify them using the `CC` and `CXX` environment
variables. For example, if you want to use the Intel compilers:

    CC=icc CXX=icpc ./install.sh

After installation, you will need to set the `HCLIB_ROOT` environment variable
to point to your
HClib installation directory. You can automatically set this variable after
installation by sourcing the `hclib_setup_env.sh` script generated by install.sh. For example, assuming
HClib was installed with `INSTALL_PREFIX=/opt/local`:

    source /opt/local/bin/hclib_setup_env.sh

HClib Modules
---------------------------------------------

To support unified scheduling on heterogeneous computation and communication
resources HClib uses a module system. Loading a given module automatically adds
support the the HClib API and runtime for accessing some new resource
in HClib. By default, HClib only supports executing parallel
programs on multi-core x86 platforms. However, with existing modules this can be
extended to include support for GPU execution and communication over MPI,
OpenSHMEM, and UPC++.

While you can write a fully functioning HClib program without any additional
modules, many of the test programs saved in this repo load the `system` module.
The `system` module adds basic OS-related routines, such as asynchronous memory
allocation and deallocation. If you would like to run any of these basic tests,
the instructions below guide you through the process of building and installing
the 'ssytem' module.

If install.sh is used to create an HClib installation, the `system` module will
automatically be built and installed. However, if you wish to configure and
install HClib manually you will also need to build and install the `system`
module manually. Once you have completed your HClib install, navigate to the
`hclib/modules/system` directory and run:

    make install

Ensure that you have `HCLIB_ROOT` set in your environment first.

Dependencies
---------------------------------------------

* automake
* gcc >= 4.8.4, or clang >= 3.5
  (must support -std=c++11 and -std=c11)

Tutorial
---------------------------------------------

If you are new to HClib then take a look of `hclib/tutorial` directory.
It contains presentations and simple examples that appeared in our
past tutorials on HClib. You can follow the README inside sub-directories
there to build and run those examples.


Testing
---------------------------------------------

The main regression tests for HClib are in the test/c and test/cpp folders. The
`test_all.sh` scripts in each of those folders will automatically build and run
all test cases.


Static Checks
---------------------------------------------

As part of the development workflow for HClib, any newly committed code should
be checked using standard static checking tools.

In particular, run cppcheck on all modified files. cppcheck is available online
at [1]. cppcheck should be run by cd-ing to tools/cppcheck and executing the
run.sh script from there (this assumes cppcheck is on your path). Any new errors
printed by cppcheck should be addressed before committing.

You should also run astyle on all modified files. astyle is a source code
auto-formatter. Simply cd to tools/astyle and execute the run.sh script from
there. This assumes you have astyle installed and it is on your path.

[1] https://sourceforge.net/projects/cppcheck/
