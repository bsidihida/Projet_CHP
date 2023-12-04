NEOS
====

[![pipeline status](https://gitlab.inria.fr/memphis/neos/badges/develop/pipeline.svg)](https://gitlab.inria.fr/memphis/neos/commits/develop)
[![coverage report](https://gitlab.inria.fr/memphis/neos/badges/develop/coverage.svg)](https://gitlab.inria.fr/memphis/neos/commits/develop)

Developper Documentation
------------------------

[Documentation is Here](https://memphis.gitlabpages.inria.fr/neos/)

Get Neos
--------
To use last development state of Neos, please clone the master branch.

To get sources please use these commands:

      git clone git@gitlab.inria.fr:memphis/neos.git
      cd neos
      git submodule init
      git submodule update


My first code with Neos
-----------------------

Neos is available as a Docker Container.

- On your Operating System, install [Docker Desktop](https://docs.docker.com/get-docker/).

- Launch Docker Desktop

- In a Terminal Shell get Neos Docker Image:

```
$ docker pull registry.gitlab.inria.fr/memphis/neos
```

- Create a working directory (for example, the directory will be in your home dir, in a subfolder called Neos):

```
$ mkdir -p ~/Neos/first_example
```

- Get into Neos, bind your working directory with a similar one in Docker:

```
$ docker run -v~/Neos/first_example:/builds/first_example -it registry.gitlab.inria.fr/memphis/neos /bin/bash
To run a command as administrator (user "root"), use "sudo <command>".
See "man sudo_root" for details.

gitlab@8cc432b5546a:/builds$ ls
neos first_example
```

- Go in your working directory, copy template files

```
gitlab@8cc432b5546a:/builds$ cd first_example
gitlab@8cc432b5546a:/builds/first_example$
gitlab@8cc432b5546a:/builds/first_example$ cp ../neos/first_example/* .
gitlab@8cc432b5546a:/builds/first_example$ ls
CMakeLists.txt  levelsept.cpp
```

- Compile your code in CMake way:

```
gitlab@8cc432b5546a:/builds/first_example$ mkdir build && cd build
gitlab@8cc432b5546a:/builds/first_example$ cmake .. && make
...
```

- Run it

```
gitlab@8cc432b5546a:/builds/first_example$ ./levelset
...
```

Build and install
-----------------
Neos requires to have some library dependencies:
   * Lapacke or Eigen3 (Eigen is enable by default)
   * MPI (optional)
   * PetsC (optional: for laplacian)
   * Bitpit (optional but strongly recommended for CORE and GEOMETRY features) (http://optimad.github.io/bitpit/)
   * mpi4py and pybind11 pip modules (optional: for Python bindings)

The main options to configure Neos build are:
   * BUILD_CORE: Enable levelset, transport...
   * BUILD_GEOMETRY: Enable somme standard geometries
   * BUILD_MATH:
      * BUILD_INTERPOLATOR: Enable distance weighted, Polynomial and RBF interpolation.
      * BUILD_LAPLACIAN: Enable laplacian

# Installation example

```
# install bitpit
cd $HOME
git clone https://github.com/optimad/bitpit.git
cd bitpit
export BITPIT_ROOT=$PWD
mkdir build && cd build
cmake .. -DPETSC_DIR=/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/ -DPETSC_ARCH='' -DBUILD_DOCUMENTATION=ON -DBUILD_EXAMPLES=ON -DENABLE_MPI=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$BITPIT_ROOT

cd $HOME
git clone --recursive git@gitlab.inria.fr:memphis/neos.git
cd neos
mkdir build ; cd build
export CMK=" -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD../install"
export CMK2=" -DBUILD_CORE=ON -DBUILD_GEOMETRY=ON -DBUILD_LAPLACIAN=ON -DENABLE_MPI=ON -DGEN_PYTHON=ON -DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON"
export CMK3=" -DPETSC_DIR=/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/ -DPETSC_ARCH=''"
cmake ..  $CMK $CMK1 $CMK2 $CMK3 -DBUILD_DOCUMENTATION=ON
make -j4

examples/time-poisson
```

## To customize examples

with your preferred editor, edit `../examples/time-poisson.cpp`

```
make
examples/time-poisson
```

## To install Python module

```
cd python
pip  install --user --upgrade ./pyneos
# or
pip3 install --user --upgrade ./pyneos
```

# CI Tools

  * CI/VM: https://ci.inria.fr/project/neos/show
  * SonarQube: https://sonarqube.inria.fr/sonarqube/dashboard?id=memphis%3Aneos%3Agit%3Adevelop
