Posture Generator
=================
[![Build Status](https://travis-ci.org/jorisv/PG.svg?branch=multi_robot)](https://travis-ci.org/jorisv/PG)

![LGPL-3](https://www.gnu.org/graphics/lgplv3-88x31.png)

Posture Generator (PG) is a library that computes the inverse kinematics of
kinematic trees.


## Dependencies

* C++ compiler with C++11 support
* Boost
* CMake
* Eigen 3
* RBDyn
* SpaceVecAlg
* roboptim-core
* sch-core

For Python bindings:

* Python 2 or 3
* PyBindGen


## Install

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=[your install prefix]
make
make install
```

Note: the usual install prefix is `/usr`, but it defaults to `/usr/local`.

If you generate Python bindings on Debian-based distributions (e.g. Ubuntu),
you may want to choose the Debian layout for Python packages (e.g.
`/usr/lib/python2.7/dist-packages/`) instead of the default one (e.g.
`/usr/lib/python2.7/site-packages/`). If this is the case, add the
`-DPYTHON_DEB_LAYOUT=ON` flag to the CMake command.

