# 🧬🧰 PyNCBItk [![Stars](https://img.shields.io/github/stars/althonos/pyncbitk.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyncbitk/stargazers)

*(Unofficial) [Cython](https://cython.org/) bindings and Python interface to the [NCBI C++ Toolkit](https://www.ncbi.nlm.nih.gov/toolkit).*

This package contains the runtime components of PyNCBItk to allow reusing them
across versions of the Cython bindings without having to rebuild. See the main
[`pyncbitk`](https://pypi.org/project/pyncbitk) package for more information.


## 🔧 Installing

The PyNCBItk runtime supports pre-built binaries for Linux x86-64 (based on `manylinux_2_34`),
and for MacOS x86-64 and Aarch64 (for MacOS 13.3+). At the moment, the package is only tested
on these platforms.

Compiling from scratch uses CMake and the [Conan C/C++ package manager](https://docs.conan.io/2/)
to build the official [`ncbi-cxx-toolkit-public` recipe](https://github.com/ncbi/ncbi-cxx-toolkit-conan)
and install the artifacts into a self-contained Python package. If compilation fails for your
platform, please open an issue on the [issue tracker](https://github.com/althonos/pyncbitk/issues).


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The NCBI C++ Toolkit is a "United States Government Work" and therefore lies in
the public domain, but may be subject to copyright by the U.S. in foreign
countries. Some restrictions apply, see the
[NCBI C++ Toolkit license](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/doc/public/LICENSE).
