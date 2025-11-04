# 🧬🧰 PyNCBItk [![Stars](https://img.shields.io/github/stars/althonos/pyncbitk.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyncbitk/stargazers)

*(Unofficial) [Cython](https://cython.org/) bindings and Python interface to the [NCBI C++ Toolkit](https://www.ncbi.nlm.nih.gov/toolkit).*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyncbitk/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyncbitk/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyncbitk?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyncbitk/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![PyPI](https://img.shields.io/pypi/v/pyncbitk.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyncbitk)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyncbitk?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyncbitk)
[![AUR](https://img.shields.io/aur/version/python-pyncbitk?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyncbitk)
[![Wheel](https://img.shields.io/pypi/wheel/pyncbitk.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyncbitk/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyncbitk.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyncbitk/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyncbitk.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyncbitk/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyncbitk/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square)](https://git.lumc.nl/mflarralde/pyncbitk/)
[![Issues](https://img.shields.io/github/issues/althonos/pyncbitk.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyncbitk/issues)
[![Docs](https://img.shields.io/readthedocs/pyncbitk/latest?style=flat-square&maxAge=600)](https://pyncbitk.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyncbitk/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyncbitk?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyncbitk)

***⚠️ This package is a work-in-progress and in a very experimental state. Expect segmentation faults, compilation issues, missing features, incomplete documentation.***

## 🗺️ Overview

The [NCBI C++ Toolkit](https://ncbi.github.io/cxx-toolkit/) is a framework of
C++ libraries to work with biological sequence data developed at the
[National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/).
It features a flexible object model for representing sequences of various
origin, including composite or virtual sequences; a resource manager
to easily manipulate heterogeneous data sources; and a comprehensive API to the
various BLAST algorithms[\[1\]](#ref1) developed at the NBCI.

PyNCBItk is a Python library that provides bindings to the NCBI C++ Toolkit
data model and BLAST+ interface using [Cython](https://cython.org). It exposes
the internals of the C++ Toolkit, allowing BLAST queries to be run directly
from the Python interpreter without external I/O.

## 📋 Roadmap

The package is in a very experimental state, and only a few core features are
supported at the moment:

- [x] Loading sequences from a FASTA file.
- [x] Creating basic sequences through the Python API.
- [x] Running BLAST searches with default parameters.
- [ ] Thorough BLAST configuration.
- [ ] Error and warning management.
- [ ] Support for all kinds of sequence storage.
- [ ] Multi-threading for database searches using Python threads.
- [ ] Advanced interface for the object manager.
- [ ] Interface for all sequence and alignment types.

## 🔧 Installing

PyNCBItk is available for all modern Python (3.7+). Compilation is done
through [CMake](https://cmake.org) using [Scikit-build-core](https://scikit-build-core.readthedocs.io).

To install an alpha release, use `pip` with the `--pre` flag:
```console
$ pip install pyncbitk --pre
```

The `pyncbitk` package requires additional runtime libraries that are distributed
in the `pyncbitk-runtime` package. These libraries should be available in a pre-compiled
wheel for Linux and MacOS platforms. Otherwise, they can be compiled on setup 
with the [Conan C/C++ package manager](https://docs.conan.io/2/)
to handle compilation of the NCBI C++ Toolkit. *The project will take ages to
compile the first time, but afterwards only the Cython code will have to be
recompiled.*

## 💡 Example

```python
from pyncbitk.objects.seqset import BioSeqSet
from pyncbitk.objtools import DatabaseReader, FastaReader
from pyncbitk.algo.blast import BlastN

# read the sequences from FASTA-formatted files
queries = BioSeqSet(FastaReader("queries.fna", split=False))
subjects = BioSeqSet(FastaReader("subjects.fna", split=False))

# run `blastn` with default parameters
blastn = BlastN()
results = blastn.run(queries, subjects)
```

The result is a `SearchResultsSet` which contains one `SearchResults` object
per query/subject pair. The `SearchResults` object summarizes the result
and contains the hit alignments in a `SeqAlignSet`.

See the [Examples section](https://pyncbitk.readthedocs.io/en/latest/examples/index.html) 
in the [online documentation](https://pyncbitk.readthedocs.io/en/latest/examples/index.html)
for more information.

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the
[GitHub issue tracker](https://github.com/althonos/pyncbitk/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyncbitk/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyncbitk/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The NCBI C++ Toolkit is a "United States Government Work" and therefore lies in
the public domain, but may be subject to copyright by the U.S. in foreign
countries. Some restrictions apply, see the
[NCBI C++ Toolkit license](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/doc/public/LICENSE).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the NCBI or any associated entity. It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD
at the [Leiden University Medical Center](https://www.lumc.nl/en/) in
the [Zeller team](https://github.com/zellerlab).*

## 📚 References

- <a id="ref1">\[1\]</a> Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. *Journal of molecular biology*, 215(3), 403–410. [doi:10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)
