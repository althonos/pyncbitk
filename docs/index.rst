PyNCBItk |Stars|
================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyncbitk.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyncbitk/stargazers
   :class: dark-light

*Cython bindings and Python interface to the NCBI C++ Toolkit.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyncbitk/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyncbitk/actions
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyncbitk?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyncbitk/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/pyncbitk.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyncbitk
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyncbitk?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyncbitk
   :class: dark-light

.. |AUR| image:: https://img.shields.io/aur/version/python-pyncbitk?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyncbitk
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyncbitk?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyncbitk/#files
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyncbitk.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyncbitk/#files
   :class: dark-light

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyncbitk.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyncbitk/#files
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyncbitk/
   :class: dark-light

.. |Mirror| image:: https://img.shields.io/badge/mirror-LUMC-003EAA.svg?style=flat-square&maxAge=3600
   :target:https://git.lumc.nl/mflarralde/pyncbitk/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyncbitk.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyncbitk/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/pyncbitk?style=flat-square&maxAge=3600
   :target: http://pyncbitk.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyncbitk/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyncbitk?style=flat-square&color=303f9f&maxAge=3600&label=downloads
   :target: https://pepy.tech/project/pyncbitk
   :class: dark-light


Overview
--------

The NCBI C++ Toolkit is a framework of C++ libraries to work with biological
sequence data developed at the 
`National Center for Biotechnology Information <https://www.ncbi.nlm.nih.gov/>`_. 
It features a flexible object model for representing sequences of various 
origin, including composite or virtual sequences; a resource manager
to easily manipulate heterogeneous data sources; and a comprehensive API to the
various BLAST algorithms developed at the NBCI.

PyNCBItk is a Python library that provides bindings to the NCBI C++ Toolkit 
data model and BLAST+ interface using `Cython <https://cython.org>`_:

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyncbitk`` as a ``pip`` or ``conda`` dependency, no need
      for the BLAST+ binaries or any external dependency.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible

      Load a `BioSeq` from a :wiki:`FASTA format` file or create 
      it programmatically through the :doc:`Python API <api/index>`.

   .. grid-item-card:: :fas:`gears` Practical

      Retrieve results as they become available as dedicated 
      `SearchResults` objects.


Setup
-----

PyNCBItk is available for all modern Python versions (3.7+).

Run ``pip install pyncbitk`` in a shell to download the latest release from PyPI,
or have a look at the :doc:`Installation page <guide/install>` to find other ways 
to install PyNCBItk.

Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   Examples <examples/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst

License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
The NCBI C++ Toolkit is a "United States Government Work" and therefore lies in 
the public domain. Some restrictions apply, see the upstream 
`license files <https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/doc/public/LICENSE>`_.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the NCBI or any associated entity. It was developed
by* `Martin Larralde <https://github.com/althonos/>`_ *during his PhD project
at the* `Leiden University Medical Center <https://www.lumc.nl/en/>`_ *in
the* `Zeller team <https://github.com/zellerlab>`_.
