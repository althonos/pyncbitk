Installation
============

.. note::

    Wheels are provided for Linux x86-64 platforms. Other machines will have
    to build the wheel from the source distribution. Building ``pyncbitk``
    involves compiling the NCBI C++ Toolkit, which requires a C++ compiler
    to be available.


.. PyPi
.. ^^^^

.. ``pyncbitk`` is hosted on GitHub, but the easiest way to install it is to download
.. the latest release from its `PyPi repository <https://pypi.python.org/pypi/pyncbitk>`_.
.. It will install all dependencies then install ``pyncbitk`` either from a wheel if
.. one is available, or from source after compiling the Rust code :

.. .. code:: console

.. 	$ pip install --user pyncbitk


.. Conda
.. ^^^^^

.. PyNCBItk is also available as a `recipe <https://anaconda.org/bioconda/pyncbitk>`_
.. in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
.. use the ``conda`` installer:

.. .. code:: console

..    $ conda install -c bioconda pyncbitk


.. Arch User Repository
.. ^^^^^^^^^^^^^^^^^^^^

.. A package recipe for Arch Linux can be found in the Arch User Repository
.. under the name `python-pyncbitk <https://aur.archlinux.org/packages/python-pyncbitk>`_.
.. It will always match the latest release from PyPI.

.. Steps to install on ArchLinux depend on your `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
.. (``yaourt``, ``aura``, ``yay``, etc.). For ``aura``, you'll need to run:

.. .. code:: console

..     $ aura -A python-pyncbitk


.. BioArchLinux
.. ^^^^^^^^^^^^

.. The `BioArchLinux <https://bioarchlinux.org>`_ project provides pre-compiled packages
.. based on the AUR recipe. Add the BioArchLinux package repository to ``/etc/pacman.conf``:

.. .. code:: ini

..     \[bioarchlinux\]
..     Server = https://repo.bioarchlinux.org/$arch

.. Then install the latest version of the package and its dependencies with ``pacman``:

.. .. code:: console

..     $ pacman -S python-pyncbitk


.. Piwheels
.. ^^^^^^^^

.. ``pyncbitk`` works on Raspberry Pi computers, and pre-built wheels are compiled
.. for `armv7l` platforms on piwheels. Run the following command to install these
.. instead of compiling from source:

.. .. code:: console

..    $ pip3 install pyncbitk --extra-index-url https://www.piwheels.org/simple

.. Check the `piwheels documentation <https://www.piwheels.org/faq.html>`_ for
.. more information.


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

   $ git clone --recursive https://github.com/althonos/pyncbitk
   $ pip install --user ./pyncbitk

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``build``
^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
compile manually, although you will need to install the ``build`` package
first:

.. code:: console

   $ git clone --recursive https://github.com/althonos/pyncbitk
   $ cd pyncbitk
   $ python -m build . --wheel
   # python -m installer pyncbitk*.whl

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.
