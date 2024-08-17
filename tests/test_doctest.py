# coding: utf-8
"""Test doctest contained tests in every file of the module.
"""

import doctest
import importlib
import os
import sys
import types
import warnings
import pkgutil

import pyncbitk


def _load_tests_from_module(tests, module, globs, setUp=None, tearDown=None):
    """Load tests from module, iterating through submodules."""
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith("_")):
        if isinstance(attr, types.ModuleType):
            suite = doctest.DocTestSuite(
                attr,
                globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
            tests.addTests(suite)
    return tests


def load_tests(loader, tests, ignore):
    """`load_test` function used by unittest to find the doctests."""

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)

    # doctests are not compatible with `green`, so we may want to bail out
    # early if `green` is running the tests
    if sys.argv[0].endswith("green"):
        return tests

    # recursively traverse all library submodules and load tests from them
    packages = [None, pyncbitk]

    for pkg in iter(packages.pop, None):
        for (_, subpkgname, subispkg) in pkgutil.walk_packages(pkg.__path__):
            # do not import __main__ module to avoid side effects!
            if subpkgname.startswith(("__main__", "__test__", "tests")):
                continue
            # import the submodule and add it to the tests
            module = importlib.import_module(".".join([pkg.__name__, subpkgname]))
            globs = dict(
                **module.__dict__
            )
            # remove some duplicate tests declared by Cython
            if hasattr(module, "__test__") and hasattr(module, "__reduce_cython__"):
                module.__test__.clear()
            tests.addTests(
                doctest.DocTestSuite(
                    module,
                    globs=globs,
                    setUp=setUp,
                    tearDown=tearDown,
                    optionflags=+doctest.ELLIPSIS,
                )
            )
            # if the submodule is a package, we need to process its submodules
            # as well, so we add it to the package queue
            if subispkg and subpkgname != "tests":
                packages.append(module)
    return tests
