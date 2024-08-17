# noqa: D104

from . import (
    test_doctest,
)


def load_tests(loader, suite, pattern):
    test_doctest.load_tests(loader, suite, pattern)
    return suite
