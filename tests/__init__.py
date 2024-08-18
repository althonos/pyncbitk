# noqa: D104

from . import (
    test_doctest,
    test_objects
)


def load_tests(loader, suite, pattern):
    test_doctest.load_tests(loader, suite, pattern)
    suite.addTests(loader.loadTestsFromModule(test_objects))
    return suite
