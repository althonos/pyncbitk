# noqa: D104

from . import (
    test_seqdata,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_seqdata))
    return suite
