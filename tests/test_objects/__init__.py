# noqa: D104

from . import (
    test_general,
    test_seqdata,
    test_seqdesc,
    test_seqid,
    test_seqinst,
    test_seqloc,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_general))
    suite.addTests(loader.loadTestsFromModule(test_seqdata))
    suite.addTests(loader.loadTestsFromModule(test_seqdesc))
    suite.addTests(loader.loadTestsFromModule(test_seqid))
    suite.addTests(loader.loadTestsFromModule(test_seqinst))
    suite.addTests(loader.loadTestsFromModule(test_seqloc))
    return suite
