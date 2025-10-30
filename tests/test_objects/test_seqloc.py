import unittest
import pickle

from pyncbitk.objects.seqinst import *
from pyncbitk.objects.seqdata import *
from pyncbitk.objects.seqloc import *
from pyncbitk.objects.seqid import LocalId
from pyncbitk.objects.general import ObjectId


class TestSeqLoc:
    pass


class TestWholeSeqLoc(unittest.TestCase, TestSeqLoc):

    def test_init(self):
        id_ = LocalId(ObjectId(1))
        loc = WholeSeqLoc(id_)
        self.assertEqual(loc.sequence_id, id_)

    def test_pickle(self):
        loc = WholeSeqLoc(LocalId(ObjectId(1)))
        cpy = pickle.loads(pickle.dumps(loc))
        self.assertEqual(loc, cpy)

    def test_repr(self):
        id_ = LocalId(ObjectId(1))
        loc = WholeSeqLoc(id_)
        self.assertEqual(repr(loc), f"WholeSeqLoc({id_!r})")
        self.assertEqual(loc.sequence_id, id_)



class TestSeqIntervalLoc(unittest.TestCase, TestSeqLoc):

    def test_init(self):
        id_ = LocalId(ObjectId(1))
        loc = SeqIntervalLoc(id_, 0, 1)
        self.assertEqual(loc.sequence_id, id_)
        self.assertEqual(loc.start, 0)
        self.assertEqual(loc.stop, 1)

    def test_init_value_error(self):
        id_ = LocalId(ObjectId(1))
        with self.assertRaises(ValueError):
            loc = SeqIntervalLoc(id_, 10, 2)

    def test_pickle(self):
        id_ = LocalId(ObjectId(1))
        loc = SeqIntervalLoc(id_, 0, 1)
        cpy = pickle.loads(pickle.dumps(loc))
        self.assertEqual(loc, cpy)

    def test_repr(self):
        id_ = LocalId(ObjectId(1))
        loc = SeqIntervalLoc(id_, 0, 1)
        self.assertEqual(repr(loc), f"SeqIntervalLoc({id_!r}, start=0, stop=1)")
        self.assertEqual(loc.sequence_id, id_)

    