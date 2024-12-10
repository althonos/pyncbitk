import unittest
import pickle

from pyncbitk.objects.general import ObjectId
from pyncbitk.objects.seqid import *


class TestSeqId:
    pass


class TestLocalId(unittest.TestCase):

    def test_init(self):
        obj_id = ObjectId(1)
        seq_id = LocalId(obj_id)
        self.assertEqual(seq_id.object_id, obj_id)

    def test_pickle(self):
        seq_id = LocalId(ObjectId(1))
        seq_id2 = pickle.loads(pickle.dumps(seq_id))
        self.assertEqual(seq_id, seq_id2)

    def test_repr(self):
        obj_id = ObjectId(1)
        seq_id = LocalId(obj_id)
        self.assertEqual(repr(seq_id), f"LocalId({obj_id!r})")
        