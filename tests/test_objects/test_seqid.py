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
    
    def test_eq(self):
        id1 = LocalId(ObjectId(1))
        id2 = LocalId(ObjectId(1))
        id3 = LocalId(ObjectId(2))
        self.assertEqual(id1, id1)
        self.assertEqual(id1, id2)
        self.assertNotEqual(id1, id3)

    def test_hash(self):
        id1 = LocalId(ObjectId(1))
        id2 = LocalId(ObjectId(1))
        id3 = LocalId(ObjectId(2))
        self.assertEqual(hash(id1), hash(id1))
        self.assertEqual(hash(id1), hash(id2))
        self.assertNotEqual(hash(id1), hash(id3))