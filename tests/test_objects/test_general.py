import unittest
import pickle

from pyncbitk.objects.general import *

class TestObjectId(unittest.TestCase):

    def test_init_int(self):
        oid = ObjectId(1)
        self.assertEqual(oid.value, 1)
    
    def test_init_str(self):
        oid = ObjectId("sequence 1")
        self.assertEqual(oid.value, "sequence 1")

    def test_init_bytes(self):
        oid = ObjectId(b"sequence 1")
        self.assertEqual(oid.value, "sequence 1")

    def test_init_err(self):
        self.assertRaises(TypeError, ObjectId, [])

    def test_pickle(self):
        oid1 = ObjectId(1)
        cpy1 = pickle.loads(pickle.dumps(oid1))
        self.assertEqual(oid1, cpy1)

        oid2 = ObjectId("sequence 1")
        cpy2 = pickle.loads(pickle.dumps(oid2))
        self.assertEqual(oid2, cpy2)

    def test_repr(self):
        oid1 = ObjectId(1)
        self.assertEqual(repr(oid1), "ObjectId(1)")
        oid2 = ObjectId("sequence 1")
        self.assertEqual(repr(oid2), "ObjectId('sequence 1')")