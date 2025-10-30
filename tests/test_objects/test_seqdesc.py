import unittest
import pickle

from pyncbitk.objects.seqdesc import *


class TestSeqDesc:
    pass


class TestNameDesc(unittest.TestCase):

    def test_init(self):
        name = "something"
        desc = NameDesc(name)
        self.assertEqual(str(desc), name)

    def test_init_error(self):
        self.assertRaises(TypeError, NameDesc, 1)
        self.assertRaises(TypeError, NameDesc, [])
        self.assertRaises(TypeError, NameDesc, object())

    def test_pickle(self):
        n1 = NameDesc("a")
        n2 = pickle.loads(pickle.dumps(n1))
        self.assertEqual(n1, n2)

    def test_repr(self):
        name = "something"
        desc = NameDesc(name)
        self.assertEqual(repr(desc), f"NameDesc({name!r})")
    
    def test_eq(self):
        n1 = NameDesc("desc 1")
        n2 = NameDesc("desc 1")
        n3 = NameDesc("desc 3")
        self.assertEqual(n1, n1)
        self.assertEqual(n1, n2)
        self.assertNotEqual(n1, n3)
        self.assertNotEqual(n2, n3)

    def test_hash(self):
        n1 = NameDesc("desc 1")
        n2 = NameDesc("desc 1")
        n3 = NameDesc("desc 3")
        self.assertEqual(hash(n1), hash(n1))
        self.assertEqual(hash(n1), hash(n2))
        self.assertNotEqual(hash(n1), hash(n3))
        self.assertNotEqual(hash(n2), hash(n3))