import unittest
import pickle

from pyncbitk.objects.seqinst import *
from pyncbitk.objects.seqdata import *
from pyncbitk.objects.seqloc import WholeSeqLoc
from pyncbitk.objects.seqid import LocalId
from pyncbitk.objects.general import ObjectId


class TestSeqInst:
    pass


class TestContinuousInst(unittest.TestCase, TestSeqInst):

    def test_init(self):
        data = IupacNaData("ATGC")
        inst = ContinuousInst(data)
        self.assertEqual(inst.data, data)

    def test_init_molecule_detect(self):
        i1 = ContinuousInst(IupacNaData("ATGC"))
        self.assertEqual(i1.molecule, "dna")
        i2 = ContinuousInst(NcbiEAaData("ATGC"))
        self.assertEqual(i2.molecule, "protein")

    def test_pickle(self):
        i1 = ContinuousInst(IupacNaData("ATGC"))
        p1 = pickle.loads(pickle.dumps(i1))
        self.assertEqual(i1, p1)

    def test_repr(self):
        data = IupacNaData("ATGC")
        i1 = ContinuousInst(data)
        self.assertEqual(repr(i1), f"ContinuousInst({data!r}, strandedness='double', molecule='dna', length=4)")


class TestRefInst(unittest.TestCase, TestSeqInst):

    def test_init(self):
        seqloc = WholeSeqLoc(LocalId(ObjectId(1)))
        inst = RefInst(seqloc)
        self.assertEqual(inst.seqloc, seqloc)

    def test_pickle(self):
        seqloc = WholeSeqLoc(LocalId(ObjectId(1)))
        i1 = RefInst(seqloc)
        p1 = pickle.loads(pickle.dumps(i1))
        self.assertEqual(i1, p1)

    def test_repr(self):
        seqloc = WholeSeqLoc(LocalId(ObjectId(1)))
        i1 = RefInst(seqloc)
        self.assertEqual(repr(i1), f"RefInst({seqloc!r})")


class TestDeltaInst(unittest.TestCase, TestSeqInst):
    
    def test_init_empty(self):
        inst = DeltaInst()
        self.assertEqual(len(inst), 0)
        self.assertIs(inst.length, None)

    def test_init(self):
        d1 = LocDelta(WholeSeqLoc(LocalId(ObjectId(1))))
        d2 = LiteralDelta(4, IupacNaData("ATGC"))
        inst = DeltaInst([d1, d2])
        self.assertEqual(len(inst), 2)

    def test_pickle(self):
        d1 = LocDelta(WholeSeqLoc(LocalId(ObjectId(1))))
        d2 = LiteralDelta(4, IupacNaData("ATGC"))
        inst = DeltaInst([d1, d2])
        inst2 = pickle.loads(pickle.dumps(inst))
        self.assertEqual(inst, inst2)

    def test_repr(self):
        d1 = LocDelta(WholeSeqLoc(LocalId(ObjectId(1))))
        d2 = LiteralDelta(4, IupacNaData("ATGC"))
        inst = DeltaInst([d1, d2])
        self.assertEqual(repr(inst), f"DeltaInst({[d1, d2]!r})" )