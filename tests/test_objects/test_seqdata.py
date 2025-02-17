import unittest
import pickle

from pyncbitk.objects.seqdata import *


class TestSeqData:
    datatype = None
    
    def test_encode_decode(self):
        data = self.datatype.encode("ATGC")
        self.assertEqual(data.decode(), "ATGC")

    def test_length(self):
        data = self.datatype.encode("ATGC")
        self.assertEqual(data.length, 4)
        data = self.datatype.encode("ATG")
        self.assertEqual(data.length, 3)


class TestIupacNaData(TestSeqData, unittest.TestCase):
    datatype = IupacNaData 

    def test_init_str(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.decode(), "ATGC")

    def test_init_bytes(self):
        data = IupacNaData(b"ATGC")
        self.assertEqual(data.decode(), "ATGC")

    def test_memoryview(self):
        data = IupacNaData(b"ATGC")
        mem = memoryview(data)
        self.assertEqual(mem.tobytes(), b"ATGC")
        self.assertEqual(mem.shape[0], 4)

    def test_init_error(self):
        with self.assertRaises(TypeError):
            data = IupacNaData(1)
        with self.assertRaises(TypeError):
            data = IupacNaData(None)
        with self.assertRaises(ValueError):
            data = IupacNaData("12345")
    
    def test_length(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.length, 4)
        data = IupacNaData("ATG")
        self.assertEqual(data.length, 3)

    def test_data(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.data, b"ATGC")

    def test_decode(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.decode(), "ATGC")

    def test_cmp(self):
        d1 = IupacNaData("ATGC")
        d2 = IupacNaData("ATGC")
        self.assertEqual(d1, d2)
        self.assertLessEqual(d1, d2)
        self.assertGreaterEqual(d1, d2)

        d3 = IupacNaData("ATG")
        self.assertLess(d3, d1)
        self.assertLessEqual(d3, d1)
        self.assertGreater(d1, d3)
        self.assertGreaterEqual(d1, d3)
        self.assertNotEqual(d3, d1)

    def test_pickle(self):
        d1 = IupacNaData("ATGC")
        d2 = pickle.loads(pickle.dumps(d1))
        self.assertEqual(d1, d2)
        d3 = pickle.loads(pickle.dumps(d1, protocol=pickle.HIGHEST_PROTOCOL))
        self.assertEqual(d1, d3)


# class TestNcbi2NaData(TestSeqData, unittest.TestCase):
#     datatype = Ncbi2NaData


# class TestNcbi4NaData(TestSeqData, unittest.TestCase):
#     datatype = Ncbi4NaData


# class TestNcbi8NaData(TestSeqData, unittest.TestCase):
#     datatype = Ncbi8NaData


# class TestNcbiPNaData(TestSeqData, unittest.TestCase):
#     datatype = NcbiPNaData


# class TestIupacAaData(TestSeqData, unittest.TestCase):
#     datatype = IupacAaData

# class TestNcbi8AaData(TestSeqData, unittest.TestCase):
#     datatype = Ncbi8AaData


class TestNcbiEAaData(TestSeqData, unittest.TestCase):
    datatype = NcbiEAaData


# class TestNcbiPAaData(TestSeqData, unittest.TestCase):
#     datatype = NcbiPAaData


# class TestNcbiStdAa(TestSeqData, unittest.TestCase):
#     datatype = NcbiStdAa

