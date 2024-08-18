import unittest

from pyncbitk.objects.seqdata import *


class TestSeqData:
    pass


class TestIupacNaData(TestSeqData, unittest.TestCase):

    def test_init_str(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.decode(), "ATGC")
        self.assertEqual(memoryview(data).tobytes(), b"ATGC")

    def test_init_bytes(self):
        data = IupacNaData(b"ATGC")
        self.assertEqual(data.decode(), "ATGC")
        self.assertEqual(memoryview(data).tobytes(), b"ATGC")

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

        data = IupacNaData("ATGC")
        self.assertEqual(data.length, 4)

    def test_decode(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.decode(), "ATGC")

    def test_data(self):
        data = IupacNaData("ATGC")
        self.assertEqual(data.data, "ATGC")

