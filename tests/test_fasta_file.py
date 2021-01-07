import unittest
from libs.fasta_file import FastaFile


class TestFastaFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_files_dir = "test_files"

    def test_read_fast_file(self):
        expected = [
            {"sequence": "aaaatttt", "title": ">random sequence 1"},
            {"sequence": "ttttgggg", "title": ">random sequence 2"},
            {"sequence": "ggggaaaa", "title": ">random sequence 3"}
        ]
        filepath = "{}/read.fasta".format(self.test_files_dir)
        result = FastaFile.read(filepath)
        self.assertEqual(len(expected), len(result))

        for i in range(len(expected)):
            self.assertEqual(expected[i], result[i])


if __name__ == '__main__':
    unittest.main()
