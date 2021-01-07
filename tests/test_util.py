import unittest
from libs import util


class TestUtil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_files_dir = "test_files"

    def test_one(self):
        source = "ggggatcagggg"
        target = "ggggattagggg"
        res = util.get_fuzzy_token_ratio(source, target)
        print(res)




if __name__ == '__main__':
    unittest.main()
