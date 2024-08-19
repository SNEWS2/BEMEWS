from importlib import import_module
import unittest
import sys

class TestInit(unittest.TestCase):
    def test_import(self):
        import_module('BEMEWS')
        self.assertTrue('BEMEWS' in sys.modules)