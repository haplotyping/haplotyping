import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import unittest
from haplotyping.general import *

class GeneralTestCase(unittest.TestCase):
    
    def setUp(self):
        self.exampleKmer = "AAAACCCGGT"
        self.reverseComplementExampleKmer = "ACCGGGTTTT"
        self.canonicalExampleKmer = "AAAACCCGGT"
        self.exampleInvalidKmer = "AAAAqCCGGT"
        
        
    def test_reverse_complement(self):
        self.assertEqual(General.reverse_complement(self.exampleKmer),
                         self.reverseComplementExampleKmer,"incorrect reverse complement")
        self.assertEqual(General.reverse_complement(self.reverseComplementExampleKmer),
                         self.exampleKmer,"incorrect reverse complement")
        self.assertRaises(Exception,General.reverse_complement, self.exampleInvalidKmer)
        
    def test_canonical(self):
        self.assertEqual(General.canonical(self.exampleKmer),
                         self.canonicalExampleKmer,"incorrect canonical")
        self.assertEqual(General.canonical(self.reverseComplementExampleKmer),
                         self.canonicalExampleKmer,"incorrect canonical")
        self.assertRaises(Exception,General.canonical, self.exampleInvalidKmer)
        
    def tearDown(self):
        pass
    
