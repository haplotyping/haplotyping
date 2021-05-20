import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import unittest, tempfile, logging, h5py
from haplotyping.index.database import *

class DatabaseTestCase(unittest.TestCase):
    
    def setUp(self):
        
        logging.basicConfig(format="%(asctime)s | %(name)s |  %(levelname)s: %(message)s", datefmt="%m-%d-%y %H:%M:%S")
        logging.getLogger("haplotyping.index.database").setLevel(logging.ERROR)
        
        try:
            #set variables
            self.k = 31
            self.name = "Test hdf5 database"
            self.minimumFrequency = 2
            #create temporary index
            self.tmpDirectory = tempfile.TemporaryDirectory()
            self.dataLocation = os.path.abspath("./tests/index/testdata/") 
            self.tmpIndexLocation = self.tmpDirectory.name+"/kmer.data.h5"
            self.reads = [self.dataLocation+"/reads.fastq.gz"]
            self.pairedReads = [[self.dataLocation+"/reads_R1.fastq.gz", self.dataLocation+"/reads_R2.fastq.gz"]]
            self.database = haplotyping.index.Database(31, self.name, self.tmpDirectory.name+"/kmer.data", 
                                          self.dataLocation+"/kmer.list.sorted.gz", self.reads, self.pairedReads,
                                          minimumFrequency=self.minimumFrequency)
        except:
            self.tmpDirectory.cleanup()
                
    def test_configuration(self):
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            self.assertTrue("/config" in h5file,"couldn't find configuration group")
            self.assertEqual(h5file["/config"].attrs["k"],self.k,"unexpected k-mer size")
            self.assertEqual(h5file["/config"].attrs["name"],self.name,"unexpected index name")
            self.assertEqual(h5file["/config"].attrs["minimumFrequency"],
                             self.minimumFrequency,"unexpected minimum frequency")
            self.assertEqual(h5file["/config"].attrs["readLengthMinimum"],151,"unexpected minimum read length")
            self.assertEqual(h5file["/config"].attrs["readLengthMaximum"],151,"unexpected maximum read length")
            self.assertEqual(h5file["/config"].attrs["readUnpairedTotal"],2500,"unexpected unpaired number of reads")
            self.assertEqual(h5file["/config"].attrs["readPairedTotal"],5000,"unexpected paired number of reads")
            self.assertEqual(h5file["/config"].attrs["readTotal"],7500,"unexpected number of reads")
            
    def test_ckmer(self):
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            #test for empty
            self.assertTrue(h5file["/split/ckmer"].shape[0]>0,"no splitting k-mers")
            previousCkmer = None
            for row in np.array(h5file["/split/ckmer"]):
                if not previousCkmer == None:
                    self.assertTrue(row[0]>previousCkmer,"k-mers not properly sorted")
                previousCkmer = row[0]
                self.assertTrue(row[2]>=row[3][0][0],"inconsistent direct left distinct")
                self.assertTrue(row[2]>=row[3][0][1],"inconsistent direct left number")
                self.assertTrue(row[2]>=row[3][1][0],"inconsistent direct right distinct")
                self.assertTrue(row[2]>=row[3][1][1],"inconsistent direct right number")
                self.assertTrue(row[2]>=row[4][0],"inconsistent connected")
                self.assertTrue(row[2]>=row[5][0],"inconsistent cycle")
                self.assertTrue(row[2]>=row[6][0],"inconsistent reversal")
            
    def test_direct(self):
        directList = []
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            #test for empty
            self.assertTrue(h5file["/relations/direct"].shape[0]>0,"no direct relations")
            #create list with keys
            for row in np.array(h5file["/relations/direct"]):
                directList.append((row[0],row[1],row[2],row[3],row[4],row[5],))
            directSet = set(directList)
            #test for unique information
            self.assertEqual(len(directSet),len(directList),"duplication in direct relations")
            #test for symmetry
            for row in np.array(h5file["/relations/direct"]):
                self.assertTrue((row[2],row[3],row[0],row[1],row[4],row[5],) 
                                in directList,"direct relations not symmetric")            
        
    def tearDown(self):
        if self.tmpDirectory:
            self.tmpDirectory.cleanup()
    
