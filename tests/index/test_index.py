import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import unittest, tempfile, logging, h5py, gzip, csv, shutil, pytest
import numpy as np
from haplotyping.index.database import *

class IndexTestCase(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        
        logging.basicConfig(format="%(asctime)s | %(name)s |  %(levelname)s: %(message)s", datefmt="%m-%d-%y %H:%M:%S")
        logging.getLogger("haplotyping.index.database").setLevel(logging.ERROR)
        
        try:
            self.serviceDataLocation = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../data/testdata/")
            self.dataLocation = os.path.join(os.path.abspath(os.path.dirname(__file__)), "testdata/")
            self.sortedListLocation = os.path.join(self.dataLocation,"kmer.list.sorted.gz")
            #get k from sorted list
            with gzip.open(self.sortedListLocation, "rt") as f: 
                reader = csv.reader(f, delimiter="\t")
                for line in reader:
                    self.k = len(line[0])
                    break
            self.name = "Test hdf5 database"
            self.minimumFrequency = 2
            #find reads
            (self.unpairedReadFiles, self.pairedReadFiles, 
                 self.allReadFiles) = haplotyping.index.Database.detectReadFiles(self.dataLocation)
            #compute read statistics directly
            self.readUnpairedTotal=0
            self.readPairedTotal=0
            self.minReadLength = float("inf")
            self.maxReadLength = 0
            for filename in self.unpairedReadFiles:
                with gzip.open(filename, "rb") as f:
                    for i, l in enumerate(f):
                        if i%4==1:
                            self.maxReadLength=max(self.maxReadLength,len(l.decode().strip()))
                            self.minReadLength=min(self.minReadLength,len(l.decode().strip()))
                    self.readUnpairedTotal+=round((i+1)/4)
            for pairedReadFile in self.pairedReadFiles:
                for filename in pairedReadFile:
                    with gzip.open(filename, "rb") as f:
                        for i, l in enumerate(f):
                            if i%4==1:
                                self.maxReadLength=max(self.maxReadLength,len(l.decode().strip()))
                                self.minReadLength=min(self.minReadLength,len(l.decode().strip()))
                        self.readPairedTotal+=round((i+1)/4)
            #create temporary index
            self.tmpDirectory = tempfile.TemporaryDirectory()            
            self.tmpIndexLocation = self.tmpDirectory.name+"/kmer.data.h5"
            #create database
            self.database = haplotyping.index.Database(self.k, self.name, self.tmpDirectory.name+"/kmer.data", 
                                          self.sortedListLocation , self.unpairedReadFiles, self.pairedReadFiles,
                                          minimumFrequency=self.minimumFrequency)
            #use as service testdata
            if os.path.isfile(self.tmpIndexLocation):
                os.makedirs(os.path.join(self.serviceDataLocation,"kmer/demoset/test/agria"), exist_ok=True)
                shutil.copy2(os.path.join(self.dataLocation,"kmer.kmc.kmc_pre"),
                             os.path.join(self.serviceDataLocation,"kmer/demoset/test/agria"))
                shutil.copy2(os.path.join(self.dataLocation,"kmer.kmc.kmc_suf"),
                             os.path.join(self.serviceDataLocation,"kmer/demoset/test/agria"))
                shutil.copy2(self.tmpIndexLocation,
                             os.path.join(self.serviceDataLocation,"kmer/demoset/test/agria/kmer.data.h5"))
        except:
            self.tmpDirectory.cleanup()
                
    def test_configuration(self):
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            self.assertTrue("/config" in h5file,"couldn't find configuration group")
            self.assertEqual(h5file["/config"].attrs["k"],self.k,"unexpected k-mer size")
            self.assertEqual(h5file["/config"].attrs["name"],self.name,"unexpected index name")
            self.assertEqual(h5file["/config"].attrs["minimumKmerFrequencies"],
                             self.minimumFrequency,"unexpected minimum frequency")
            self.assertEqual(h5file["/config"].attrs["minimumReadLength"],self.minReadLength,
                             "unexpected minimum read length")
            self.assertEqual(h5file["/config"].attrs["maximumReadLength"],self.maxReadLength,
                             "unexpected maximum read length")
            self.assertEqual(h5file["/config"].attrs["numberReadsUnpaired"],self.readUnpairedTotal,
                             "unexpected unpaired number of reads")
            self.assertEqual(h5file["/config"].attrs["numberReadsPaired"],self.readPairedTotal,
                             "unexpected paired number of reads")
            self.assertEqual(h5file["/config"].attrs["numberReads"],self.readPairedTotal+self.readUnpairedTotal,
                             "unexpected number of reads")
            
    def test_ckmer(self):
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            ckmers = h5file.get("/split/ckmer")
            bases = h5file.get("/split/base")
            ckmerTotal = ckmers.shape[0]
            #test for empty
            self.assertTrue(ckmerTotal>0,"no splitting k-mers")
            previousCkmer = None
            for row in ckmers:
                if not previousCkmer == None:
                    self.assertTrue(row[0]>previousCkmer,"k-mers not properly sorted")
                previousCkmer = row[0]
                if row[1].decode()=="l" or row[1].decode()=="b":
                    baseLeft = bases[row[3][0]]
                    self.assertEqual(haplotyping.General.reverse_complement(row[0].decode())[:-1],baseLeft[0].decode(),
                                     "wrong base left split")
                if row[1].decode()=="r" or row[1].decode()=="b":
                    baseRight = bases[row[3][1]]
                    self.assertEqual(row[0].decode()[:-1],baseRight[0].decode(),
                                     "wrong base right split")
                self.assertTrue(row[2]>=row[4][1][0],"inconsistent direct left distinct")
                self.assertTrue(row[2]>=row[4][1][1],"inconsistent direct left number")
                self.assertTrue(row[2]>=row[4][2][0],"inconsistent direct right distinct")
                self.assertTrue(row[2]>=row[4][2][1],"inconsistent direct right number")
                self.assertTrue(row[2]>=row[6][0],"inconsistent cycle")
                self.assertTrue(row[2]>=row[7][0],"inconsistent reversal")
            #test cycle
            for row in h5file["/relations/cycle"]:
                self.assertTrue(row[0]<ckmerTotal,"linking to non-existing entry")
                self.assertTrue(row[1]>0,"length should be positive")
                self.assertTrue(row[2]>0,"number should be positive")
                self.assertTrue(ckmers[row[0]][6][0]==row[2],"k-mer labelled with incorrect number as cycle")
            #test reversal
            for row in h5file["/relations/reversal"]:
                self.assertTrue(row[0]<ckmerTotal,"linking to non-existing entry")
                self.assertTrue(row[1]>0,"length should be positive")
                self.assertTrue(row[2]>0,"number should be positive")
                self.assertTrue(ckmers[row[0]][7][0]==row[2],"k-mer labelled with incorrect number as reversal")
                
    def test_base(self):
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            ckmers = h5file.get("/split/ckmer")
            bases = h5file.get("/split/base")
            baseTotal = bases.shape[0]
            #test for empty
            self.assertTrue(baseTotal>0,"no bases")
            previousBase = None
            for row in bases:
                if not previousBase == None:
                    self.assertTrue(row[0]>previousBase,"bases not properly sorted")
                previousBase = row[0]
                self.assertTrue(row[1]>0,"base without k-mers")
                letterTotal = 0
                for i in range(len(haplotyping.index.Database.letters)):
                    letter = haplotyping.index.Database.letters[i]
                    if row[2][i][0]>0:
                        ckmerRow = ckmers[row[2][i][1]]
                        letterTotal += row[2][i][0]
                        ckmer = haplotyping.General.canonical(row[0].decode()+letter)
                        self.assertEqual(ckmer,ckmerRow[0].decode(),"branch should equal linked k-mer")
                        self.assertEqual(row[2][i][0],ckmerRow[2],"branch number should equal number of linked k-mer")
                self.assertEqual(row[1],letterTotal,"base number should equal sum of branch numbers")                    
            
    def test_direct(self):
        directList = []
        with h5py.File(self.tmpIndexLocation,"r") as h5file:
            ckmerTotal = h5file["/split/ckmer"].shape[0]
            #test for empty
            self.assertTrue(h5file["/relations/direct"].shape[0]>0,"no direct relations")
            #create list with keys, skip problematic
            for row in h5file["/relations/direct"]:
                directList.append((row[0][0],row[0][1],row[1][0],row[1][1],row[2],row[3],))
            directSet = set(directList)
            #test for unique information
            self.assertEqual(len(directSet),len(directList),"duplication in direct relations")
            #test for symmetry, skip problematic
            for row in h5file["/relations/direct"]:
                self.assertTrue((row[1][0],row[1][1],row[0][0],row[0][1],row[2],row[3],) 
                                    in directList,"direct relations not symmetric")              
                    
    @classmethod
    def tearDownClass(self):
        if self.tmpDirectory:            
            self.tmpDirectory.cleanup()
    
