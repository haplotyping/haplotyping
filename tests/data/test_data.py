import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import unittest, tempfile, logging, shutil, json, pytest
import haplotyping.data

class ServiceDataTestCase(unittest.TestCase):
    
    splitDatabase = "kmer.data.h5"
    kmerDatabase = "kmer.kmc"
    indexName = "index.xlsx"
    postfix = ".xlsx"
        
    @classmethod
    def setUpClass(self):
        
        logging.basicConfig(format="%(asctime)s | %(name)s |  %(levelname)s: %(message)s", datefmt="%m-%d-%y %H:%M:%S")
        logging.getLogger("haplotyping.index.database").setLevel(logging.ERROR)   
        
        try:            
            self.resourceNames = []
            self.serviceLocation = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../service/testdata/")
            self.dataLocation = os.path.join(os.path.abspath(os.path.dirname(__file__)), "testdata/")
            #create temporary directory
            self.tmpDirectory = tempfile.TemporaryDirectory()
            os.mkdir(os.path.join(self.tmpDirectory.name,"data"))
            os.mkdir(os.path.join(self.tmpDirectory.name,"service"))
            #temporary locations
            self.tmpDatabaseLocation = os.path.join(self.tmpDirectory.name,"service","db.sqlite")
            self.tmpIdentifiersLocation = os.path.join(self.tmpDirectory.name,"service","identifiers.json")
            self.tmpKmerLocation = os.path.join(self.tmpDirectory.name,"service","kmer")
            self.tmpMarkerLocation = os.path.join(self.tmpDirectory.name,"service","marker")
            #copy excel-files and k-mer databases to temporary directory 
            self.indexFound = False
            for filename in os.listdir(self.dataLocation):
                if os.path.isfile(os.path.join(self.dataLocation,filename)):
                    if filename==self.indexName:
                        shutil.copy(os.path.join(self.dataLocation,filename), os.path.join(self.tmpDirectory.name,"data"))
                        self.indexFound = True
                    elif filename.endswith(self.postfix):                    
                        shutil.copy(os.path.join(self.dataLocation,filename), os.path.join(self.tmpDirectory.name,"data")) 
                        self.resourceNames.append(filename) 
            os.makedirs(self.tmpKmerLocation)
            os.makedirs(self.tmpMarkerLocation)
            shutil.copytree(os.path.join(self.dataLocation,"kmer/"),
                            self.tmpKmerLocation, 
                            dirs_exist_ok=True)
        except Exception as e:
            self.tmpDirectory.cleanup()
            
    def _packageName(self,filename):
        return os.path.splitext(filename)[0]+".package.json"
    
    def _reportName(self,filename):
        return os.path.splitext(filename)[0]+".report.txt"

    @pytest.mark.run(order=1)
    def test_validate_index(self):
        self.assertTrue(self.indexFound, "no index found") 
        validator = haplotyping.data.ValidateIndex(os.path.join(self.tmpDirectory.name,"data"),self.indexName)
        self.assertTrue(validator.valid, "index did not validate") 
        #store report and package
        with open(os.path.join(self.tmpDirectory.name,"data",self._reportName(self.indexName)), "w") as f:
            f.write(validator.createTextReport())
        validator.createPackageJSON(os.path.join(self.tmpDirectory.name,"data",
                                                 self._packageName(self.indexName)))
        
        
    @pytest.mark.run(order=2)
    def test_validate_resources(self):
        self.assertTrue(len(self.resourceNames)>0, "no resource(s) found") 
        for filename in self.resourceNames:
            validator = haplotyping.data.ValidateData(os.path.join(self.tmpDirectory.name,"data"),
                                  filename,self._packageName(self.indexName))
            self.assertTrue(validator.valid, "resource '{}' did not validate".format(filename)) 
            validator.createPackageJSON(os.path.join(self.tmpDirectory.name,"data",self._packageName(filename)))
            
    @pytest.mark.run(order=3)
    def test_validation(self):
        self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",self.indexName)), "no index file") 
        self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",
                                                    self._packageName(self.indexName))), "no index package") 
        self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",
                                                    self._reportName(self.indexName))), "no index report")   
        for filename in self.resourceNames:
            self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",filename)), 
                            "no copied resource file '{}'".format(filename)) 
            self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",
                                                    self._packageName(filename))), 
                            "no resource package '{}'".format(filename))        
            
    @pytest.mark.run(order=4)
    def test_construct_database(self):
        #create database
        if os.path.isfile(self.tmpDatabaseLocation):
            os.remove(self.tmpDatabaseLocation)
        haplotyping.data.ConstructDatabase(os.path.join(self.tmpDirectory.name,"data"),
                                               os.path.join(self.tmpDirectory.name,"service"))
        self.assertTrue(os.path.isfile(self.tmpDatabaseLocation), "no database")   
        self.assertTrue(os.path.isfile(self.tmpIdentifiersLocation), "no identifiers")
        
    @pytest.mark.run(order=5)
    def test_check_database(self):  
        #check database
        haplotyping.data.CheckDatabase(os.path.join(self.tmpDirectory.name,"service"))
        #use as service testdata
        if os.path.isfile(self.tmpDatabaseLocation):
            shutil.rmtree(os.path.join(self.serviceLocation,"kmer"), ignore_errors=True)
            shutil.rmtree(os.path.join(self.serviceLocation,"marker"), ignore_errors=True)
            os.makedirs(os.path.join(self.serviceLocation,"kmer"))
            os.makedirs(os.path.join(self.serviceLocation,"marker"))
            shutil.copytree(self.tmpKmerLocation,
                            os.path.join(self.serviceLocation,"kmer/"), 
                            dirs_exist_ok=True)
            shutil.copytree(self.tmpMarkerLocation,
                            os.path.join(self.serviceLocation,"marker/"), 
                            dirs_exist_ok=True)
            #copy database
            shutil.copy2(self.tmpDatabaseLocation,os.path.join(self.serviceLocation,"db.sqlite"))            
                
    @classmethod
    def tearDownClass(self):
        if self.tmpDirectory:
            self.tmpDirectory.cleanup()            
    
