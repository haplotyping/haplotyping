import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import unittest, tempfile, logging, shutil, json
import haplotyping.data

class DatabaseTestCase(unittest.TestCase):
    
    pedigreeName = "pedigree.xlsx"
    postfix = ".xlsx"
    
    def setUp(self):
        
        logging.basicConfig(format="%(asctime)s | %(name)s |  %(levelname)s: %(message)s", datefmt="%m-%d-%y %H:%M:%S")
        logging.getLogger("haplotyping.index.database").setLevel(logging.ERROR)                
        
        self.resourceNames = []        
        self.dataLocation = os.path.abspath("./tests/data/testdata/") 
        
        try:            
            #create temporary directory
            self.tmpDirectory = tempfile.TemporaryDirectory()
            os.mkdir(os.path.join(self.tmpDirectory.name,"data"))
            os.mkdir(os.path.join(self.tmpDirectory.name,"service"))
            #copy excel-files to temporary directory 
            pedigreeFound = False
            for filename in os.listdir(self.dataLocation):
                if os.path.isfile(os.path.join(self.dataLocation,filename)):
                    if filename==self.pedigreeName:
                        shutil.copy(os.path.join(self.dataLocation,filename), os.path.join(self.tmpDirectory.name,"data"))
                        pedigreeFound = True
                    elif filename.endswith(self.postfix):                    
                        shutil.copy(os.path.join(self.dataLocation,filename), os.path.join(self.tmpDirectory.name,"data")) 
                        self.resourceNames.append(filename)

            #validate pedigree
            self.assertTrue(pedigreeFound, "no pedigree found") 
            validator = haplotyping.data.ValidatePedigree(os.path.join(self.tmpDirectory.name,"data"),self.pedigreeName)
            self.assertTrue(validator.valid, "pedigree did not validate") 
            #store report and package
            with open(os.path.join(self.tmpDirectory.name,"data",self._reportName(self.pedigreeName)), "w") as f:
                f.write(validator.createTextReport())
            validator.createPackageJSON(os.path.join(self.tmpDirectory.name,"data",
                                                     self._packageName(self.pedigreeName)))   

            #validate resource(s)
            self.assertTrue(len(self.resourceNames)>0, "no resource(s) found") 
            for filename in self.resourceNames:
                validator = haplotyping.data.ValidateData(os.path.join(self.tmpDirectory.name,"data"),
                                      filename,self._packageName(self.pedigreeName))
                self.assertTrue(validator.valid, "resource '{}' did not validate".format(filename)) 
                validator.createPackageJSON(os.path.join(self.tmpDirectory.name,"data",self._packageName(filename)))
                
            #create database
            haplotyping.data.ConstructDatabase(os.path.join(self.tmpDirectory.name,"data"),
                                               os.path.join(self.tmpDirectory.name,"service"))
        except:
            self.tmpDirectory.cleanup()
            
    def _packageName(self,filename):
        return os.path.splitext(filename)[0]+".package.json"
    
    def _reportName(self,filename):
        return os.path.splitext(filename)[0]+".report.txt"
                
    def test_validation(self):
        self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",self.pedigreeName)), "no pedigree file") 
        self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",
                                                    self._packageName(self.pedigreeName))), "no pedigree package") 
        self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",
                                                    self._reportName(self.pedigreeName))), "no pedigree report")   
        for filename in self.resourceNames:
            self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",filename)), 
                            "no copied resource file '{}'".format(filename)) 
            self.assertTrue(os.path.isfile(os.path.join(self.tmpDirectory.name,"data",
                                                    self._packageName(filename))), 
                            "no resource package '{}'".format(filename)) 
            
    def test_database(self):
        for filename in os.listdir(os.path.join(self.tmpDirectory.name,"service")):
            if os.path.isfile(os.path.join(self.tmpDirectory.name,"service",filename)):
                shutil.copy(os.path.join(self.tmpDirectory.name,"service",filename),"/Users/matthijs/Desktop/tmp")
        
                
    def tearDown(self):
        if self.tmpDirectory:
            self.tmpDirectory.cleanup()            
    
