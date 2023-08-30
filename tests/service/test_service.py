import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import unittest, tempfile, json, pathlib, urllib, socket, configparser

import haplotyping, haplotyping.service

class ServiceTestCase(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        
        #create temporary directory
        self.tmpDirectory = tempfile.TemporaryDirectory()
        self.tmpCacheDirectory = os.path.join(self.tmpDirectory.name,"cache")
        os.mkdir(self.tmpCacheDirectory)
        
        #create config.ini
        self.tmpConfigFile = os.path.join(self.tmpDirectory.name,"config.ini")
        config_file = configparser.ConfigParser()
        config_file["api"] = {
            "port": self._get_free_port(),
            "host": "::",
            "threads": 10,
            "debug": True,
            "proxy_prefix": True
        }
        config_file["settings"] = {
            "data_location_marker": os.path.join(os.path.abspath(os.path.dirname(__file__)),"testdata/marker"),
            "data_location_kmer": os.path.join(os.path.abspath(os.path.dirname(__file__)),"testdata/kmer"),
            "sqlite_db": "db.sqlite"
        }
        config_file["cache"] = {
            "type": "FileSystemCache",
            "dir": self.tmpCacheDirectory,
            "threshold": 100000,
            "timeout": 0
        }
        
        #try to find kmc_analysis (test using library problematic)
#         kmc_query_library = os.path.join(os.path.abspath(os.path.dirname(__file__)),
#                                                       "../../../kmc_analysis/lib/kmc_python.so")     
        kmc_query_binary_location = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                                      "../../../kmc_analysis/bin")
#         if os.path.isfile(kmc_query_library):
#             config_file["settings"]["kmc_query_library"] = kmc_query_library 
#             self.kmcFound = True
        if os.path.isdir(kmc_query_binary_location):
            config_file["settings"]["kmc_query_binary_location"] = kmc_query_binary_location  
            self.kmcFound = True
        else:
            self.kmcFound = False
            
        with open(self.tmpConfigFile,"w") as f:
            config_file.write(f)
        
        #data        
        self.dataLocation = os.path.join(os.path.abspath(os.path.dirname(__file__)),"testdata")        
        #start service
        api = haplotyping.service.API(self.dataLocation, self.tmpConfigFile, False)
        app = api.process_api_messages()
        app.testing = True
        self.client = app.test_client()
        
    def _get_free_port():
        sock = socket.socket()
        sock.bind(('', 0))
        return sock.getsockname()[1]
        
    def _test_paging(self, name, identifier):
        #get total
        response = self.client.get("/api/"+name+"/?start=0&number=0", headers={"accept": "application/json"})    
        self.assertEqual(response.status_code,200,"problem "+name+" number")
        data = json.loads(response.data)
        self.assertEqual(len(data["list"]),0,"unexpected length")
        number = 1
        total = data["total"]
        self.assertTrue(data["total"]>0,"unexpected size (empty)")
        uidSet = set()
        #get countries
        for start in range(0,total,number):
            response = self.client.get("/api/"+name+"/?start="+str(start)+"&number="+str(number), 
                                       headers={"accept": "application/json"})
            self.assertEqual(response.status_code,200,"problem "+name+" result")
            data = json.loads(response.data)
            self.assertEqual(data["start"],start,"incorrect start")
            self.assertEqual(data["number"],number,"incorrect number")
            self.assertTrue(len(data["list"])>0,"expected results")
            for item in data["list"]:
                uidSet.add(item[identifier])
        #check total
        self.assertEqual(len(uidSet),total,"incorrect total number")
        
    def test_api_tools(self):
        kmer = "TCCATCTGTGATAAAGGATCAAGTAAGCCCT"
        #canonical
        response = self.client.get("/api/tools/canonical/"+kmer, headers={"accept": "application/json"})
        self.assertEqual(response.status_code,200,"problem canonical k-mer")
        data = json.loads(response.data)
        self.assertEqual(data,haplotyping.General.canonical(kmer),"incorrect canonical k-mer")
        #reverse complement
        response = self.client.get("/api/tools/reverse-complement/"+kmer, headers={"accept": "application/json"})
        self.assertEqual(response.status_code,200,"problem reverse complement k-mer")
        data = json.loads(response.data)
        self.assertEqual(data,haplotyping.General.reverse_complement(kmer),"incorrect reverse complement k-mer")
        
    def test_api_country(self):
        #test paging
        self._test_paging("country","uid")        

    def test_api_collection(self):
        #test paging
        self._test_paging("collection","name") 
        
    def test_api_variety(self):
        #test paging
        self._test_paging("variety","uid") 
        #get test entries - total number of varieties
        response = self.client.get("/api/variety/?number=0", headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(data["total"]>0,"expected results")
        n_varieties = data["total"]
        #get test entries - variety
        response = self.client.get("/api/variety/?start=0&number="+str(n_varieties), 
                                   headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])>0,"expected results")
        varieties = data["list"]
        #get test entries - country
        response = self.client.get("/api/country/?start=0&number="+str(n_varieties), 
                                   headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])>0,"expected results")
        countries = data["list"]
        #get test entries - collection
        response = self.client.get("/api/collection/?start=0&number="+str(n_varieties), 
                                   headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])>0,"expected results")
        collections = data["list"]
        #get test entries - dataset
        response = self.client.get("/api/dataset/?start=0&number="+str(n_varieties), 
                                   headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])>0,"expected results")
        datasets = data["list"]
        #initialise year examples
        year_examples = {"single":set(),"maximum":set(),"minimum":set(),"range":set()}
        #search by name and check for year examples
        n = 0
        for i in range(len(varieties)):
            response = self.client.get("/api/variety/?"+
                                       urllib.parse.urlencode({"start":0,"number":10,"name":varieties[i]["name"]}), 
                                       headers={"accept": "application/json"})    
            data = json.loads(response.data)
            n+=data["total"]
            for j in range(len(data["list"])):
                if "year" in data["list"][j]:
                    if (data["list"][j]["year"]["min"]==None):
                        year_examples["maximum"].add(data["list"][j]["year"]["max"])
                    elif (data["list"][j]["year"]["max"]==None):
                        year_examples["minimum"].add(data["list"][j]["year"]["min"])
                    elif (data["list"][j]["year"]["min"]==data["list"][j]["year"]["max"]):
                        year_examples["single"].add(data["list"][j]["year"]["min"])
                    else:
                        year_examples["range"].add((data["list"][j]["year"]["min"],data["list"][j]["year"]["max"],))
            self.assertTrue(len(data["list"])==1,"expected results for name "+varieties[i]["name"])
            self.assertTrue(data["list"][0]["name"]==varieties[i]["name"],"wrong name")
        self.assertTrue(n_varieties==n,"incorrect total for varieties searched by name")
        #search by origin
        n = 0
        country_list = []
        for i in range(len(countries)):
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"origin":countries[i]["uid"]}), 
                                   headers={"accept": "application/json"})    
            data = json.loads(response.data)
            n+=data["total"]
            country_list.append(countries[i]["uid"])
            self.assertTrue(len(data["list"])==1,"expected results for origin "+countries[i]["uid"])
            self.assertTrue(data["list"][0]["origin"]["uid"]==countries[i]["uid"],"wrong origin")
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"origin":",".join(country_list)}), 
                                   headers={"accept": "application/json"})   
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])==1,"expected results for multiple origins")
        self.assertTrue(data["total"]==n,"incorrect total for multiple origins: "+",".join(country_list))
        #search by collection
        n = 0
        collection_list = []
        for i in range(len(collections)):
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":collections[i]["uid"]}), 
                                   headers={"accept": "application/json"})    
            data = json.loads(response.data)
            n+=data["total"]
            collection_list.append(collections[i]["uid"])
            self.assertEqual(len(data["list"]),1,"expected results for collection "+collections[i]["name"])
            matches = 0
            for j in range(len(data["list"][0]["datasets"])):
                if(("collection" in data["list"][0]["datasets"][j]) and
                   (data["list"][0]["datasets"][j]["collection"]["uid"]==collections[i]["uid"])) :
                    matches+=1
            self.assertTrue(matches>0,"no datasets with collection "+collections[i]["uid"])
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":",".join(collection_list)}), 
                                   headers={"accept": "application/json"})   
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])==1,"expected results for multiple collections")
        self.assertTrue(data["total"]<=n,"incorrect total for multiple collections: "+",".join(collection_list))
        #search by year
        for year in year_examples["single"]:
            search_year = str(year)
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"year":search_year}), 
                                   headers={"accept": "application/json"})  
            data = json.loads(response.data)
            self.assertTrue(data["total"]>0,"expected results when searching for "+search_year)
            for j in range(len(data["list"])):
                self.assertTrue(not((data["list"][j]["year"]["min"]==None) and 
                                    (data["list"][j]["year"]["max"]==None)) ,
                                "incorrect result when searching for "+search_year)
                self.assertTrue((data["list"][j]["year"]["max"]==None) or 
                                (data["list"][j]["year"]["max"]>=year),
                                "incorrect maximum year when searching for "+search_year)
                self.assertTrue((data["list"][j]["year"]["min"]==None) or 
                                (data["list"][j]["year"]["min"]<=year),
                                "incorrect minimum year when searching for "+search_year)
        for year in year_examples["minimum"]:
            search_year = ">"+str(year-1)
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"year":search_year}), 
                                   headers={"accept": "application/json"})  
            data = json.loads(response.data)
            self.assertTrue(data["total"]>0,"expected results when searching for "+search_year)
            for j in range(len(data["list"])):
                self.assertTrue(not((data["list"][j]["year"]["min"]==None) and 
                                    (data["list"][j]["year"]["min"]==None)) ,
                                "incorrect result when searching for "+search_year)
                self.assertTrue((not (data["list"][j]["year"]["min"]==None)) and 
                                (data["list"][j]["year"]["min"]>=year),
                                "incorrect minimum year when searching for "+search_year)
        for year in year_examples["maximum"]:
            search_year = "<"+str(year+1)
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"year":search_year}), 
                                   headers={"accept": "application/json"})  
            data = json.loads(response.data)
            self.assertTrue(data["total"]>0,"expected results when searching for "+search_year)
            for j in range(len(data["list"])):
                self.assertTrue(not((data["list"][j]["year"]["max"]==None) and 
                                    (data["list"][j]["year"]["min"]==None)) ,
                                "incorrect result when searching for "+search_year)
                self.assertTrue((not (data["list"][j]["year"]["max"]==None)) and 
                                (data["list"][j]["year"]["max"]<=year),
                                "incorrect maximum year when searching for "+search_year)
        for year_range in year_examples["range"]:
            search_year = str(year_range[0])+"-"+str(year_range[1])
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"year":search_year}), 
                                   headers={"accept": "application/json"})  
            data = json.loads(response.data)
            self.assertTrue(data["total"]>0,"expected results when searching for "+search_year)
            for j in range(len(data["list"])):
                self.assertTrue(not((data["list"][j]["year"]["max"]==None) and 
                                    (data["list"][j]["year"]["min"]==None)) ,
                                "incorrect result when searching for "+search_year)
                self.assertTrue((data["list"][j]["year"]["min"]==None) or
                                (data["list"][j]["year"]["min"]>=year_range[0]),
                                "incorrect minimum year when searching for "+search_year)
                self.assertTrue((data["list"][j]["year"]["max"]==None) or
                                (data["list"][j]["year"]["max"]<=year_range[1]),
                                "incorrect maximum year when searching for "+search_year)
        #search by parents       
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"hasParents": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["parents"])>0,
                            "unexpected entry without parents")
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"hasParents": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["parents"])==0,
                            "unexpected entry with parents")
        self.assertEqual(n_parents_true+n_parents_false, n_varieties,
                            "unexpected total for condition on parents")
        #search by offspring     
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"hasOffspring": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["offspring"])>0,
                            "unexpected entry without offspring")
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"hasOffspring": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["offspring"])==0,
                            "unexpected entry with offspring")
        self.assertEqual(n_parents_true+n_parents_false, n_varieties,
                            "unexpected total for condition on offspring")
        #search by dataset
        for dataType in ["any","none","kmer","split","marker"]:
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"dataType": dataType}), 
                                   headers={"accept": "application/json"}) 
            data = json.loads(response.data)
            if(dataType=="any"):
                n_any = data["total"]
            elif(dataType=="none"):
                n_none = data["total"]
            elif(dataType=="kmer"):
                n_kmer = data["total"]
            elif(dataType=="split"):
                n_split = data["total"]
            for i in range(len(data["list"])):
                if(dataType=="any"):
                    self.assertTrue(len(data["list"][i]["datasets"])>0,
                            "unexpected entry with no datasets when dataType condition is {}".format(dataType))
                elif(dataType=="none"):
                    self.assertTrue(len(data["list"][i]["datasets"])==0,
                            "unexpected entry with datasets when dataType condition is {}".format(dataType))
                elif(dataType=="kmer"):
                    for j in range(len(data["list"][i]["datasets"])):
                        self.assertTrue((data["list"][i]["datasets"][j]["type"]=="kmer") or
                                        (data["list"][i]["datasets"][j]["type"]=="split"),
                                "unexpected entry with dataType other then '{}' when condition is {}".format(
                                    dataType,dataType))
                elif(dataType=="split"):
                    for j in range(len(data["list"][i]["datasets"])):
                        self.assertTrue(data["list"][i]["datasets"][j]["type"]==dataType,
                                "unexpected entry with dataType other then '{}' when condition is {}".format(
                                    dataType,dataType))
        self.assertTrue(n_kmer<=n_any,
                            "unexpected number with 'kmer' for condition dataType on datasets")
        self.assertTrue(n_split<=n_any,
                            "unexpected number with 'split' for condition dataType on datasets")
        self.assertEqual(n_any+n_none, n_varieties,
                            "unexpected total for condition dataType on datasets")
        #get variety by uid
        for i in range(len(varieties)):
            response = self.client.get("/api/variety/"+varieties[i]["uid"], 
                                       headers={"accept": "application/json"})  
            self.assertEqual(response.status_code,200,"problem variety "+varieties[i]["uid"])
            data = json.loads(response.data)
            self.assertEqual(varieties[i]["uid"],data["uid"],"unexpected reponse variety "+varieties[i]["uid"])
        #get multiple varieties by uids
        uids = [item["uid"] for item in varieties]
        response = self.client.post("/api/variety/", 
                                       headers={"accept": "application/json"},json={"uids":uids})  
        self.assertEqual(response.status_code,200,"problem getting multiple varieties by uids")
        data = json.loads(response.data)
        self.assertEqual(data["total"], len(uids),
                            "unexpected total when getting multiple varieties by uid")
        self.assertEqual(len(data["list"]), len(uids),
                            "unexpected total when getting multiple varieties by uid")
        #check handling non-existing
        response = self.client.get("/api/variety/nonexistinguid", 
                                       headers={"accept": "application/json"})  
        self.assertEqual(response.status_code,404,"problem non existing uid")
        
        
    def test_api_dataset(self):
        #test paging
        self._test_paging("dataset","uid") 
        #get test entries - total number of datasets
        response = self.client.get("/api/dataset/?number=0", headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(data["total"]>0,"expected results")
        n_datasets = data["total"]
        #get test entries - dataset
        response = self.client.get("/api/dataset/?start=0&number="+str(n_datasets), 
                                   headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])>0,"expected results")
        datasets = data["list"]
        #get test entries - collection
        response = self.client.get("/api/collection/?start=0&number="+str(n_datasets), 
                                   headers={"accept": "application/json"})    
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])>0,"expected results")
        collections = data["list"]
        #search by collection
        n = 0
        collection_list = []
        for i in range(len(collections)):
            response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":collections[i]["uid"]}), 
                                   headers={"accept": "application/json"})    
            data = json.loads(response.data)
            n+=data["total"]
            if data["total"]>0:
                collection_list.append(collections[i]["name"])
                self.assertTrue(len(data["list"])==1,"expected results for collection {}".format(collections[i]["uid"]))
                self.assertTrue("collection" in data["list"][0],"no collection for dataset {}".format(collections[i]["uid"]))
                self.assertEqual(data["list"][0]["collection"]["uid"],collections[i]["uid"],
                                 "no dataset with collection {}".format(collections[i]["uid"]))
        response = self.client.get("/api/collection/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":",".join(collection_list)}), 
                                   headers={"accept": "application/json"})   
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])==1,"expected results for multiple collections")
        self.assertTrue(data["total"]<=n,"incorrect total for multiple collections: "+",".join(collection_list))
        #search by variety 
        response = self.client.get("/api/collection/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"hasVariety": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_variety_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue((not "variety" in data["list"][j]) or (data["list"][j]["variety"]==None),
                            "unexpected entry without variety")
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"hasVariety": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_variety_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue((not "variety" in data["list"][j]["variety"]) or (data["list"][j]["variety"]==None),
                            "unexpected entry with variety")
        self.assertEqual(n_variety_true+n_variety_false, n_datasets,
                            "unexpected total for condition on variety {}".format(urllib.parse.urlencode({"start":0,"number":n_datasets,"variety": False})))
        #get dataset by uid
        for i in range(len(datasets)):
            response = self.client.get("/api/dataset/"+datasets[i]["uid"], 
                                       headers={"accept": "application/json"})  
            self.assertEqual(response.status_code,200,"problem dataset "+datasets[i]["uid"])
            data = json.loads(response.data)
            self.assertEqual(datasets[i]["uid"],data["uid"],"unexpected reponse dataset "+datasets[i]["uid"])
        #get multiple datasets by uids
        uids = [item["uid"] for item in datasets]
        response = self.client.post("/api/dataset/", 
                                       headers={"accept": "application/json"},json={"uids":uids})  
        self.assertEqual(response.status_code,200,"problem getting multiple datasets by uids")
        data = json.loads(response.data)
        self.assertEqual(data["total"], len(uids),
                            "unexpected total when getting multiple datasets by uid")
        self.assertEqual(len(data["list"]), len(uids),
                            "unexpected total when getting multiple datasets by uid")
        #check handling non-existing
        response = self.client.get("/api/dataset/nonexistinguid", 
                                       headers={"accept": "application/json"})  
        self.assertEqual(response.status_code,404,"problem non existing uid")

    def test_html(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code,200,"site not available")
        self.assertTrue(b"<html" in response.data and b"</html>" in response.data,"no html")
        
    def test_api_kmer(self):
        #get dataset uid
        response = self.client.get("/api/dataset/?hasVariety=true&dataType=split", headers={"accept": "application/json"})  
        self.assertEqual(response.status_code,200,"get dataset with kmer database")
        data = json.loads(response.data)
        self.assertTrue(("list" in data) and (len(data["list"])>0),"no dataset with kmer database found")
        uid = data["list"][0]["uid"]
        #first check for kmc
        self.assertTrue(self.kmcFound,"kmc not found, couldn't do these tests")        
        #get frequency with mismatches (specific tests for this dataset, update if changed!)
        kmer1 = "AGTATATACATCAACTACGAAAATGGAAAAG"
        freq1 = 11
        kmer2 = "AGTATATACATCAACTAGGAAAATGGAAAAG"
        freq2 = 1
        response = self.client.get("/api/kmer/{}/{}?mismatches=2".format(uid,kmer1), 
                                   headers={"accept": "application/json"}) 
        self.assertEqual(response.status_code,200,"get k-mer frequency")
        data = json.loads(response.data)
        self.assertTrue(("kmers" in data) and (kmer1 in data["kmers"]),
                        "k-mer {} not found in result".format(kmer1))
        self.assertEqual(data["kmers"][kmer1],freq1,"frequency {} expected to be {}".format(kmer1,freq1))
        self.assertTrue(kmer2 in data["kmers"], "mismatch k-mer {} not found in result".format(kmer1))
        self.assertEqual(data["kmers"][kmer2],freq2,"frequency {} expected to be {}".format(kmer1,freq1))
        #get frequency for invalid k-mer
        invalidKmer = "AAAA"
        response = self.client.get("/api/kmer/{}/{}?mismatches=2".format(uid,invalidKmer), 
                                   headers={"accept": "application/json"}) 
        self.assertEqual(response.status_code,200,"get k-mer frequency")
        data = json.loads(response.data)
        self.assertTrue(("kmers" in data) and (len(data["kmers"])==0), "no k-mers expected")
        #get multiple frequencies
        response = self.client.post("/api/kmer/{}".format(uid), 
                                   headers={"accept": "application/json"},json={"kmers":[kmer1,kmer2]}) 
        self.assertEqual(response.status_code,200,"problem getting k-mers")
        data = json.loads(response.data)
        self.assertTrue("kmers" in data,"k-mers not found in result")
        self.assertTrue(kmer1 in data["kmers"], "k-mer {} not found in result".format(kmer1))
        self.assertEqual(data["kmers"][kmer1],freq1,"frequency {} expected to be {}".format(kmer1,freq1))
        self.assertTrue(kmer2 in data["kmers"], "k-mer {} not found in result".format(kmer2))
        self.assertEqual(data["kmers"][kmer2],freq2,"frequency {} expected to be {}".format(kmer1,freq1))
        #get k-mers from sequence
        sequence = "TATCAACCAAAAATAACCTTTCTTTGAAAACAAGTAATAAGATATGAAGATAAATGACAAGAGTGATATTCTGTATATACTTACAGCT"
        response = self.client.post("/api/kmer/{}/sequence".format(uid), 
                                   headers={"accept": "application/json"},json={"sequence":sequence})  
        self.assertEqual(response.status_code,200,"problem getting k-mers from sequence")
        data = json.loads(response.data)
        k = data["info"]["kmer_length"]
        self.assertTrue(("kmers" in data) and (len(data["kmers"])>0),"no k-mers in response")
        self.assertEqual(len(data["kmers"]),1+len(sequence)-k,"unexpected number of k-mers")
        #get path
        kmer1 = sequence[0:k]
        kmer2 = sequence[-k:]
        response = self.client.get("/api/kmer/{}/{}/path/{}".format(uid,kmer1,kmer2), 
                                   headers={"accept": "application/json"}) 
        self.assertEqual(response.status_code,200,"get k-mer path {} - {}".format(kmer1,kmer2))
        data = json.loads(response.data)
        self.assertEqual(len(data),len(sequence),"unexpected path size")
        #find nearest split
        kmer = "AGTATATACATCAACTACGAAAATGGAAAAG"
        response = self.client.get("/api/kmer/{}/{}/split".format(uid,kmer), 
                                   headers={"accept": "application/json"}) 
        self.assertEqual(response.status_code,200,"get k-mer nearest split {}".format(kmer))    
        data = json.loads(response.data)
        self.assertTrue(("distance" in data) and (data["distance"]>0),"distance not found in result")
        self.assertTrue(("splittingKmers" in data) and (len(data["splittingKmers"])>0),"splitting k-mers not found in result")
        self.assertTrue(("pathKmers" in data) and (len(data["pathKmers"])==data["distance"]),"path k-mers not found in result")
            
            
        
        

    @classmethod
    def tearDownClass(self):
        if self.tmpDirectory:
            self.tmpDirectory.cleanup()    
     
        
        
        
        