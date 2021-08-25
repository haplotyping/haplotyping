import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import unittest, json, pathlib, urllib

import haplotyping, haplotyping.service

class APITestCase(unittest.TestCase):
    
    def setUp(self):
        location = str(pathlib.Path().absolute())
        location = location + pathlib.os.sep + "tests" + pathlib.os.sep + "service"
        api = haplotyping.service.API(location, False)
        app = api.process_api_messages(False)
        app.config["TESTING"] = True
        self.client = app.test_client()
        
    def _test_paging(self, name, identifier):
        #get total
        response = self.client.get("/api/"+name+"/?start=0&number=0", headers={"accept": "application/json"})    
        self.assertEqual(response.status_code,200,"problem "+name+" number")
        data = json.loads(response.data)
        self.assertEqual(len(data["list"]),0,"unexpected length")
        number = 1
        total = data["total"]
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
                if (data["list"][j]["year"]["min"]==None) and (data["list"][j]["year"]["max"]==None):
                    pass
                elif (data["list"][j]["year"]["min"]==None):
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
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":collections[i]["name"]}), 
                                   headers={"accept": "application/json"})    
            data = json.loads(response.data)
            n+=data["total"]
            collection_list.append(collections[i]["name"])
            self.assertTrue(len(data["list"])==1,"expected results for collection "+collections[i]["name"])
            matches = 0
            for j in range(len(data["list"][0]["datasets"])):
                if(data["list"][0]["datasets"][j]["collection"]==collections[i]["name"]) :
                    matches+=1
            self.assertTrue(matches>0,"no datasets with collection "+collections[i]["name"])
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
                                    (data["list"][j]["year"]["min"]==None)) ,
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
                                    (data["list"][j]["year"]["max"]==None)) ,
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
                                    (data["list"][j]["year"]["max"]==None)) ,
                                "incorrect result when searching for "+search_year)
                self.assertTrue((data["list"][j]["year"]["min"]==None) or
                                (data["list"][j]["year"]["min"]>=year_range[0]),
                                "incorrect minimum year when searching for "+search_year)
                self.assertTrue((data["list"][j]["year"]["max"]==None) or
                                (data["list"][j]["year"]["max"]<=year_range[1]),
                                "incorrect maximum year when searching for "+search_year)
        #search by parents       
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"parents": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["parents"])>0,
                            "unexpected entry without parents")
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"parents": False}), 
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
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"offspring": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["offspring"])>0,
                            "unexpected entry without offspring")
        response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"offspring": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_parents_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(len(data["list"][j]["offspring"])==0,
                            "unexpected entry with offspring")
        self.assertEqual(n_parents_true+n_parents_false, n_varieties,
                            "unexpected total for condition on offspring")
        #search by dataset
        for dataset in ["any","none","kmer","split"]:
            response = self.client.get("/api/variety/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_varieties,"dataset": dataset}), 
                                   headers={"accept": "application/json"}) 
            data = json.loads(response.data)
            if(dataset=="any"):
                n_any = data["total"]
            elif(dataset=="none"):
                n_none = data["total"]
            elif(dataset=="kmer"):
                n_kmer = data["total"]
            elif(dataset=="split"):
                n_split = data["total"]
            for i in range(len(data["list"])):
                if(dataset=="any"):
                    self.assertTrue(len(data["list"][i]["datasets"])>0,
                            "unexpected entry with no datasets when dataset condition is "+dataset)
                elif(dataset=="none"):
                    self.assertTrue(len(data["list"][i]["datasets"])==0,
                            "unexpected entry with datasets when dataset condition is "+dataset)
                elif(dataset=="kmer"):
                    for j in range(len(data["list"][i]["datasets"])):
                        self.assertTrue(not (data["list"][i]["datasets"][j]["kmer"]==None),
                                "unexpected entry with dataset without 'kmer' when dataset condition is "+dataset)
                elif(dataset=="split"):
                    for j in range(len(data["list"][i]["datasets"])):
                        self.assertTrue(not (data["list"][i]["datasets"][j]["split"]==None),
                                "unexpected entry with dataset without 'split' when dataset condition is "+dataset)
        self.assertTrue(n_kmer<=n_any,
                            "unexpected number with 'kmer' for condition on datasets")
        self.assertTrue(n_split<=n_any,
                            "unexpected number with 'split' for condition on datasets")
        self.assertEqual(n_any+n_none, n_varieties,
                            "unexpected total for condition on datasets")
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
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":collections[i]["name"]}), 
                                   headers={"accept": "application/json"})    
            data = json.loads(response.data)
            n+=data["total"]
            collection_list.append(collections[i]["name"])
            self.assertTrue(len(data["list"])==1,"expected results for collection "+collections[i]["name"])
            self.assertEqual(data["list"][0]["collection"],collections[i]["name"],
                             "no dataset with collection "+collections[i]["name"])
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":1,"collection":",".join(collection_list)}), 
                                   headers={"accept": "application/json"})   
        data = json.loads(response.data)
        self.assertTrue(len(data["list"])==1,"expected results for multiple collections")
        self.assertTrue(data["total"]<=n,"incorrect total for multiple collections: "+",".join(collection_list))
        #search by variety 
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"variety": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_variety_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(not(data["list"][j]["variety"]==None),
                            "unexpected entry without variety")
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"variety": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_variety_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(data["list"][j]["variety"]==None,
                            "unexpected entry with variety")
        self.assertEqual(n_variety_true+n_variety_false, n_datasets,
                            "unexpected total for condition on variety")
        #search by kmer
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"kmer": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_kmer_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(data["list"][j]["kmer"],
                            "unexpected entry without kmer")
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"kmer": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_kmer_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(not data["list"][j]["kmer"],
                            "unexpected entry with kmer")
        self.assertEqual(n_kmer_true+n_kmer_false, n_datasets,
                            "unexpected total for condition on kmer")
        #search by split
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"split": True}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_split_true = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(data["list"][j]["split"],
                            "unexpected entry without split")
        response = self.client.get("/api/dataset/?"+
                                   urllib.parse.urlencode({"start":0,"number":n_datasets,"split": False}), 
                                   headers={"accept": "application/json"}) 
        data = json.loads(response.data)
        n_split_false = data["total"]
        for j in range(len(data["list"])):
            self.assertTrue(not data["list"][j]["split"],
                            "unexpected entry with split")
        self.assertEqual(n_split_true+n_split_false, n_datasets,
                            "unexpected total for condition on split")
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
     
        
        
        
        