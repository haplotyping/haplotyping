import logging,requests
import haplotyping.graph
import haplotyping.graph.baseGraph as baseGraph

class API():
    def __init__(self, url: str):
        #logger
        self._api_logger = logging.getLogger(__name__)
        #store call parameters
        self._baseUrl : str = str(url)
            
    def _getList(self, request:str, key: str, number: int = 1000):
        start = 0
        result = {}
        while True:
            fullRequest = "{}{}{}start={}&number={}".format(self._baseUrl,
                request,"&" if "?" in request else "?",start,number)
            response = requests.get(fullRequest)
            if response.ok:
                data = response.json()
                total = data.get("total",None)
                items = data.get("list",[])
                for item in items:
                    uid = item.get(key,None)
                    if uid:
                        del item[key]
                        result[uid] = item
                if start+number>=total:
                    return result
                else:
                    start+=number
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None
    
    def _getById(self, request, uid, key="uid"):
        fullRequest = "{}{}{}".format(self._baseUrl,request,requests.utils.quote(uid))
        response = requests.get(fullRequest)
        if response.ok:
            data = response.json()
            if key in data:
                del data[key]
                return data
            else:
                return None
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def _getByIds(self, request, uids, key1="uids", key2="uid"):
        fullRequest = "{}{}".format(self._baseUrl,request)
        response = requests.post(fullRequest, json = {key1: list(uids)})
        if response.ok:
            data = response.json()
            result = {}
            items = data.get("list",[])
            for item in items:
                uid = item.get(key2,None)
                if uid:
                    del item[key2]
                    result[uid] = item
            return result
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    
    """
    get a dict with all countries
    """
    def getCountries(self): 
        return self._getList("country/","uid")
        
                
    """
    get a dict with all collections
    """
    def getCollections(self):
        return self._getList("collection/","uid")
    
    """
    get a dict with (all) datasets
    """
    def getDatasets(self, **kwargs):
        parameters = {
            "collection": None,
            "dataType": None,
            "hasVariety": None
        }
        request = "dataset/"
        for key, value in kwargs.items():
            if key in parameters.keys():
                if isinstance(value,bool):
                    parameters[key] = "true" if value else "false"
                elif isinstance(value,list):
                    parameters[key] = ",".join(value)
                else:
                    parameters[key] = value
            if not parameters[key]==None:
                request = "{}{}{}={}".format(request,
                         "&" if "?" in request else "?",key,requests.utils.quote(parameters[key]))
        return self._getList(request,"uid")
    
    def getDatasetById(self, uid):
        return self._getById("dataset/",uid)
    
    def getDatasetsById(self, uids):
        return self._getByIds("dataset/",uids)
    
    """
    get a dict with (all) varieties
    """
    def getVarieties(self, **kwargs):
        parameters = {
            "name": None,
            "origin": None,
            "year": None,
            "collection": None,
            "dataType": None,
            "hasParents": None,
            "hasOffspring": None
        }
        request = "variety/"
        for key, value in kwargs.items():
            if key in parameters.keys():
                if isinstance(value,bool):
                    parameters[key] = "true" if value else "false"
                elif isinstance(value,list):
                    parameters[key] = ",".join(value)
                else:
                    parameters[key] = value
            if not parameters[key]==None:
                request = "{}{}{}={}".format(request,
                         "&" if "?" in request else "?",key,requests.utils.quote(parameters[key]))
        print(request)
        return self._getList(request,"uid")
    
    def getVarietyById(self, uid):
        return self._getById("variety/",uid)
    
    def getVarietiesById(self, uids):
        return self._getByIds("variety/",uids)
    
    def getKmerFrequency(self, datasetUid, kmer, mismatches=0):
        fullRequest = "{}kmer/{}/{}{}".format(self._baseUrl, datasetUid, kmer, 
                                              "?mismatches={}".format(int(mismatches)) if mismatches>0 else "")
        response = requests.get(fullRequest)
        if response.ok:
            data = response.json()
            result = data.get("kmers",{})
            if not kmer in result.keys():
                result[kmer] = 0
            return result
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
    
    def getKmerFrequencies(self, datasetUid, kmers, mismatches=0):
        fullRequest = "{}kmer/{}".format(self._baseUrl, datasetUid)        
        response = requests.post(fullRequest, json={"kmers": list(kmers),"mismatches": mismatches})
        if response.ok:
            data = response.json()
            result = data.get("kmers",{})
            for kmer in kmers:
                if not kmer in result.keys():
                    result[kmer] = 0
            return result
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getKmerSequence(self, datasetUid, sequence, mismatches=0):
        fullRequest = "{}kmer/{}/sequence".format(self._baseUrl, datasetUid)        
        response = requests.post(fullRequest, json={"sequence": sequence,"mismatches": mismatches})
        if response.ok:
            data = response.json()
            result = data.get("kmers",{})
            return result
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerData(self, datasetUid):
        fullRequest = "{}marker/{}/data".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerInfo(self, datasetUid):
        fullRequest = "{}marker/{}/info".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerMapping(self, datasetUid):
        fullRequest = "{}marker/{}/mapping".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
    
    
class APIGraph(baseGraph.Graph):
    """basic or minimal version of the De Bruijn graph using the dataset 
        defined by uid via API"""
    
    def __init__(self, url: str, uid: str):
        """use the API at provided url"""
        #call parent constructor
        super(haplotyping.graph.api.APIGraph, self).__init__()
        #logger
        self._api_logger = logging.getLogger(__name__)
        #store call parameters
        self._baseUrl : str = str(url)
        self._uid : str = str(uid)
        #initialise
        self._variety : str = None
        self._collection : str = None
        #report to logger
        self._api_logger.debug("using dataset {}".format(self._uid))
        self._api_logger.debug("using API at {}".format(self._baseUrl))
        #checks
        self._api_initial_checks()
        
    def _api_initial_checks(self):
        """
        get dataset and check for k-mer and split
        """
        self._api_logger.debug("get dataset information from API")        
        #get dataset info
        response = requests.post(self._baseUrl+"dataset/", json = {"uids": [self._uid]})
        if response.ok:
            data = response.json()
            if data.get("total",0)==1:
                dataset = data["list"][0]
                if not (dataset["type"]=="kmer" or dataset["type"]=="split"):
                    self._api_logger.error("dataset {} has no k-mer database".format(self._uid))
                if not (dataset["type"]=="split"):
                    self._api_logger.error("dataset {} has no splitting k-mer database".format(self._uid))   
                if dataset["variety"]:
                    self._variety = str(dataset["variety"]["name"])
                    self._api_logger.debug("set variety to '{}'".format(self._variety))
                if dataset["collection"]:
                    self._collection = str(dataset["collection"])
                    self._api_logger.debug("set collection to '{}'".format(self._collection))
                #get split info
                response = requests.get(self._baseUrl+"split/"+self._uid+"/info")
                if response.ok:
                    data = response.json()
                    self._k = int(data["k"])
                    self._api_logger.debug("set k-mer size to {}".format(self._k))
                else:
                    self._api_logger.error("request to {} didn't succeed".format(self._baseUrl)) 
            else:
                self._api_logger.error("dataset {} not found".format(self._uid))   
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            
    def getRightNeighbours(self, kmer: str):
        kmerList = []
        for b in ["A","C","G","T"]:
            kmerList.append(kmer[1:]+b)
        response = requests.post(self._baseUrl+"kmer/"+self._uid, 
                                 json = {"kmers": sorted(kmerList)})
        if response.ok:
            data = response.json()
            return data
        else:
            self._logger.error("request to "+self._baseUrl+" didn't succeed")
            
    def findRightConnection(self, kmer: str, maximumDistance: int = 100, minimumFrequency: int = 2):
        newKmer = kmer
        for i in range(maximumDistance):
            rightNeighbours = set()
            leftSplitters = set()
            for rb in ["A","C","G","T"]:
                rightNeighbours.add(newKmer[1:]+rb)
            if i>0:
                for lb in ["A","C","G","T"]:
                    if not lb==newKmer[0]:
                        leftSplitters.add(lb+newKmer[1:])
            kmerList = list(rightNeighbours.union(leftSplitters))
            response = requests.post(self._baseUrl+"kmer/"+self._uid, 
                                 json = {"kmers": sorted(kmerList)})
            if response.ok:
                data = response.json()
                kmerFound = [k for k in data["kmers"] if data["kmers"][k]>=minimumFrequency]
                if len(leftSplitters.intersection(kmerFound))>0:
                    return [newKmer]
                else:
                    rightSplitters = list(rightNeighbours.intersection(kmerFound))
                    if len(rightSplitters)>1:
                        return sorted(rightSplitters)
                    elif len(rightSplitters)==0:
                        return None
                    else:
                        newKmer = rightSplitters[0]   
                        print(kmerFound)
            else:
                self._logger.error("request to {} didn't succeed".format(self._baseUrl))
        
      