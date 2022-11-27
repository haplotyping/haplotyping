import logging,requests
import haplotyping.graph
import haplotyping.graph.baseGraph as baseGraph
    
class APIGraph(baseGraph.Graph):
    """basic or minimal version of the De Bruijn graph using the dataset 
        defined by uid via API"""
    
    def __init__(self, url: str, uid: str, username: str = None, password: str = None):
        """use the API at provided url"""
        #call parent constructor
        super(haplotyping.graph.api.APIGraph, self).__init__()
        #logger
        self._api_logger = logging.getLogger(__name__)
        #store call parameters
        self._baseUrl : str = str(url)
        self._uid : str = str(uid)
        self._username : str = username
        self._password : str = password
        #initialise
        self._variety : str = None
        self._collection : str = None
        #report to logger
        self._api_logger.debug("using dataset {}".format(self._uid))
        self._api_logger.debug("using API at {}".format(self._baseUrl))
        #set authentication
        if self._username and self._password:
            self._apiAuth = requests.auth.HTTPBasicAuth(self._username, self._password)
        else:
            self._apiAuth = None
        #set headers
        self._apiHeaders = {"accept": "application/json"}
        #checks
        self._api_initial_checks()
        
    def _api_initial_checks(self):
        """
        get dataset and check for k-mer and split
        """
        self._api_logger.debug("get dataset information from API")        
        #get dataset info
        response = requests.post(self._baseUrl+"dataset/", json = {"uids": [self._uid]}, 
                                 auth=self._apiAuth, headers=self._apiHeaders)
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
                    self._collection = str(dataset["collection"]["name"])
                    self._api_logger.debug("set collection to '{}'".format(self._collection))
                if self._variety:
                    self._name = "'{}'".format(self._variety)
                    if self._collection:
                        self._name = "{} from '{}' collection".format(self._name,self._collection)
                #get split info
                response = requests.get(self._baseUrl+"split/"+self._uid+"/info", 
                                        auth=self._apiAuth, headers=self._apiHeaders)
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
                                 json = {"kmers": sorted(kmerList)}, 
                                 auth=self._apiAuth, headers=self._apiHeaders)
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
            kmerList = sorted(rightNeighbours.union(leftSplitters))
            response = requests.post(self._baseUrl+"kmer/"+self._uid, 
                                 json = {"kmers": sorted(kmerList)}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                kmerFound = [k for k in data["kmers"] if data["kmers"][k]>=minimumFrequency]
                if len(leftSplitters.intersection(kmerFound))>0:
                    return [newKmer]
                else:
                    rightSplitters = sorted(rightNeighbours.intersection(kmerFound))
                    if len(rightSplitters)>1:
                        return sorted(rightSplitters)
                    elif len(rightSplitters)==0:
                        return None
                    else:
                        newKmer = rightSplitters[0]   
                        print(kmerFound)
            else:
                self._logger.error("request to {} didn't succeed".format(self._baseUrl))
        
      