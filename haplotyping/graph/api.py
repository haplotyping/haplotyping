import logging,requests
import haplotyping

class APIGraph(haplotyping.baseGraph.Graph):
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
        self._api_logger.debug("using dataset "+str(self._uid))
        self._api_logger.debug("using API at "+str(self._baseUrl))
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
                if not dataset["kmer"]:
                    self._api_logger.error("dataset "+self._uid+" has no k-mer database")
                if not dataset["split"]:
                    self._api_logger.error("dataset "+self._uid+" has no split k-mer database")   
                if dataset["variety"]:
                    self._variety = str(dataset["variety"]["name"])
                    self._api_logger.debug("set variety to '%s'" % (self._name))
                if dataset["collection"]:
                    self._collection = str(dataset["collection"])
                    self._api_logger.debug("set collection to '%s'" % (self._collection))
                #get split info
                response = requests.get(self._baseUrl+"split/"+self._uid+"/info")
                if response.ok:
                    data = response.json()
                    self._k = int(data["k"])
                    self._api_logger.debug("set k-mer size to %d" % (self._k))
                else:
                    self._api_logger.error("request to "+self._baseUrl+" didn't succeed") 
            else:
                self._api_logger.error("dataset "+self._uid+" not found")   
        else:
            self._api_logger.error("request to "+self._baseUrl+" didn't succeed")
        
        
    
    