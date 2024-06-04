import logging,requests
import haplotyping.graph
import haplotyping.api
import haplotyping.graph.baseGraph as baseGraph
import pandas as pd
    
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
        self._varietyUid : str = None
        self._username : str = username
        self._password : str = password
        #initialise
        self._variety : str = None
        self._collection : str = None
        self._datasetFrequencyCkmers = set()
        self._datasetFrequencies = {}
        self._datasetVarieties = {}
        self._datasetUids = set()
        #report to logger
        self._api_logger.debug("using dataset {}".format(self._uid))
        self._api_logger.debug("using API at {}".format(self._baseUrl))
        #api
        self._api = haplotyping.api.API(url, username=username, password=password)        
        #checks
        self._api_initial_checks()
        
    def _api_initial_checks(self):
        """
        get dataset and check for k-mer and split
        """
        self._api_logger.debug("get dataset information from API")        
        #get dataset info
        try:
            dataset = self._api.getDatasetById(self._uid)
            if not (dataset["type"]=="kmer" or dataset["type"]=="split"):
                self._api_logger.error("dataset {} has no k-mer database".format(self._uid))
            if not (dataset["type"]=="split"):
                self._api_logger.error("dataset {} has no splitting k-mer database".format(self._uid))   
            if dataset["variety"]:
                self._variety = str(dataset["variety"]["name"])
                self._varietyUid = str(dataset["variety"]["uid"])
                self._api_logger.debug("set variety to '{}'".format(self._variety))
            if dataset["collection"]:
                self._collection = str(dataset["collection"]["name"])
                self._api_logger.debug("set collection to '{}'".format(self._collection))
            if self._variety:
                self._name = "'{}'".format(self._variety)
                if self._collection:
                    self._name = "{} from '{}' collection".format(self._name,self._collection)
            #get split info
            try:
                data = self._api.getSplitInfo(self._uid)
                self._k = int(data["k"])
                self._api_logger.debug("set k-mer size to {}".format(self._k))
            except:
                self._api_logger.error("k-mer database {} not found".format(self._uid)) 
        except:
            self._api_logger.error("dataset {} not found".format(self._uid))
                    
    def getDatasetVarieties(self):
        return self._datasetVarieties
    
    def getDatasetFrequencies(self, datasetUids):
        #compute dataset uids
        if datasetUids==None:
            datasetUids = set(self._api.getDatasets(hasVariety=True, dataType="kmer").keys())
        candidateAndConnectionKmers = [k[0] for k in self.getCandidates()]
        for connection in self._connections:
            candidateAndConnectionKmers = candidateAndConnectionKmers + [k[0] for k in connection.getOrientatedCkmers()]
        recompute = False
        if len(set(datasetUids).difference(self._datasetUids))>0:
            recompute = True
        elif len(set(candidateAndConnectionKmers).difference(self._datasetFrequencyCkmers))>0:
            recompute = True
        if recompute:
            self._logger.debug("recompute frequency for {} k-mers in {} datasets".format(
                len(candidateAndConnectionKmers),len(datasetUids)))
            self._datasetFrequencyCkmers = set(candidateAndConnectionKmers)
            self._datasetFrequencies = {}
            self._datasetVarieties = {}
            self._datasetUids = set()
            #get dataset info
            data = self._api.getDatasetById(datasetUids)
            for dsUid,item in data.items():
                if "variety" in item.keys():
                    self._datasetVarieties[dsUid] = item["variety"]
            #get frequencies
            for dsUid in datasetUids:
                data = self._api.getKmer(dsUid,candidateAndConnectionKmers, mismatches=0)
                self._datasetFrequencies[dsUid] = data.get("kmers",{})
                nfound = len([x for x in self._datasetFrequencies[dsUid].values() if x>0])
                name = self._datasetVarieties.get(dsUid,{"name": None})["name"]
                for k in candidateAndConnectionKmers:
                    if not k in self._datasetFrequencies[dsUid].keys():
                        self._datasetFrequencies[dsUid][k] = 0
                self._logger.debug("get frequency for {} of {} k-mers in {}{}".format(
                    nfound,len(candidateAndConnectionKmers),dsUid," ({})".format(name) if name else ""))
            #register datasets
            self._datasetUids.update(datasetUids)
        #return result
        return pd.DataFrame(self._datasetFrequencies).transpose()


