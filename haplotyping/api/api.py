import logging,requests, html
import re, haplotyping.general
from graphviz import Digraph

class API():
    def __init__(self, url: str, username: str = None, password: str = None):
        #logger
        self._api_logger = logging.getLogger(__name__)
        #store call parameters
        self._baseUrl : str = str(url)
        self._username = username
        self._password = password
        #set authentication
        if self._username and self._password:
            self._apiAuth = requests.auth.HTTPBasicAuth(self._username, self._password)
        else:
            self._apiAuth = None
        #set headers
        self._apiHeaders = {"accept": "application/json"}
        
    def _addParametersToRequest(self, request, parameters, kwargs):
        for key, value in kwargs.items():
            if key in parameters.keys():
                if isinstance(value,bool):
                    parameters[key] = "true" if value else "false"
                elif isinstance(value,list):
                    parameters[key] = ",".join(value)
                else:
                    parameters[key] = str(value)
            if not parameters[key]==None:
                request = "{}{}{}={}".format(request,
                         "&" if "?" in request else "?",key,requests.utils.quote(parameters[key]))
        return request
    
    def _getList(self, request:str, key: str, number: int = 1000):
        start = 0
        result = {}
        while True:
            fullRequest = "{}{}{}start={}&number={}".format(self._baseUrl,
                request,"&" if "?" in request else "?",start,number)
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
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
        fullRequest = "{}{}{}".format(self._baseUrl,request,requests.utils.quote(str(uid)))
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
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
        response = requests.post(fullRequest, json = {key1: list(uids)}, auth=self._apiAuth, headers=self._apiHeaders)
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
       
    # ----------------
    
    def getToolsCanonical(self, kmer: str):
        """
            get canonical representation of the k-mer
            preferably use the local haplotyping.General.canonical method
            
            :param str kmer: k-mer
            :return: canonical
            :rtype: str
        """
        fullRequest = "{}{}{}".format(self._baseUrl,"tools/canonical/",requests.utils.quote(str(kmer)))
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getToolsReverseComplement(self, kmer: str):
        """
            get the reverse complement of the k-mer
            preferably use the local haplotyping.General.reverse_complement method
        
            :param str kmer: k-mer
            :return: reverse complement
            :rtype: str
        """
        fullRequest = "{}{}{}".format(self._baseUrl,"tools/reverse-complement/",requests.utils.quote(kmer))
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    # ----------------
    
    def getCountries(self): 
        """
            get a list with all countries
        
            :return: all countries
            :rtype: list
        """
        return self._getList("country/","uid")
        
    # ----------------
    
    def getCollections(self):
        """
            get all collections
        
            :return: a list of collections
            :rtype: list
        """
        return self._getList("collection/","uid")
    
    def getCollectionById(self, collectionUid: str):
        """
            get a collection by uid
        
            :param str collectionUid: the collection uid
            :return: collection
            :rtype: dict
        """
        if not isinstance(datasetUid,str):
            self._api_logger.error("requesting multiple collections not supported")
        else:
            return self._getById("collection/",collectionUid)
    
    # ----------------
    
    def getDatasets(self, **kwargs):
        """
            get a list with (all) datasets

            :param str collection: optional, comma separated list of collection uids
            :param str dataType: optional, one of "marker", "kmer", "split"
            :param bool hasVariety: optional, has a linked variety
            :return: (all) datasets
            :rtype: list
        """
        parameters = {
            "collection": None,
            "dataType": None,
            "hasVariety": None
        }
        request = "dataset/"
        request = self._addParametersToRequest(request, parameters, kwargs)
        return self._getList(request,"uid")
    
    def getDatasetById(self, datasetUid):
        """
            get one or multiple datasets by uid
        
            :param Union[str,list] datasetUid: dataset uid(s)
            :return: dataset(s)
            :rtype: dict
        """
        if not isinstance(datasetUid,str):
            return self._getByIds("dataset/",datasetUid)
        else:
            return self._getById("dataset/",datasetUid)
    
    # ----------------
    
    def getVarieties(self, **kwargs):
        """
            get a dict with (all) varieties

            :param str name: optional, name or synonym
            :param str nameContains: optional, part of name or synonym
            :param str origin: optional, comma separated list of country codes
            :param str year: optional, year of variety (e.g. '1995', '<1995', '>1995', '1990-1995')
            :param str collection: optional, comma separated list of collection uids
            :param str dataType: optional, one of "marker", "kmer", "split"
            :param bool hasParents: optional, has defined parents
            :param bool hasOffspring: optional, has defined offspring
            :return: (all) varieties
            :rtype: dict
        """
        parameters = {
            "name": None,
            "nameContains": None,
            "origin": None,
            "year": None,
            "collection": None,
            "dataType": None,
            "hasParents": None,
            "hasOffspring": None
        }
        request = "variety/"
        request = self._addParametersToRequest(request, parameters, kwargs)
        return self._getList(request,"uid")
    
    def getVarietyById(self, varietyUid):
        """
            get one or multiple varieties by uid
        
            :param Union[str,list] varietyUid: variety uid(s)
            :return: varieties
            :rtype: dict
        """
        if not isinstance(varietyUid,str):
            return self._getByIds("variety/",varietyUid)
        else:
            return self._getById("variety/",varietyUid)
    
    # ----------------
    
    def getKmer(self, datasetUid, kmer, mismatches=0):
        """
            get k-mer(s)
        
            :param str datasetUid: dataset uid
            :param Union[str,list] kmer: one or multiple k-mers
            :param int mismatches: optional, number of mismatches
            :return: k-mer data
            :rtype: dict
        """
        if not isinstance(mismatches,int) or (mismatches<0):
            mismatches = 0
        if not isinstance(kmer,str):
            fullRequest = "{}kmer/{}".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer),"mismatches": mismatches}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}kmer/{}/{}{}".format(self._baseUrl, datasetUid, kmer, 
                                              "?mismatches={}".format(int(mismatches)) if mismatches>0 else "")
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getKmerFrequency(self, datasetUid, kmer, mismatches=0):
        """
            get k-mer frequency (based on output getKmer)
        
            :param str datasetUid: dataset uid
            :param Union[str,list] kmer: one or multiple k-mers
            :param int mismatches: optional, number of mismatches
            :return: k-mer frequencies
            :rtype: dict
        """
        data = self.getKmer(datasetUid, kmer, mismatches)
        #adjust output from getKmer
        result = data.get("kmers",{})
        if not isinstance(kmer,str):
            for k in kmer:
                if not k in result.keys():
                    result[k] = 0
        else:
            if not kmer in result.keys():
                result[kmer] = 0
        return result


    def getKmerSequence(self, datasetUid, sequence, mismatches=None):
        """
            get k-mers in sequence
        
            :param str datasetUid: dataset uid
            :param Union[str,list] sequence: one or multiple sequences
            :param Union[int,None] mismatches: optional, number of mismatches
            :return: k-mer frequencies
            :rtype: dict
        """
        if (mismatches is None) or (not isinstance(mismatches,int)) or (mismatches<0):
            mismatches = 0
        if not isinstance(sequence,str):
            fullRequest = "{}kmer/{}/sequence".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"sequences": sequence, "mismatches": mismatches}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}kmer/{}/sequence".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"sequence": sequence, "mismatches": mismatches}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getKmerSequenceFrequency(self, datasetUid, sequence, mismatches=None):
        """
            get frequencies for k-mers in sequence (based on getKmerSequence)
        
            :param str datasetUid: dataset uid
            :param Union[str,list] sequence: one or multiple sequences
            :param Union[int,None] mismatches: optional, number of mismatches
            :return: k-mer frequencies
            :rtype: dict
        """
        data = self.getKmerSequence(datasetUid, sequence, mismatches)
        #adjust output from getKmerSequence
        result = data.get("kmers",{})
        return result
        
    def getKmerPath(self, datasetUid, kmerFrom, kmerTo, **kwargs):
        """
            get path
        
            :param str datasetUid: dataset uid
            :param str kmerFrom: starting k-mer
            :param str kmerTo: ending k-mer
            :param Union[int,None] minimumFrequency: optional, minimum frequency
            :param Union[int,None] distance: optional, maximum distance
            :return: path
            :rtype: dict
        """
        parameters = {
            "minimumFrequency": None,
            "distance": None
        }
        fullRequest = "{}kmer/{}/{}/path/{}".format(self._baseUrl, datasetUid, kmerFrom, kmerTo)
        fullRequest = self._addParametersToRequest(fullRequest, parameters, kwargs)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            result = response.json()
            return result
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getKmerSplit(self, datasetUid, kmer, **kwargs):
        """
            get nearest splitting k-mer
        
            :param str datasetUid: dataset uid
            :param str kmer: k-mer
            :param int minimumFrequency: optional, minimum frequency
            :param int distance: optional, maximum distance
            :return: splitting k-mer
            :rtype: dict
        """
        parameters = {
            "minimumFrequency": None,
            "distance": None
        }
        fullRequest = "{}kmer/{}/{}/split".format(self._baseUrl, datasetUid, kmer)
        fullRequest = self._addParametersToRequest(fullRequest, parameters, kwargs)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            print(response,fullRequest)
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getSplitInfo(self, datasetUid):
        """
            get information splitting k-mer database
        
            :param str datasetUid: dataset uid
            :return: information splitting k-mer database
            :rtype: dict
        """
        fullRequest = "{}split/{}/info".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getSplitBase(self, datasetUid:str, base):
        """
            get splitting k-mer base(s)
        
            :param str datasetUid: dataset uid
            :param Union[str,list] base: splitting k-mer base(s)
            :return: splitting k-mer base(s)
            :rtype: Union[dict,list]
        """ 
        if not isinstance(base,str):
            fullRequest = "{}split/{}/base".format(self._baseUrl, datasetUid)
            response = requests.post(fullRequest, json={"bases": sorted(base)}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}split/{}/base/{}".format(self._baseUrl, datasetUid, base)
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getSplitKmer(self, datasetUid:str, kmer):
        """
            get splitting k-mer(s)
        
            :param str datasetUid: dataset uid
            :param Union[str,list] kmer: splitting k-mer(s)
            :return: splitting k-mer(s)
            :rtype: Union[dict,list]
        """ 
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer".format(self._baseUrl, datasetUid)
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}split/{}/kmer/{}".format(self._baseUrl, datasetUid, kmer)
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
    
    def getSplitSequence(self, datasetUid, sequence):
        """
            get splitting k-mers in a sequence
        
            :param str datasetUid: dataset uid
            :param str sequence: sequence
            :return: splitting k-mers
            :rtype: list
        """ 
        fullRequest = "{}split/{}/kmer/sequence".format(self._baseUrl, datasetUid)        
        response = requests.post(fullRequest, json={"sequence": sequence}, 
                                 auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getPositionedSplitSequence(self, datasetUid, sequence):
        """
            get positioned and orientated splitting k-mers in a sequence (based on output getSplitSequence)
        
            :param str datasetUid: dataset uid
            :param str sequence: sequence
            :return: splitting k-mers
            :rtype: dict
        """ 
        splitList = self.getSplitSequence(datasetUid, sequence)
        #adjust output from getSplitSequence
        result = []
        for entry in splitList:
            ckmer = entry["ckmer"]
            reverseCkmer = haplotyping.General.reverse_complement(ckmer)
            result.extend([{"position": m.start(), "orientation": "forward", "data": entry} 
                  for m in re.finditer(entry["ckmer"], sequence)])
            result.extend([{"position": m.start(), "orientation": "forward", "data": entry} 
                  for m in re.finditer(reverseCkmer, sequence)])
        return result
        
    def getSplitDirect(self, datasetUid, kmer):
        """
            get direct neighbours
        
            :param str datasetUid: dataset uid
            :param Union[str,list] kmer: splitting k-mer(s)
            :return: direct neighbours
            :rtype: Union[dict,list]
        """ 
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer/direct".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)},
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}split/{}/kmer/direct/{}".format(self._baseUrl, datasetUid, kmer)        
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None    
            
    def getSplitReads(self, datasetUid, kmer):
        """
            get read information
        
            :param str datasetUid: dataset uid
            :param Union[str,list] kmer: splitting k-mer(s)
            :return: read information
            :rtype: list
        """        
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer/read".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)},
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}split/{}/kmer/read/{}".format(self._baseUrl, datasetUid, kmer)        
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
            
    def getSplitPaired(self, datasetUid, kmer):
        """
            get paired splitting k-mers
        
            :param str datasetUid: dataset uid
            :param Union[str,list] kmer: splitting k-mer(s)
            :return: paired splitting k-mers
            :rtype: Union[list,dict]
        """  
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer/paired".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)},
                                     auth=self._apiAuth, headers=self._apiHeaders)
        else:
            fullRequest = "{}split/{}/kmer/paired/{}".format(self._baseUrl, datasetUid, kmer)        
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerData(self, datasetUid):
        """
            get marker data
        
            :param str datasetUid: dataset uid
            :return: marker data
            :rtype: dict
        """
        fullRequest = "{}marker/{}/data".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerInfo(self, datasetUid):
        """
            get marker info
        
            :param str datasetUid: dataset uid
            :return: marker info
            :rtype: dict
        """
        fullRequest = "{}marker/{}/info".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerMapping(self, datasetUid):
        """
            get marker mapping
        
            :param str datasetUid: dataset uid
            :return: marker mapping
            :rtype: dict
        """
        fullRequest = "{}marker/{}/mapping".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getPedigree(self, varietyUids, ancestorSteps: int = 0, offspringSteps: int = 0, 
                    ancestorOffspringSteps: int = 0, offspringAncestorSteps: int = 0):
        """
            get pedigree
        
            :param Union[str,list] uids: variety uid(s)
            :param int ancestorSteps: optional, number of ancestor steps
            :param int offspringSteps: optional, number of offspring steps
            :param int ancestorOffspringSteps: optional, number of offspring steps for ancestors
            :param int offspringAncestorSteps: optional, number of ancestor steps for offspring
            :return: varieties
            :rtype: dict
        """
        try:
            if isinstance(varietyUids,str):
                varietyUids = [varietyUids]
            def flattenParents(items):
                parents = set()
                for item in items:
                    for entry in item.get("parents",[]):
                        if "uid" in entry:
                            parents.add(entry["uid"])
                        else:
                            parents.update(flattenParents([entry]))
                return parents
            def flattenOffspring(items):
                offspring = set()
                for item in items:
                    for entry in item.get("offspring",[]):
                        if "uid" in entry:
                            offspring.add(entry["uid"])
                return offspring
            def getAncestors(varietyUids, steps):
                if steps==0 or len(varietyUids)==0:
                    return {}
                else:
                    subResult = self.getVarietyById(varietyUids)
                    subResult.update(getAncestors(flattenParents(subResult.values()),steps-1))
                    return subResult        
            def getOffspring(varietyUids, steps):
                if steps==0 or len(varietyUids)==0:
                    return {}
                else:
                    subResult = self.getVarietyById(varietyUids)
                    subResult.update(getOffspring(flattenOffspring(subResult.values()),steps-1))
                    return subResult        
            result = self.getVarietyById(varietyUids)
            ancestorResult = getAncestors(flattenParents(result.values()), ancestorSteps)
            ancestorOffspringResult = getOffspring(flattenOffspring(ancestorResult.values()), ancestorOffspringSteps)
            offspringResult = getOffspring(flattenOffspring(result.values()), offspringSteps)
            offspringAncestorResult = getAncestors(flattenParents(offspringResult.values()), offspringAncestorSteps)
            result.update(ancestorResult)
            result.update(ancestorOffspringResult)
            result.update(offspringResult)
            result.update(offspringAncestorResult)
            return result
        except:
            self._api_logger.error("getting pedigree didn't succeed")
            
    def visualizeVarieties(self, varieties: dict, selected: list = None, coloring: dict = None, name=None):
        """
            get visualization for set of varieties (pedigree)
        
            :param dict varieties: varieties
            :param list selected: optional, list of varieties to be marked as selected
            :param dict coloring: optional, for each color a list of varieties to be marked with this color
            :return: visualization
            :rtype: graphviz object
        """
        varietyColoring = {}
        selectedColor = "lightblue"
        defaultColor = "azure"
        try:
            if not coloring==None:
                duplicateError = False
                for color in coloring:
                    for uid in coloring[color]:
                        if not isinstance(uid,str):
                            continue
                        if uid in varietyColoring:
                            duplicateError = True
                        else:
                            varietyColoring[uid] = color
                if duplicateError:
                    self._api_logger.error("duplicate definition in coloring")
            if not selected==None:
                for uid in selected:
                    if not isinstance(uid,str):
                        continue
                    if not uid in varietyColoring:
                        varietyColoring[uid] = selectedColor
        except:
            self._api_logger.error("incorrect definition selected or coloring")
        pedigree_title = "Pedigree"
        if not name==None:
            pedigree_title = name
        pedigree_title = "<<font point-size=\"14\">{}</font>>".format(html.escape(pedigree_title))
        
        g = Digraph("Pedigree")
        g.attr(label=pedigree_title, labelloc="t", nodesep="0", ranksep="1")
        for uid in varieties.keys():
            variety = varieties[uid]
            node_title="<"
            node_title+="<font point-size=\"8\">{}</font><br/>".format(html.escape(variety["name"]))
            if "origin" in variety:
                node_title+="<font point-size=\"6\">{}</font>".format(html.escape(variety["origin"]["uid"]))
                if "year" in variety:
                    node_title+="<font point-size=\"6\">&nbsp;</font>"
            if "year" in variety:
                node_title+="<font point-size=\"6\">{}</font>".format(html.escape(variety["year"]["description"]))
            node_title+=">"
            node_fillcolor = varietyColoring.get(uid,defaultColor)
            g.node("variety_{}".format(uid),label=node_title, color="black", style="filled", fillcolor=node_fillcolor)
        def flattenParents(entry):
            uids = []
            if "parents" in entry:
                for parent in entry["parents"]:
                    if "uid" in parent:
                        uids.append(parent["uid"])
                    else:
                        uids.extend(flattenParents(parent))
            return uids
        def linkParents(g,node_key,entry,counter=0):   
            if len(set(flattenParents(entry)).intersection(varieties.keys()))>0:
                if "parents" in entry:
                    for parent in entry["parents"]:
                        if "uid" in parent:
                            if parent["uid"] in varieties.keys():
                                parent_key = "variety_{}".format(parent["uid"])
                                g.edge(parent_key,node_key,label="",style="solid", color="black",dir="forward")
                        else:
                            parent_key = "parent_{}_{}".format(uid,counter)
                            g.node(parent_key,label="", shape="diamond", color="black", style="filled", fillcolor="azure")
                            g.edge(parent_key,node_key,label="",style="solid", color="black",dir="forward")
                            linkParents(g,parent_key,parent,counter+1)
        for uid in varieties.keys():
            linkParents(g,"variety_{}".format(uid),varieties[uid])
        return g
