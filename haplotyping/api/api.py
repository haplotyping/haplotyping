import logging,requests, html
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
        fullRequest = "{}{}{}".format(self._baseUrl,request,requests.utils.quote(uid))
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
        
    
    def getCountries(self): 
        """
            get a dict with all countries
        
            :return: all countries
            :rtype: dict
        """
        return self._getList("country/","uid")
        
                
    def getCollections(self):
        """
            get a dict with all collections
        
            :return: all collections
            :rtype: dict
        """
        return self._getList("collection/","uid")
    
    def getDatasets(self, **kwargs):
        """
            get a dict with (all) datasets

            :param str collection: optional, comma separated list of collection uids
            :param str dataType: optional, one of "marker", "kmer", "split"
            :param bool hasVariety: optional, has a linked variety
            :return: (all) datasets
            :rtype: dict
        """
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
        """
            get dataset
        
            :param str/list uid: dataset uid(s)
            :return: dataset
            :rtype: dict
        """
        if not isinstance(uid,str):
            return self._getByIds("dataset/",uid)
        else:
            return self._getById("dataset/",uid)
    
    def getVarieties(self, **kwargs):
        """
            get a dict with (all) varieties

            :param str name: optional, name or synonym
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
        return self._getList(request,"uid")
    
    def getVarietyById(self, uid):
        """
            get variety
        
            :param str/list uid: variety uid(s)
            :return: variety
            :rtype: dict
        """
        if not isinstance(uid,str):
            return self._getByIds("variety/",uid)
        else:
            return self._getById("variety/",uid)
    
    def getKmerFrequency(self, uid, kmer, mismatches=0):
        """
            get k-mer frequency
        
            :param str uid: dataset uid
            :param str/list kmer: k-mer(s)
            :param int mismatches: optional, number of mismatches
            :return: k-mer frequencies
            :rtype: dict
        """
        if not isinstance(kmer,str):
            fullRequest = "{}kmer/{}".format(self._baseUrl, uid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer),"mismatches": mismatches}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                result = data.get("kmers",{})
                for k in kmer:
                    if not k in result.keys():
                        result[k] = 0
                return result
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None
        else:
            fullRequest = "{}kmer/{}/{}{}".format(self._baseUrl, uid, kmer, 
                                              "?mismatches={}".format(int(mismatches)) if mismatches>0 else "")
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                result = data.get("kmers",{})
                if not kmer in result.keys():
                    result[kmer] = 0
                return result
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None
    
    def getKmerSequence(self, uid, sequence, mismatches=0):
        """
            get frequencies for k-mers in sequence
        
            :param str uid: dataset uid
            :param str sequence: sequence
            :param int mismatches: optional, number of mismatches
            :return: k-mer frequencies
            :rtype: dict
        """
        fullRequest = "{}kmer/{}/sequence".format(self._baseUrl, uid)        
        response = requests.post(fullRequest, json={"sequence": sequence,"mismatches": mismatches}, 
                                 auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            result = data.get("kmers",{})
            return result
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getSplitInfo(self, datasetUid):
        fullRequest = "{}split/{}/info".format(self._baseUrl, datasetUid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getSplitKmer(self, datasetUid, kmer):
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer".format(self._baseUrl, datasetUid)
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)}, 
                                     auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                return data
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None
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
        fullRequest = "{}split/{}/kmer/sequence".format(self._baseUrl, datasetUid)        
        response = requests.post(fullRequest, json={"sequence": sequence}, 
                                 auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getSplitDirect(self, datasetUid, kmer):
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer/direct".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)},
                                     auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                return data
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None
        else:
            fullRequest = "{}split/{}/kmer/direct/{}".format(self._baseUrl, datasetUid, kmer)        
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                return data
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None        
        
    def getSplitConnected(self, datasetUid, kmer):
        if not isinstance(kmer,str):
            fullRequest = "{}split/{}/kmer/connected".format(self._baseUrl, datasetUid)        
            response = requests.post(fullRequest, json={"kmers": sorted(kmer)},
                                     auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                return data
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None
        else:
            fullRequest = "{}split/{}/kmer/connected/{}".format(self._baseUrl, datasetUid, kmer)        
            response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                return data
            else:
                self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
                return None        
        
    def getMarkerData(self, uid):
        """
            get marker data
        
            :param str uid: dataset uid
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
        
    def getMarkerInfo(self, uid):
        """
            get marker info
        
            :param str uid: dataset uid
            :return: marker info
            :rtype: dict
        """
        fullRequest = "{}marker/{}/info".format(self._baseUrl, uid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getMarkerMapping(self, uid):
        """
            get marker mapping
        
            :param str uid: dataset uid
            :return: marker mapping
            :rtype: dict
        """
        fullRequest = "{}marker/{}/mapping".format(self._baseUrl, uid)
        response = requests.get(fullRequest, auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            self._api_logger.error("request to {} didn't succeed".format(self._baseUrl))
            return None
        
    def getPedigree(self, uids, ancestorSteps: int = 0, offspringSteps: int = 0, 
                    ancestorOffspringSteps: int = 0, offspringAncestorSteps: int = 0):
        """
            get pedigree
        
            :param str/list uids: dataset uid(s)
            :param int ancestorSteps: optional, number of ancestor steps
            :param int offspringSteps: optional, number of offspring steps
            :param int ancestorOffspringSteps: optional, number of offspring steps for ancestors
            :param int offspringAncestorSteps: optional, number of ancestor steps for offspring
            :return: varieties
            :rtype: dict
        """
        try:
            if isinstance(uids,str):
                uids = [uids]
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
            def getAncestors(uids, steps):
                if steps==0 or len(uids)==0:
                    return {}
                else:
                    subResult = self.getVarietyById(uids)
                    subResult.update(getAncestors(flattenParents(subResult.values()),steps-1))
                    return subResult        
            def getOffspring(uids, steps):
                if steps==0 or len(uids)==0:
                    return {}
                else:
                    subResult = self.getVarietyById(uids)
                    subResult.update(getOffspring(flattenOffspring(subResult.values()),steps-1))
                    return subResult        
            result = self.getVarietyById(uids)
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
            
    def visualizeVarieties(self, varieties: dict, selected: list = None, coloring: dict = None):
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
        g = Digraph("Pedigree")
        g.attr(label="Pedigree", labelloc="t", nodesep="0", ranksep="1")
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
