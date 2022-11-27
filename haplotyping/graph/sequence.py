import logging,requests,re
from haplotyping.graph.api import APIGraph
from haplotyping.general import General
import haplotyping

#TODO: minimum frequency and maximum distance centrally configured

class SequenceGraph(APIGraph):
    """constructing the De Bruijn graph"""
    
    def __init__(self, url: str, uid: str, sequence: str, expandSteps: int = 10, 
                 ignoreMissingPath: bool = True, username: str = None, password: str = None):
        """
        construct the De Bruijn graph from a reference sequence using the dataset 
        defined by uid via the API at provided url
        """
        #call parent constructor
        super(haplotyping.graph.sequence.SequenceGraph, self).__init__(url, uid, username, password)
        #logger
        self._logger = logging.getLogger(__name__)
        #store call parameters
        self._sequence = str(sequence).strip()
        self._expandSteps = expandSteps
        self._ignoreMissingPath = ignoreMissingPath
        self._logger.info("construct De Bruijn graph for sequence of length "+str(len(self._sequence)))
        #construct by expanding
        self._sequenceConstructGraph()
        #try to fix missing connections in the graph
        self._fixMissingConnections()
        #try to glue missing connections
        self._glueMissingConnections()
        #fix start and end based on connected candidates
        self._expandConnectedStartEndCandidates()
        
    def __repr__(self):
        text = super(haplotyping.graph.sequence.SequenceGraph, self).__repr__()
        text = text + " for sequence of length %d" % (len(self._sequence))
        return text
    
    def getOrientatedCkmers(self, minimumPosition : int = None, maximumPosition : int = None):
        result = super(haplotyping.graph.sequence.SequenceGraph, self).getOrientatedCkmers()
        if minimumPosition==None and maximumPosition==None:
            return result
        else:
            filteredResult = set()
            for orientatedCkmer in result:
                entry = self._orientatedCkmers[orientatedCkmer]
                if entry._order==None:
                    continue
                else:
                    if not minimumPosition==None and entry._order<minimumPosition:
                        continue
                    if not maximumPosition==None and entry._order>maximumPosition:
                        continue
                filteredResult.add(orientatedCkmer)                                            
            return filteredResult

    def _sequenceConstructGraph(self):
        """
        main function to construct the De Bruijn Graph
        """
        self._logger.debug("start constructing the De Bruijn Graph")  
        #initialise
        self._connectedCkmers = {}
        #get main k-mers from sequence
        response = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/sequence", json = {"sequence": self._sequence},
                                 auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            if len(data)==0:
                self._logger.error("no splitting k-mers found in sequence")
            else:
                #pre-compute positions
                forward_positions = {}
                backward_positions = {}
                #initialise minimum and maximum
                positionMin = None
                positionMax = None
                startCandidate = None
                endCandidate = None
                dataIndex = {item["ckmer"]: item for item in data}
                #create initial k-mers
                for i in range(len(self._sequence)-(self._k-1)):
                    kmer = self._sequence[i:i+self._k]
                    ckmer = haplotyping.General.canonical(kmer)
                    if kmer==ckmer:
                        if ckmer in forward_positions.keys():
                            forward_positions[ckmer].append(i)
                        else:
                            forward_positions[ckmer]=[i]
                        direction="forward"
                    else:
                        if ckmer in backward_positions.keys():
                            backward_positions[ckmer].append(i)
                        else:
                            backward_positions[ckmer]=[i]
                        direction="backward"
                    if ckmer in dataIndex.keys():
                        orientatedCkmerKey = (ckmer, direction)
                        if not orientatedCkmerKey in self._orientatedCkmers.keys():
                            newEntry = self._createOrientatedCkmer(dataIndex[ckmer], direction)
                            newEntry._setPosition(i)
                        else:
                            newEntry = self._orientatedCkmers[orientatedCkmerKey]
                        newEntry._setLevel(0)
                        newEntry._setCandidate()
                        self._selected.add(orientatedCkmerKey)
                        if positionMin==None or i<positionMin:
                            startCandidate = orientatedCkmerKey
                            positionMin = i
                        if positionMax==None or i>positionMax:
                            endCandidate = orientatedCkmerKey
                            positionMax = i
                self._logger.debug("found {} splitting k-mers in sequence".format(len(self._orientatedCkmers)))
                
                #compute (and possibly add) best start candidate(s)
                searchFinished = False
                for i in range(positionMin):
                    startKmer = self._sequence[i:self._k+i]
                    response = requests.get(self._baseUrl+"kmer/"+self._uid+"/"+str(startKmer)+"?mismatches=2",
                                            auth=self._apiAuth, headers=self._apiHeaders)
                    if response.ok:
                        data = response.json()
                        #found a hit
                        if len(data["kmers"])>0:
                            #expand hits to the first splitting k-mer
                            for kmer in data["kmers"].keys():
                                splitData = self._findSplit(kmer,"right")
                                if len(splitData["orientatedCkmers"])>0:
                                    searchFinished = True
                                    self._logger.debug("process {} k-mers as potential start candidate(s)".format(
                                        len(splitData["orientatedCkmers"]))) 
                                    #alternative route
                                    if (i+splitData["distance"])>positionMin:
                                        #mark found splitting k-mers as connected on incoming side
                                        for orientatedCkmerKey in splitData["orientatedCkmers"]:
                                            if orientatedCkmerKey in self._orientatedCkmers.keys():
                                                self._orientatedCkmers[orientatedCkmerKey]._setConnectedViaCandidate(
                                                    "incoming") 
                                            else:
                                                response = requests.get(
                                                    self._baseUrl+"split/"+self._uid+"/kmer/"+str(orientatedCkmerKey[0]),
                                                    auth=self._apiAuth, headers=self._apiHeaders)
                                                if response.ok:
                                                    data = response.json()
                                                    if data:
                                                        orientatedCkmer = self._createOrientatedCkmer(data, orientatedCkmerKey[1])
                                                        orientatedCkmer._setConnectedViaCandidate("incoming")
                                    else:
                                        #mark found splitting k-mers as candidates
                                        for orientatedCkmerKey in splitData["orientatedCkmers"]:
                                            if orientatedCkmerKey in self._orientatedCkmers.keys():
                                                self._orientatedCkmers[orientatedCkmerKey]._setCandidate()
                                                self._orientatedCkmers[orientatedCkmerKey]._setLevel(1) 
                                                self._orientatedCkmers[orientatedCkmerKey]._setPosition(
                                                   i+splitData["distance"])
                                            else:
                                                response = requests.get(
                                                    self._baseUrl+"split/"+self._uid+"/kmer/"+str(orientatedCkmerKey[0]),
                                                    auth=self._apiAuth, headers=self._apiHeaders)
                                                if response.ok:
                                                    data = response.json()
                                                    if data:
                                                        orientatedCkmer = self._createOrientatedCkmer(data, orientatedCkmerKey[1])
                                                        orientatedCkmer._setCandidate()
                                                        orientatedCkmer._setLevel(1) 
                                                        orientatedCkmer._setPosition(i+splitData["distance"])
                                                else:
                                                    self._logger.debug("request for {} didn't succeed".format(
                                                        orientatedCkmerKey[0]))
                                            self._setStart(orientatedCkmerKey)                            
                                    if (i+splitData["distance"])>=positionMin:
                                        #set first splitting k-mer from sequence as start candidate  
                                        self._setStart(startCandidate)                            
                        if searchFinished:
                            break
                    else:
                        self._logger.error("request for {} didn't succeed".format(startKmer))                     
                #at least one start candidate
                if len(self._start)==0:
                    #set first splitting k-mer from sequence as start candidate  
                    self._setStart(startCandidate)
                self._logger.debug("marked {} orientated k-mers as start candidate(s)".format(len(self._start)))         
                        
                #compute (and possibly add) best end candidate(s)
                searchFinished = False                
                for i in range(len(self._sequence)-self._k,positionMax-1,-1):
                    endKmer = self._sequence[i:self._k+i]
                    response = requests.get(self._baseUrl+"kmer/"+self._uid+"/"+str(endKmer)+"?mismatches=2",
                                            auth=self._apiAuth, headers=self._apiHeaders)
                    if response.ok:
                        data = response.json()
                        #found a hit
                        if len(data["kmers"])>0:
                            #expand hits to the first splitting k-mer
                            for kmer in data["kmers"].keys():
                                splitData = self._findSplit(kmer,"left")
                                if len(splitData["orientatedCkmers"])>0:
                                    searchFinished = True
                                    self._logger.debug("process {} k-mers as potential end candidate(s)".format(
                                        len(splitData["orientatedCkmers"]))) 
                                    #alternative route
                                    if (i-splitData["distance"])<positionMax:
                                        #mark found splitting k-mers as connected on incoming side
                                        for orientatedCkmerKey in splitData["orientatedCkmers"]:
                                            if orientatedCkmerKey in self._orientatedCkmers.keys():
                                                self._orientatedCkmers[orientatedCkmerKey]._setConnectedViaCandidate(
                                                    "outgoing") 
                                            else:
                                                response = requests.get(
                                                    self._baseUrl+"split/"+self._uid+"/kmer"+str(orientatedCkmerKey[0]),
                                                    auth=self._apiAuth, headers=self._apiHeaders)
                                                if response.ok:
                                                    data = response.json()
                                                    if data:
                                                        orientatedCkmer = self._createOrientatedCkmer(data, orientatedCkmerKey[1])
                                                        orientatedCkmer._setConnectedViaCandidate("outgoing")
                                    else:
                                        #mark found splitting k-mers as candidates
                                        for orientatedCkmerKey in splitData["orientatedCkmers"]:
                                            if orientatedCkmerKey in self._orientatedCkmers.keys():
                                                self._orientatedCkmers[orientatedCkmerKey]._setCandidate()
                                                self._orientatedCkmers[orientatedCkmerKey]._setLevel(1) 
                                                self._orientatedCkmers[orientatedCkmerKey]._setPosition(
                                                   i-splitData["distance"])
                                            else:
                                                response = requests.get(
                                                    self._baseUrl+"split/"+self._uid+"/kmer/"+str(orientatedCkmerKey[0]),
                                                    auth=self._apiAuth, headers=self._apiHeaders)
                                                if response.ok:
                                                    data = response.json()
                                                    if data:
                                                        orientatedCkmer = self._createOrientatedCkmer(data, orientatedCkmerKey[1])
                                                        orientatedCkmer._setCandidate()
                                                        orientatedCkmer._setLevel(1) 
                                                        orientatedCkmer._setPosition(i-splitData["distance"])
                                                else:
                                                    self._logger.debug("request for {} didn't succeed".format(
                                                        orientatedCkmerKey[0]))
                                            self._setEnd(orientatedCkmerKey)                            
                                    if (i-splitData["distance"])<=positionMax:
                                        #set last splitting k-mer from sequence as end candidate  
                                        self._setEnd(endCandidate)                                 
                        if searchFinished:
                            break
                    else:
                        self._logger.error("request for {} didn't succeed".format(endKmer))                     
                #at least one end candidate
                if len(self._end)==0:
                    #set last splitting k-mer from sequence as end candidate  
                    self._setEnd(endCandidate)
                self._logger.debug("marked {} orientated k-mers as end candidate(s)".format(len(self._end))) 
                
                #expand k-mer set from the connected set
                candidate_kmer_keys = [k for k in self._orientatedCkmers.keys() 
                                       if self._orientatedCkmers[k].candidate()]
                candidate_kmers = set([k[0] for k in candidate_kmer_keys])
                #expand k-mers
                for level in range(self._expandSteps):                    
                    self._expandOrientatedCkmers(level)                    
        else:
            self._logger.error("request to {} didn't succeed".format(self._baseUrl)) 

    def _expandOrientatedCkmers(self, level: int):
        """
        expand k-mers from defined level by finding direct connections
        used in _sequenceConstructGraph()
        """
        assert level>=0
        #define list of k-mers to check
        kmer_keys = [k for k in self._orientatedCkmers.keys() 
                     if (not self._orientatedCkmers[k]._level==None and self._orientatedCkmers[k]._level<=level) and 
                     not self._orientatedCkmers[k]._expanded]
        #define orientated right split bases containing a candidate
        oRightSplitBases_candidate = set([b for b in self._orientatedBases.keys() if self._orientatedBases[b].candidate()])
        #define orientated canonical k-mers with candidate base
        ocKmers_candidateBases = [self._orientatedBases[b]._orientatedCkmers for b in oRightSplitBases_candidate]
        ocKmers_candidateBases = set([item for sublist in ocKmers_candidateBases for item in sublist])
        #extend kmer-keys
        kmer_keys.extend([c for c in ocKmers_candidateBases if not self._orientatedCkmers[c]._expanded])
        #process
        kmers = set([k[0] for k in kmer_keys])
        initial_candidates = sum(self._orientatedCkmers[ckmer_key].candidate() for ckmer_key in self._orientatedCkmers)
        self._logger.debug("from %d k-mers with %d candidates, expand %d at level %d" % 
                           (len(self._orientatedCkmers), initial_candidates, len(kmers), level)) 
        #get direct connections
        response = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/direct", 
                                 json = {"kmers": sorted(kmers)}, auth=self._apiAuth, headers=self._apiHeaders)
        counter_direct_connections = 0
        counter_skipped_connections = 0
        set_direct_kmers = set()
        set_new_kmers = set()
        if response.ok:
            data = response.json()
            for item in data:
                ckmer1 = item["ckmer"]
                for orientation1 in ["forward","backward"]:
                    ckmerKey1 = (ckmer1,orientation1)
                    if ckmerKey1 in self._orientatedCkmers.keys():
                        for direction1 in item["direct"]:
                            assert direction1 in ["left","right"]
                            for directionItem in item["direct"][direction1]:                                
                                ckmer2 = directionItem["ckmer"]
                                #administration
                                counter_direct_connections+=1
                                set_direct_kmers.add(ckmer2)
                                direction2 = directionItem["connection"]["direction"]
                                assert direction2 in ["left","right"]
                                #inward/outward
                                if ((orientation1=="forward" and direction1=="right") or
                                    (orientation1=="backward" and direction1=="left")):
                                    #outward
                                    orientation2 = "backward" if direction2=="right" else "forward"
                                    ckmerKey2 = (ckmer2,orientation2) 
                                    if not ckmerKey2 in self._orientatedCkmers.keys():
                                        self._createOrientatedCkmer(directionItem, orientation2)
                                        set_new_kmers.add(ckmer2)
                                    self._orientatedCkmers[ckmerKey1]._setOutgoing(ckmerKey2, 
                                                                         directionItem["connection"]["distance"], 
                                                                         directionItem["connection"]["number"],
                                                                         directionItem["connection"]["problem"]>0)
                                    self._orientatedCkmers[ckmerKey2]._setIncoming(ckmerKey1, 
                                                                         directionItem["connection"]["distance"], 
                                                                         directionItem["connection"]["number"],
                                                                         directionItem["connection"]["problem"]>0)            
                                else:
                                    #inward
                                    orientation2 = "forward" if direction2=="right" else "backward"
                                    ckmerKey2 = (ckmer2,orientation2) 
                                    if not ckmerKey2 in self._orientatedCkmers.keys():
                                        self._createOrientatedCkmer(directionItem, orientation2)
                                        set_new_kmers.add(ckmer2)
                                    self._orientatedCkmers[ckmerKey1]._setIncoming(ckmerKey2, 
                                                                         directionItem["connection"]["distance"], 
                                                                         directionItem["connection"]["number"],
                                                                         directionItem["connection"]["problem"]>0)  
                                    self._orientatedCkmers[ckmerKey2]._setOutgoing(ckmerKey1, 
                                                                         directionItem["connection"]["distance"], 
                                                                         directionItem["connection"]["number"],
                                                                         directionItem["connection"]["problem"]>0)  
            final_candidates = sum(self._orientatedCkmers[ckmer].candidate() for ckmer in self._orientatedCkmers)
            self._logger.debug("skipped %d and found %d connections to %d k-mers: %d new and %d new candidates" % 
                               (counter_skipped_connections, counter_direct_connections,len(set_direct_kmers),
                                len(set_new_kmers),(final_candidates-initial_candidates)))
            #set all k-mers to expanded
            for k in kmer_keys:
                self._orientatedCkmers[k]._expanded = True                                    
        else:
            self._logger.error("request to {} didn't succeed".format(self._baseUrl))         
            
    def _createOrientatedCkmer(self, item: dict, orientation: str):         
        """
        create k-mer, used in _sequenceConstructGraph
        """
        assert "ckmer" in item.keys() and "number" in item.keys() and "split" in item.keys()
        assert orientation in ["forward","backward"]
        ckmer = item["ckmer"]
        orientatedCkmerKey = (ckmer, orientation)
        if not orientatedCkmerKey in self._orientatedCkmers.keys():
            orientatedCkmerEntry = self.SequenceCkmer(self, ckmer, orientation, item["number"], item["split"])                 
        else:
            self._logger.warning("orientated k-mer already exists")
        return self._orientatedCkmers[orientatedCkmerKey]
    
    def _expandConnectedStartEndCandidates(self):
        connectedSets = self.getConnectedCandidates(True)
        expandedStartCandidates = 0
        expandedEndCandidates = 0
        #get start and end entries
        startEntries = set()
        endEntries = set()
        startAlternatives = set()
        endAlternatives = set()
        removedConnectedSets = 0
        for connectedSet in connectedSets:            
            if len(connectedSet["connected"])==1:
                for orientatedCkmerKey in connectedSet["connected"]:
                    self._orientatedCkmers[orientatedCkmerKey]._unsetCandidate()
                for orientatedCkmerKey in connectedSet["start"]:
                    self._unsetStart(orientatedCkmerKey)
                for orientatedCkmerKey in connectedSet["end"]:
                    self._unsetEnd(orientatedCkmerKey)
                removedConnectedSets+=1
            else:
                connectedList = sorted(connectedSet["connected"])
                for orientatedCkmerKey in connectedSet["start"]:
                    self._setStart(orientatedCkmerKey)
                    orientatedCkmer = self._orientatedCkmers[orientatedCkmerKey]
                    assert orientatedCkmer.candidate()
                    if orientatedCkmerKey[1]=="backward":
                        startEntries.add(haplotyping.General.reverse_complement(orientatedCkmerKey[0]))
                    else:
                        startEntries.add(orientatedCkmerKey[0])                
                for orientatedCkmerKey in connectedSet["end"]:
                    self._setEnd(orientatedCkmerKey)
                    orientatedCkmer = self._orientatedCkmers[orientatedCkmerKey]
                    assert orientatedCkmer.candidate()
                    if orientatedCkmerKey[1]=="backward":
                        endEntries.add(haplotyping.General.reverse_complement(orientatedCkmerKey[0]))
                    else:
                        endEntries.add(orientatedCkmerKey[0])
        #find alternatives
        response = requests.post(self._baseUrl+"kmer/"+self._uid, 
                                 json = {"kmers": sorted(startEntries), "mismatches": 2}, 
                                 auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            startAlternatives = set(data["kmers"].keys()).difference(startEntries)
        response = requests.post(self._baseUrl+"kmer/"+self._uid, 
                                 json = {"kmers": sorted(endEntries), "mismatches": 2},
                                 auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            endAlternatives = set(data["kmers"].keys()).difference(endEntries)
        #find splits
        shortestDistances = self.getShortestDistances()
        for kmer in startAlternatives:
            foundSplits = self._findSplit(kmer,"right")
            for orientatedCkmerKey in foundSplits["orientatedCkmers"]:
                if orientatedCkmerKey in self._orientatedCkmers and orientatedCkmerKey in shortestDistances.keys():
                    connected = set([k for k in shortestDistances[orientatedCkmerKey].keys() 
                                     if shortestDistances[orientatedCkmerKey][k]!=0])
                    #only if connected to ending k-mers
                    if ((len(connected.intersection(endAlternatives))>0 or 
                        len(connected.intersection(self._end))>0) 
                        and not self._orientatedCkmers[orientatedCkmerKey].candidate()):
                        self._setStart(orientatedCkmerKey)  
                        expandedStartCandidates+=1
        for kmer in endAlternatives:
            foundSplits = self._findSplit(kmer,"left")
            for orientatedCkmerKey in foundSplits["orientatedCkmers"]:
                if orientatedCkmerKey in self._orientatedCkmers.keys() and orientatedCkmerKey in shortestDistances.keys():
                    connected = set([k for k in shortestDistances[orientatedCkmerKey].keys() 
                                     if shortestDistances[orientatedCkmerKey][k]!=0])
                    #only if connected to ending k-mers
                    if ((len(connected.intersection(startAlternatives))>0 or 
                        len(connected.intersection(self._start))>0) 
                        and not self._orientatedCkmers[orientatedCkmerKey].candidate()):
                        self._setEnd(orientatedCkmerKey)
                        expandedEndCandidates+=1
        self._logger.debug("expanded with {} start and {} end k-mers for {} connected sets of candidates".format(
           expandedStartCandidates,expandedEndCandidates,len(connectedSets)-removedConnectedSets))      
        #reset incorrect start and end
        correctedStartCandidates = 0
        correctedEndCandidates = 0
        for orientatedCkmer in self._orientatedCkmers:
            if orientatedCkmer in self._start:
                if len([k for k in self._start if shortestDistances[k][orientatedCkmer]>0])>0:
                    self._unsetStart(orientatedCkmer)
                    correctedStartCandidates+=1
            if orientatedCkmer in self._end:
                if len([k for k in self._end if shortestDistances[orientatedCkmer][k]>0])>0:
                    self._unsetEnd(orientatedCkmer)
                    correctedEndCandidates+=1
        if correctedStartCandidates>0 or correctedEndCandidates>0:
            self._logger.debug("corrected {} incorrect start and {} incorrect end k-mers".format(
               expandedStartCandidates,expandedEndCandidates,len(connectedSets)-removedConnectedSets)) 
                          
                            
    def _findPath(self, orientedCkmerFrom, orientedCkmerTo, distance:int):
        minimumFrequency = 2
        kmerFrom = (orientedCkmerFrom[0] if orientedCkmerFrom[1]=="forward" 
                    else haplotyping.General.reverse_complement(orientedCkmerFrom[0]))
        kmerTo = (orientedCkmerTo[0] if orientedCkmerTo[1]=="forward" 
                    else haplotyping.General.reverse_complement(orientedCkmerTo[0]))
        response = requests.get("{}kmer/{}/{}/path/{}?minimumFrequency={}&distance={}".format(
            self._baseUrl,self._uid,kmerFrom,kmerTo,minimumFrequency,distance),
                                auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            return data
        else:
            return None
        
    def _findSplit(self, kmer:str, direction:str):
        assert direction in ["right","left"]
        #check split
        response = requests.get(self._baseUrl+"split/"+self._uid+"/kmer/"+str(kmer),
                                auth=self._apiAuth, headers=self._apiHeaders)
        if response.ok:
            data = response.json()
            if data and data["ckmer"]:
                newKmer = kmer if direction=="right" else haplotyping.General.reverse_complement(kmer)
                newCkmer = haplotyping.General.canonical(kmer)
                if direction=="right":
                    return {"distance": 0, 
                            "orientatedCkmers":[(newCkmer,"forward" if newKmer==newCkmer else "backward")], 
                            "pathKmers": []}
                else:
                    return {"distance": 0, 
                            "orientatedCkmers":[(newCkmer,"forward" if not newKmer==newCkmer else "backward")], 
                            "pathKmers": []}
        #check neighbours
        maximumDistance = 1000
        minimumFrequency = 2        
        newKmer = kmer if direction=="right" else haplotyping.General.reverse_complement(kmer)    
        response = requests.get("{}kmer/{}/{}/split?minimumFrequency={}&distance={}".format(
            self._baseUrl,self._uid,newKmer,minimumFrequency,maximumDistance),
                                auth=self._apiAuth, headers=self._apiHeaders)
        result = {"distance": None, "orientatedCkmers":[], "pathKmers": [newKmer]}
        if response.ok:
            data = response.json()
            result["distance"] = data["distance"]
            result["pathKmers"] = data["pathKmers"]
            for splittingKmer in data["splittingKmers"]:
                splittingCkmer = haplotyping.General.canonical(splittingKmer)
                if direction=="right":
                    result["orientatedCkmers"].append(
                        (splittingCkmer, "forward" if splittingKmer==splittingCkmer else "backward"))
                else:
                    result["orientatedCkmers"].append(
                        (splittingCkmer, "forward" if not splittingKmer==splittingCkmer else "backward"))     
        return result        
    
    def _fixMissingConnections(self):
        #compute k-mers with only single orientation (not necessary to extend this to double orientation?)
        relevantKmers = set([c for c in self._ckmers if len(self._ckmers[c])==1])
        #compute distances
        shortestDistances = self.getShortestDistances()
        #define orientated canonical candidate k-mers
        ocKmers_candidate = set([c for c in self._orientatedCkmers.keys() if self._orientatedCkmers[c].candidate()])
        #define orientated right split bases containing a candidate
        oRightSplitBases_candidate = set([b for b in self._orientatedBases.keys() if self._orientatedBases[b].candidate()])
        #define orientated canonical k-mers with candidate base
        ocKmers_candidateBases = [self._orientatedBases[b]._orientatedCkmers for b in oRightSplitBases_candidate]
        ocKmers_candidateBases = set([item for sublist in ocKmers_candidateBases for item in sublist])
        #define the set of relevant open orientated k-mers: open is missing connections in at least one direction
        ocKmers_open = set()
        for c in ocKmers_candidateBases:
            c_data = self._orientatedCkmers[c]
            if (len(ocKmers_candidateBases.intersection(c_data._incoming))==0 or
               len(ocKmers_candidateBases.intersection(c_data._outgoing))==0):
                ocKmers_open.add(c)
        self._logger.debug("found %d relevant k-mers with missing connections to other candidates" % (len(ocKmers_open)))
        
        #use knowledge from connections for k-mers with same base    
        numberOfFixedMissingConnections = 0
        for c in ocKmers_candidateBases:
            c_data = self._orientatedCkmers[c]
            assert c_data._orientation in ["forward","backward"]            
            incomingRelevantSplitSide = "right" if c_data._orientation == "forward" else "left"
            outgoingRelevantSplitSide = "right" if c_data._orientation == "backward" else "left"
            #check for missing incoming
            if incomingRelevantSplitSide in c_data._orientatedBases:
                b = c_data._orientatedBases[incomingRelevantSplitSide]
                assert (b in oRightSplitBases_candidate) or not c_data.candidate()
                if b in oRightSplitBases_candidate:
                    b_data = self._orientatedBases[b]
                    assert c in b_data._orientatedCkmers
                    known_incoming_kmers = ocKmers_candidateBases.intersection(c_data._incoming)
                    potential_incoming_kmers = set()
                    for c2 in ocKmers_candidateBases.intersection(b_data._orientatedCkmers):
                        potential_incoming_kmers.update(self._orientatedCkmers[c2]._incoming.keys())
                    potential_incoming_kmers = ocKmers_candidateBases.intersection(potential_incoming_kmers)
                    if len(potential_incoming_kmers.difference(known_incoming_kmers))>0:
                        new_incoming_kmers = potential_incoming_kmers.difference(known_incoming_kmers)
                        distance = self.getShortestDistance(potential_incoming_kmers,b_data._orientatedCkmers)
                        #if there are no connections at all, just add all possible edges
                        #this can however introduce false connections
                        if len(known_incoming_kmers)==0:
                            for c3 in new_incoming_kmers:
                                numberOfFixedMissingConnections+=1
                                self._orientatedCkmers[c]._setIncoming(c3,distance,0,False) 
                                if c[0]==c3[0]:
                                    self._orientatedCkmers[c3]._setIncoming(c,distance,0,False) 
                                else:
                                    self._orientatedCkmers[c3]._setOutgoing(c,distance,0,False) 
                        #otherwise, only check if potential incoming can restore missing connections
                        else:
                            #find connected
                            connectedOrientatedKmers = set(shortestDistances.index[shortestDistances[c]!=0])
                            connectedKmers = set(c[0] for c in connectedOrientatedKmers)
                            #collect paired connected k-mers
                            pairedConnectedKmers = set()
                            for connectedCkmersEntry in self._connectedCkmers:
                                if c[0] in connectedCkmersEntry:
                                    for pairedCkmersEntry in self._connectedCkmers[connectedCkmersEntry]["paired"]:
                                        pairedConnectedKmers.update(pairedCkmersEntry)
                            missingConnections = pairedConnectedKmers.intersection(relevantKmers).difference(connectedKmers)
                            if len(missingConnections)>0:
                                for c3 in new_incoming_kmers:
                                    newConnectedOrientatedKmers = set(shortestDistances.index[shortestDistances[c3]!=0])
                                    newConnectedKmers = set(c[0] for c in newConnectedOrientatedKmers)
                                    #add this connection if it will restore missing connections
                                    if (c3[0] in missingConnections or 
                                        len(newConnectedKmers.intersection(missingConnections))>0):
                                        numberOfFixedMissingConnections+=1
                                        self._orientatedCkmers[c]._setIncoming(c3,distance,0,False) 
                                        if c[0]==c3[0]:
                                            self._orientatedCkmers[c3]._setIncoming(c,distance,0,False) 
                                        else:
                                            self._orientatedCkmers[c3]._setOutgoing(c,distance,0,False) 
            #check for missing outgoing
            if outgoingRelevantSplitSide in c_data._orientatedBases:
                b = c_data._orientatedBases[outgoingRelevantSplitSide]
                assert (b in oRightSplitBases_candidate) or not c_data.candidate()
                if b in oRightSplitBases_candidate:
                    b_data = self._orientatedBases[b]
                    assert c in b_data._orientatedCkmers
                    known_outgoing_kmers = ocKmers_candidateBases.intersection(c_data._outgoing)
                    potential_outgoing_kmers = set()
                    for c2 in ocKmers_candidateBases.intersection(b_data._orientatedCkmers):
                        potential_outgoing_kmers.update(self._orientatedCkmers[c2]._outgoing)
                    potential_outgoing_kmers = ocKmers_candidateBases.intersection(potential_outgoing_kmers)
                    if len(potential_outgoing_kmers.difference(known_outgoing_kmers))>0:
                        new_outgoing_kmers = potential_outgoing_kmers.difference(known_outgoing_kmers)
                        distance = self.getShortestDistance(b_data._orientatedCkmers,potential_outgoing_kmers)
                        #if there are no connections at all, just add all possible edges
                        #this can however introduce false connections
                        if len(known_outgoing_kmers)==0:
                            for c3 in new_outgoing_kmers:
                                numberOfFixedMissingConnections+=1
                                self._orientatedCkmers[c]._setOutgoing(c3,distance,0,False) 
                                if c[0]==c3[0]:
                                    self._orientatedCkmers[c3]._setOutgoing(c,distance,0,False) 
                                else:
                                    self._orientatedCkmers[c3]._setIncoming(c,distance,0,False) 
                        #otherwise, only check if potential incoming can restore missing connections
                        else:
                            #find connected
                            connectedOrientatedKmers = set(shortestDistances.index[shortestDistances[c]!=0])
                            connectedKmers = set(c[0] for c in connectedOrientatedKmers)
                            #collect paired connected k-mers
                            pairedConnectedKmers = set()
                            for connectedCkmersEntry in self._connectedCkmers:
                                if c[0] in connectedCkmersEntry:
                                    for pairedCkmersEntry in self._connectedCkmers[connectedCkmersEntry]["paired"]:
                                        pairedConnectedKmers.update(pairedCkmersEntry)
                            missingConnections = pairedConnectedKmers.intersection(relevantKmers).difference(connectedKmers)
                            if len(missingConnections)>0:
                                for c3 in new_outgoing_kmers:
                                    newConnectedOrientatedKmers = set(shortestDistances.index[shortestDistances[c3]!=0])
                                    newConnectedKmers = set(c[0] for c in newConnectedOrientatedKmers)
                                    #add this connection if it will restore missing connections
                                    if (c3[0] in missingConnections or 
                                        len(newConnectedKmers.intersection(missingConnections))>0):
                                        numberOfFixedMissingConnections+=1
                                        self._orientatedCkmers[c]._setOutgoing(c3,distance,0,False) 
                                        if c[0]==c3[0]:
                                            self._orientatedCkmers[c3]._setOutgoing(c,distance,0,False) 
                                        else:
                                            self._orientatedCkmers[c3]._setIncoming(c,distance,0,False) 
        
        self._logger.debug("automatically included %d direct connections" % (numberOfFixedMissingConnections))
        #reset distances, they should be recomputed with new connections
        if numberOfFixedMissingConnections>0:            
            self._resetDistances()
           
    
    def _glueMissingConnections(self):
        #compute distances
        shortestDistances = self.getShortestDistances()
        #define orientated right split bases containing a candidate
        oRightSplitBases_candidate = set([b for b in self._orientatedBases.keys() if self._orientatedBases[b].candidate()])
        #define orientated canonical k-mers with candidate base
        #flat_list = [item for sublist in t for item in sublist]
        ocKmers_candidateBases = [self._orientatedBases[b]._orientatedCkmers for b in oRightSplitBases_candidate]
        ocKmers_candidateBases = set([item for sublist in ocKmers_candidateBases for item in sublist])
        relevantCkmers = set([c[0] for c in ocKmers_candidateBases])
        #define the set of relevant open orientated k-mers:
        #- open is missing connections in at least one direction
        ocKmers_disconnected = set()
        numberOfGluedMissingConnections = 0
        for c in ocKmers_candidateBases:
            c_data = self._orientatedCkmers[c]
            #at least normal expansion should have been tried
            if not c_data._expanded:
                continue
            baseCkmers = set()
            for direction in c_data._orientatedBases:
                baseCkmers.update(self._orientatedBases[c_data._orientatedBases[direction]]._orientatedCkmers)
            if (len(c_data._incoming)==0 or len(c_data._outgoing)==0):
                ocKmers_disconnected.add(c)
        self._logger.debug("found %d relevant k-mers with missing connections to other candidates" % 
                           (len(ocKmers_disconnected)))  
        
        def glue_missing(c, side, numberOfGluedMissingConnections, maximumDistance: int = 1000, minimumFrequency: int = 2):
            assert side in ["incoming","outgoing"]
            assert c[1] in ["forward","backward"]
            #explore to the right, so get the correct orientated initial k-mer
            if (side=="outgoing" and c[1]=="forward") or (side=="incoming" and c[1]=="backward"):
                initialKmer = c[0]
            else:
                initialKmer = haplotyping.General.reverse_complement(c[0])
            #try to glue
            response = requests.get("{}kmer/{}/{}/split".format(self._baseUrl,self._uid,initialKmer), 
                                    auth=self._apiAuth, headers=self._apiHeaders)
            if response.ok:
                data = response.json()
                if len(data["splittingKmers"])>0:
                    for newKmer in data["splittingKmers"]:
                        ckmer = haplotyping.General.canonical(newKmer)
                        path = data["pathKmers"]
                        distance = len(path)  
                        if (ckmer==newKmer and side == "outgoing") or (ckmer!=newKmer and side=="incoming"):
                            c2 = (ckmer,"forward")
                        else:
                            c2 = (ckmer,"backward")
                        #for now, only process match
                        if not c2 in self._orientatedCkmers:
                            response = requests.get(self._baseUrl+"split/"+self._uid+"/kmer/"+str(c2[0]),
                                                    auth=self._apiAuth, headers=self._apiHeaders)
                            if response.ok:
                                item = response.json()
                                self._createOrientatedCkmer(item, c2[1])
                        if c2 in self._orientatedCkmers:
                            gluePath = "".join([path[0]]+[k[-1] for k in path[1:]]+[newKmer[-1]])
                            if side=="outgoing":
                                numberOfGluedMissingConnections+=1
                                self._orientatedCkmers[c]._setOutgoing(c2,distance,0,False,gluePath)
                                self._orientatedCkmers[c2]._setIncoming(c,distance,0,False,gluePath)
                            else:
                                gluePath = haplotyping.General.reverse_complement(gluePath) 
                                numberOfGluedMissingConnections+=1
                                self._orientatedCkmers[c]._setIncoming(c2,distance,0,False,gluePath)
                                self._orientatedCkmers[c2]._setOutgoing(c,distance,0,False,gluePath)
                elif len(data["pathKmers"])>0:
                    path = data["pathKmers"]
                    ckmer = haplotyping.General.canonical(path[-1])
                    distance = len(path)-1                                                             
                    if (ckmer==path[-1] and side == "outgoing") or (ckmer!=path[-1] and side=="incoming"):
                        c2 = (ckmer,"forward")
                    else:
                        c2 = (ckmer,"backward")
                    gluePath = "".join([path[0]]+[k[-1] for k in path[1:]])
                    if side=="outgoing":
                        self._orientatedCkmers[c]._setOutgoingDeadEnd(c2,distance,gluePath)
                    elif side=="incoming":
                        gluePath = haplotyping.General.reverse_complement(gluePath) 
                        self._orientatedCkmers[c]._setIncomingDeadEnd(c2,distance,gluePath)
            else:
                self._logger.error("request to {} didn't succeed".format(self._baseUrl))
                    
            return numberOfGluedMissingConnections
        
        for c in ocKmers_disconnected:
            c_data = self._orientatedCkmers[c]
            if len(c_data._incoming)==0:
                numberOfGluedMissingConnections = glue_missing(c, "incoming", numberOfGluedMissingConnections)
            if len(c_data._outgoing)==0:
                numberOfGluedMissingConnections = glue_missing(c, "outgoing", numberOfGluedMissingConnections)
                
        self._logger.debug("automatically glued %d direct connections" % (numberOfGluedMissingConnections))
        #reset distances, they should probably be recomputed with new connections
        if numberOfGluedMissingConnections>0:            
            self._resetDistances()
                    
    def _visualize_base_label(self, orientatedBase: str):
        base_label = super(haplotyping.graph.sequence.SequenceGraph, self)._visualize_base_label(orientatedBase)
        if not self._orientatedBases[orientatedBase]._order==None:
            base_label = base_label[:-1]
            base_label +="<font point-size=\"6\">"
            base_label +=" <font color=\"blue\">pos: "+str(self._orientatedBases[orientatedBase]._order)+"</font>"
            base_label +="</font>"
            base_label +=">"
        return base_label
    
    def _visualize_node_label(self, orientatedBase: str, orientatedCkmer: str):
        node_label = super(haplotyping.graph.sequence.SequenceGraph, self)._visualize_node_label(
            orientatedBase,orientatedCkmer)
        if not self._orientatedCkmers[orientatedCkmer]._level==None:
            node_label = node_label[:-1]
            node_label +="<font point-size=\"8\">"
            node_label +=" <font color=\"blue\"> "+str(self._orientatedCkmers[orientatedCkmer]._level)+"</font>"
            node_label +="</font>" + ">"
        return node_label

    class SequenceCkmer(APIGraph.OrientatedCkmer):
        
        def __init__(self, graph: APIGraph, ckmer: str, orientation: str, number: int, split: str):
            """initialise sequence k-mer"""
            #call parent constructor
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).__init__(
                graph, ckmer, orientation, number, split)
            #additional parameters
            self._level : int = None
            self._expanded : bool = False
            self._position : int = None
            self._estimatedPositions : set = set()            
            self._candidateConnected : dict = {"incoming": False, "outgoing": False}
            for direction in self._orientatedBases.keys():
                orientatedBaseKey = self._orientatedBases[direction]
                orientatedBase = self._graph._orientatedBases[orientatedBaseKey]
                for otherOrientatedCkmerKey in orientatedBase._orientatedCkmers:
                    if not self._key==otherOrientatedCkmerKey:
                        if not self._graph._orientatedCkmers[otherOrientatedCkmerKey]._position==None:
                            self._addEstimatedPositions([self._graph._orientatedCkmers[otherOrientatedCkmerKey]._position])
                        else:
                            self._addEstimatedPositions(
                                self._graph._orientatedCkmers[otherOrientatedCkmerKey]._estimatedPositions)
        
        def _setIncoming(self, orientatedCkmerKey: str, distance: int, number: int, problem: bool, path: str = None):
            #call parent method
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self)._setIncoming(
                orientatedCkmerKey,distance,number,problem,path)
            if not self._graph._orientatedCkmers[orientatedCkmerKey]._position==None:
                self._addEstimatedPositions([self._graph._orientatedCkmers[orientatedCkmerKey]._position+distance])
            else:
                self._addEstimatedPositions(
                    [p+distance for p in self._graph._orientatedCkmers[orientatedCkmerKey]._estimatedPositions])
            if not self._level == None:
                self._graph._orientatedCkmers[orientatedCkmerKey]._setLevel(self._level+1)
            if self._incoming[orientatedCkmerKey]["path"]==None:
                #check other direction
                if self._key in self._graph._orientatedCkmers[orientatedCkmerKey]._outgoing.keys():
                    self._incoming[orientatedCkmerKey]["path"] = (
                        self._graph._orientatedCkmers[orientatedCkmerKey]._outgoing[self._key]["path"])
                #try to compute path (can fail, because path not necessarily is covered by k-mer database)
                else:
                    self._incoming[orientatedCkmerKey]["path"] = self._graph._findPath(orientatedCkmerKey, self._key, distance)
            if not (self._graph._ignoreMissingPath and self._incoming[orientatedCkmerKey]["path"]==None):
                if (self._graph._orientatedCkmers[orientatedCkmerKey].candidate() or 
                    self._graph._orientatedCkmers[orientatedCkmerKey]._candidateConnected["incoming"]):
                    self._setConnectedViaCandidate("incoming")
            #consistency check
            if not self._incoming[orientatedCkmerKey]["path"]==None and not self._incoming[orientatedCkmerKey]["problem"]:
                path = self._incoming[orientatedCkmerKey]["path"]
                k = len(self._ckmer)
                assert len(orientatedCkmerKey[0])==k
                assert len(path)==k+distance
                if self._orientation=="forward":
                    assert path[-k:]==self._ckmer
                else:
                    assert path[-k:]==haplotyping.General.reverse_complement(self._ckmer)
                if orientatedCkmerKey[1]=="forward":
                    assert path[0:k]==orientatedCkmerKey[0]
                else:
                    assert path[0:k]==haplotyping.General.reverse_complement(orientatedCkmerKey[0])
            
        def _setOutgoing(self, orientatedCkmerKey: str, distance: int, number: int, problem: bool, path: str = None):
            #call parent method
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self)._setOutgoing(
                orientatedCkmerKey,distance,number,problem,path)
            if not self._graph._orientatedCkmers[orientatedCkmerKey]._position==None:
                self._addEstimatedPositions([self._graph._orientatedCkmers[orientatedCkmerKey]._position-distance])
            else:
                self._addEstimatedPositions(
                    [p-distance for p in self._graph._orientatedCkmers[orientatedCkmerKey]._estimatedPositions])
            if not self._level == None:
                self._graph._orientatedCkmers[orientatedCkmerKey]._setLevel(self._level+1)
            if self._outgoing[orientatedCkmerKey]["path"]==None:
                #check other direction
                if self._key in self._graph._orientatedCkmers[orientatedCkmerKey]._incoming.keys():
                    self._outgoing[orientatedCkmerKey]["path"] = (
                        self._graph._orientatedCkmers[orientatedCkmerKey]._incoming[self._key]["path"])
                #try to compute path (can fail, because path not necessarily is covered by k-mer database)
                else:
                    self._outgoing[orientatedCkmerKey]["path"] = self._graph._findPath(self._key, orientatedCkmerKey, distance)
            if not (self._graph._ignoreMissingPath and self._outgoing[orientatedCkmerKey]["path"]==None):
                if (self._graph._orientatedCkmers[orientatedCkmerKey].candidate() or 
                    self._graph._orientatedCkmers[orientatedCkmerKey]._candidateConnected["outgoing"]):
                    self._setConnectedViaCandidate("outgoing")
            #consistency check
            if not self._outgoing[orientatedCkmerKey]["path"]==None and not self._outgoing[orientatedCkmerKey]["problem"]:
                path = self._outgoing[orientatedCkmerKey]["path"]
                k = len(self._ckmer)
                assert len(orientatedCkmerKey[0])==k
                assert len(path)==k+distance
                if self._orientation=="forward":
                    assert path[0:k]==self._ckmer
                else:
                    assert path[0:k]==haplotyping.General.reverse_complement(self._ckmer)
                if orientatedCkmerKey[1]=="forward":
                    assert path[-k:]==orientatedCkmerKey[0]
                else:
                    assert path[-k:]==haplotyping.General.reverse_complement(orientatedCkmerKey[0])
            
        def _setExpanded(self):
            self._expanded = True
            
        def _setPosition(self, position: int):
            """set position for k-mer"""
            self._position=position
            self._estimatedPositions = set()
            self._setOrder(self._position)
                    
        def _addEstimatedPositions(self, positions: list):
            """add estimations of positions to k-mer"""
            if self._position==None:
                self._estimatedPositions.update([p for p in positions if isinstance(p, int)])
                if len(self._estimatedPositions)>0:
                    self._setOrder(min(self._estimatedPositions))
            
        def _setLevel(self, level: int):
            """update level with new information"""
            #only update if necessary
            if self._level==None or level<self._level:
                self._level=int(level)
                #update levels from connected k-mers
                for incomingOrientatedCkmerKey in self._incoming.keys():
                    self._graph._orientatedCkmers[incomingOrientatedCkmerKey]._setLevel(level+1)
                for outgoingOrientatedCkmerKey in self._outgoing.keys():
                    self._graph._orientatedCkmers[outgoingOrientatedCkmerKey]._setLevel(level+1)

        def _setConnectedViaCandidate(self, connectionDirection: str):
            assert connectionDirection in ["incoming","outgoing"]
            #only update if necessary
            if not self._candidateConnected[connectionDirection] and not self.candidate():
                self._candidateConnected[connectionDirection] = True
                if connectionDirection=="incoming":
                    if self._candidateConnected["outgoing"]:
                        self._setCandidate()
                    for orientatedCkmerKey in self._outgoing.keys():
                        if not (self._graph._ignoreMissingPath and self._outgoing[orientatedCkmerKey]["path"]==None):
                            self._graph._orientatedCkmers[orientatedCkmerKey]._setConnectedViaCandidate(connectionDirection)
                elif connectionDirection=="outgoing":
                    if self._candidateConnected["incoming"]:
                        self._setCandidate()
                    for orientatedCkmerKey in self._incoming.keys():
                        if not (self._graph._ignoreMissingPath and self._incoming[orientatedCkmerKey]["path"]==None):
                            self._graph._orientatedCkmers[orientatedCkmerKey]._setConnectedViaCandidate(connectionDirection)
                    
        def _unsetCandidate(self):
            #only update if necessary
            if self.candidate():
                super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self)._unsetCandidate()
                #note: don't undon connected_via settings (?)
            
        def _setCandidate(self):
            #only update if necessary
            if not self.candidate():
                super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self)._setCandidate()
                for orientatedCkmerKey in self._outgoing.keys(): 
                    if not (self._graph._ignoreMissingPath and self._outgoing[orientatedCkmerKey]["path"]==None):
                        self._graph._orientatedCkmers[orientatedCkmerKey]._setConnectedViaCandidate("incoming")
                for orientatedCkmerKey in self._incoming.keys():
                    if not (self._graph._ignoreMissingPath and self._incoming[orientatedCkmerKey]["path"]==None):
                        self._graph._orientatedCkmers[orientatedCkmerKey]._setConnectedViaCandidate("outgoing")
            #include additional k-mer for base
            for direction in self._orientatedBases.keys():
                orientatedBaseKey = self._orientatedBases[direction]
                orientatedBase = self._graph._orientatedBases[orientatedBaseKey]
                if len(orientatedBase._orientatedCkmers)<2:
                    response = requests.get(self._graph._baseUrl+"split/"+self._graph._uid+"/base/"+str(orientatedBaseKey[0]),
                                            auth=self._graph._apiAuth, headers=self._graph._apiHeaders)
                    if response.ok:
                        data = response.json()
                        for entry in data["ckmers"]:
                            if not entry["ckmer"]==self._ckmer:
                                for direction in entry["rightSplitBase"].keys():
                                    if entry["rightSplitBase"][direction]==orientatedBaseKey[0]:
                                        assert orientatedBaseKey[1] in ["forward","backward"]
                                        if direction=="right":
                                            entryDirection = orientatedBaseKey[1]
                                        else: 
                                            entryDirection = ("forward" if 
                                                      orientatedBaseKey[1]!="forward" else "backward") 
                                        self._graph._createOrientatedCkmer(entry,entryDirection)
                                
                
                                            
            

        
    