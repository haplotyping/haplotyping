import logging,requests,re
from haplotyping.graph.api import APIGraph
from haplotyping.general import General
import haplotyping

class SequenceGraph(APIGraph):
    """constructing the De Bruijn graph"""
    
    def __init__(self, url: str, uid: str, sequence: str):
        """
        construct the De Bruijn graph from a reference sequence using the dataset 
        defined by uid via the API at provided url
        """
        #call parent constructor
        super(haplotyping.graph.sequence.SequenceGraph, self).__init__(url, uid)
        #logger
        self._logger = logging.getLogger(__name__)
        #store call parameters
        self._sequence = str(sequence).strip()
        self._logger.info("construct De Bruijn graph for sequence of length "+str(len(self._sequence)))
        #start initialisation
        self._sequencePositions = set()
        self._sequence_construct_graph()
        #try to fix missing connections in the graph
        self._fix_missing_connections()
        #try to glue missing connections
        self._glue_missing_connections()
        
    def __repr__(self):
        text = super(haplotyping.graph.sequence.SequenceGraph, self).__repr__()
        text = text + " for sequence of length %d" % (len(self._sequence))
        return text

    def _sequence_construct_graph(self):
        """
        main function to construct the De Bruijn Graph
        """
        self._logger.debug("start constructing the De Bruijn Graph")  
        #initialise
        self._start = set()
        self._end = set()
        self._ckmers = {}
        self._orientatedCkmers = {}
        self._orientatedBases = {}
        self._connectedCkmers = {}
        #get main k-mers from sequence
        response = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/sequence", json = {"sequence": self._sequence})
        if response.ok:
            data = response.json()
            if len(data)==0:
                self._logger.error("no splitting k-mers found in sequence")
            else:
                #pre-compute positions
                forward_positions = {}
                backward_positions = {}
                for i in range(len(self._sequence)-(self._k-1)):
                    kmer = self._sequence[i:i+self._k]
                    ckmer = haplotyping.General.canonical(kmer)
                    if kmer==ckmer:
                        if ckmer in forward_positions.keys():
                            forward_positions[ckmer].append(i)
                        else:
                            forward_positions[ckmer]=[i]
                    else:
                        if ckmer in backward_positions.keys():
                            backward_positions[ckmer].append(i)
                        else:
                            backward_positions[ckmer]=[i]
                #initialise minimum and maximum
                position_min = len(self._sequence)-self._k
                position_max = 0
                #create initial k-mers
                for item in data:
                    #get positions
                    ckmer_forward_positions = set(forward_positions.get(item["ckmer"],[]))
                    ckmer_backward_positions = set(backward_positions.get(item["ckmer"],[]))
                    self._sequencePositions.update(ckmer_forward_positions)
                    self._sequencePositions.update(ckmer_backward_positions)
                    if len(ckmer_forward_positions)>0:
                        ckmerForwardEntry = self._create_ckmer(item, "forward")
                        ckmerForwardEntry.set_level(0)
                        ckmerForwardEntry.set_candidate()
                        position_min = min(position_min,min(ckmer_forward_positions))
                        position_max = max(position_max,max(ckmer_forward_positions))
                        ckmerForwardEntry.add_positions(ckmer_forward_positions)
                    if len(ckmer_backward_positions)>0:
                        ckmerBackwardEntry = self._create_ckmer(item, "backward")
                        ckmerBackwardEntry.set_level(0)
                        ckmerBackwardEntry.set_candidate()
                        position_min = min(position_min,min(ckmer_backward_positions))
                        position_max = max(position_max,max(ckmer_backward_positions))
                        ckmerBackwardEntry.add_positions(ckmer_backward_positions)
                self._logger.debug("found %d splitting k-mers in sequence" % (len(self._orientatedCkmers)))
                #find first and last splitting k-mer
                self._logger.debug("define k-mer with position %d as start k-mer" % (position_min))
                self._logger.debug("define k-mer with position %d as end k-mer" % (position_max))
                self._start.update([k for k in self._orientatedCkmers.keys() 
                                    if position_min in self._orientatedCkmers[k]._positions])
                self._end.update([k for k in self._orientatedCkmers.keys() 
                                  if position_max in self._orientatedCkmers[k]._positions])
                #expand k-mer set from the connected set
                candidate_kmer_keys = [k for k in self._orientatedCkmers.keys() 
                                       if self._orientatedCkmers[k]._candidate]
                candidate_kmers = set([k[0] for k in candidate_kmer_keys])
                connected = set(candidate_kmers)
                #local helper function to find connected k-mers
                def _find_connected(kmer_list):
                    response_connected = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/connected", 
                                         json = {"kmers": sorted(list(kmer_list))})
                    if response_connected.ok:
                        data_connected = response_connected.json()
                        for id in data_connected.keys():
                            connected.update(data_connected[id]["ckmers"])
                            connectedCkmersEntry = tuple(data_connected[id]["ckmers"])
                            if not connectedCkmersEntry in self._connectedCkmers:
                                self._connectedCkmers[connectedCkmersEntry] = {
                                    "paired": {}, 
                                    "direct": data_connected[id]["direct"]>0,
                                    "number": data_connected[id]["number"]
                                }
                            for pairedId in data_connected[id]["paired"]:
                                pairedConnectedCkmersEntry = tuple(data_connected[pairedId]["ckmers"])
                                self._connectedCkmers[connectedCkmersEntry]["paired"][pairedConnectedCkmersEntry] = {
                                    "number": data_connected[id]["paired"][pairedId]
                                } 
                    else:
                        self._logger.error("request to "+self.baseUrl+" didn't succeed") 
                #expand k-mers
                for level in range(10):                    
                    if level==0:
                        new_candidate_kmers = candidate_kmers
                    else:
                        new_candidate_kmer_keys=[k for k in self._orientatedCkmers.keys() 
                                         if self._orientatedCkmers[k]._candidate and not k[0] in candidate_kmers]
                        new_candidate_kmers = set([k[0] for k in new_candidate_kmer_keys])
                    if len(new_candidate_kmers)>0:
                        _find_connected(candidate_kmers)
                        self._logger.debug("try to connect %d k-mers" % (len(connected)))                        
                    self._expand_ckmers(level, connected)                    
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed") 

    def _expand_ckmers(self, level: int, connected: list):
        """
        expand k-mers from defined level by finding direct connections
        used in _sequence_construct_graph()
        """
        assert level>=0
        #define list of k-mers to check
        kmer_keys = [k for k in self._orientatedCkmers.keys() 
                     if self._orientatedCkmers[k]._level==level]
        kmers = set([k[0] for k in kmer_keys])
        initial_candidates = sum(self._orientatedCkmers[ckmer_key]._candidate for ckmer_key in self._orientatedCkmers)
        self._logger.debug("from %d k-mers with %d candidates, expand %d at level %d" % 
                           (len(self._orientatedCkmers), initial_candidates, len(kmers), level)) 
        #get direct connections
        response = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/direct", 
                                 json = {"kmers": sorted(list(kmers))})
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
                                if not ckmer2 in connected:
                                    counter_skipped_connections+=1
                                else:
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
                                            self._create_ckmer(directionItem, orientation2)
                                            set_new_kmers.add(ckmer2)
                                        self._orientatedCkmers[ckmerKey1].set_outgoing(ckmerKey2, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)
                                        self._orientatedCkmers[ckmerKey2].set_incoming(ckmerKey1, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)
                                        self._orientatedCkmers[ckmerKey2].set_level(level+1)
                                    else:
                                        #inward
                                        orientation2 = "forward" if direction2=="right" else "backward"
                                        ckmerKey2 = (ckmer2,orientation2) 
                                        if not ckmerKey2 in self._orientatedCkmers.keys():
                                            self._create_ckmer(directionItem, orientation2)
                                            set_new_kmers.add(ckmer2)
                                        self._orientatedCkmers[ckmerKey1].set_incoming(ckmerKey2, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)  
                                        self._orientatedCkmers[ckmerKey2].set_outgoing(ckmerKey1, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)  
                                        self._orientatedCkmers[ckmerKey2].set_level(level+1)
            final_candidates = sum(self._orientatedCkmers[ckmer]._candidate for ckmer in self._orientatedCkmers)
            self._logger.debug("skipped %d and found %d connections to %d k-mers: %d new and %d new candidates" % 
                               (counter_skipped_connections, counter_direct_connections,len(set_direct_kmers),
                                len(set_new_kmers),(final_candidates-initial_candidates))) 
                                    
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed")
            
    def _create_ckmer(self, item: dict, orientation: str):         
        """
        create k-mer, used in _sequence_construct_graph
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
    
    def _fix_missing_connections(self):
        #compute k-mers with only single orientation (not necessary to extend this to double orientation?)
        relevantKmers = set([c for c in self._ckmers if len(self._ckmers[c])==1])
        #compute distances
        shortestDistances = self.getShortestDistances()
        #define orientated canonical candidate k-mers
        ocKmers_candidate = set([c for c in self._orientatedCkmers.keys() if self._orientatedCkmers[c]._candidate])
        #define orientated right split bases containing a candidate
        oRightSplitBases_candidate = set([b for b in self._orientatedBases.keys() if self._orientatedBases[b]._candidate])
        #define orientated canonical k-mers with candidate base
        #flat_list = [item for sublist in t for item in sublist]
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
                assert (b in oRightSplitBases_candidate) or not c_data._candidate
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
                                self._orientatedCkmers[c].set_incoming(c3,distance,0,False) 
                                if c[0]==c3[0]:
                                    self._orientatedCkmers[c3].set_incoming(c,distance,0,False) 
                                else:
                                    self._orientatedCkmers[c3].set_outgoing(c,distance,0,False) 
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
                                        self._orientatedCkmers[c].set_incoming(c3,distance,0,False) 
                                        if c[0]==c3[0]:
                                            self._orientatedCkmers[c3].set_incoming(c,distance,0,False) 
                                        else:
                                            self._orientatedCkmers[c3].set_outgoing(c,distance,0,False) 
            #check for missing outgoing
            if outgoingRelevantSplitSide in c_data._orientatedBases:
                b = c_data._orientatedBases[outgoingRelevantSplitSide]
                assert (b in oRightSplitBases_candidate) or not c_data._candidate
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
                                self._orientatedCkmers[c].set_outgoing(c3,distance,0,False) 
                                if c[0]==c3[0]:
                                    self._orientatedCkmers[c3].set_outgoing(c,distance,0,False) 
                                else:
                                    self._orientatedCkmers[c3].set_incoming(c,distance,0,False) 
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
                                        self._orientatedCkmers[c].set_outgoing(c3,distance,0,False) 
                                        if c[0]==c3[0]:
                                            self._orientatedCkmers[c3].set_outgoing(c,distance,0,False) 
                                        else:
                                            self._orientatedCkmers[c3].set_incoming(c,distance,0,False) 
        
        self._logger.debug("automatically included %d direct connections" % (numberOfFixedMissingConnections))
        #reset distances, they should be recomputed with new connections
        if numberOfFixedMissingConnections>0:            
            self._resetDistances()
           
    
    def _glue_missing_connections(self):
        #compute distances
        shortestDistances = self.getShortestDistances()
        #define orientated right split bases containing a candidate
        oRightSplitBases_candidate = set([b for b in self._orientatedBases.keys() if self._orientatedBases[b]._candidate])
        #define orientated canonical k-mers with candidate base
        #flat_list = [item for sublist in t for item in sublist]
        ocKmers_candidateBases = [self._orientatedBases[b]._orientatedCkmers for b in oRightSplitBases_candidate]
        ocKmers_candidateBases = set([item for sublist in ocKmers_candidateBases for item in sublist])
        relevantCkmers = set([c[0] for c in ocKmers_candidateBases])
        #define the set of relevant open orientated k-mers:
        #- open is missing connections in at least one direction
        #- relevant if no start or end, and in reference sequence
        ocKmers_disconnected = set()
        for c in ocKmers_candidateBases:
            #skip start and end
            if c in self._start or c in self._end:
                continue            
            c_data = self._orientatedCkmers[c]
            baseCkmers = set()
            for direction in c_data._orientatedBases:
                baseCkmers.update(self._orientatedBases[c_data._orientatedBases[direction]]._orientatedCkmers)
            #skip start and end base
            if len(self._start.intersection(baseCkmers))>0 or len(self._end.intersection(baseCkmers))>0:
                continue
            if (len(c_data._incoming)==0 or len(c_data._outgoing)==0):
                self._selected.add(c)
                ocKmers_disconnected.add(c)
        self._logger.debug("found %d relevant k-mers with missing connections to other candidates" % 
                           (len(ocKmers_disconnected)))  
        
        def glue_missing(c, side, maximumDistance: int = 1000, minimumFrequency: int = 2):
            assert side in ["incoming","outgoing"]
            assert c[1] in ["forward","backward"]
            print("GLUE "+str(c))
            #explore to the right, so get the correct orientated initial k-mer
            if (side=="outgoing" and c[1]=="forward") or (side=="incoming" and c[1]=="backward"):
                initialKmer = c[0]
            else:
                initialKmer = haplotyping.General.reverse_complement(c[0])
            #try to glue
            newKmer = initialKmer
            path = [newKmer]
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
                        ckmer = haplotyping.General.canonical(newKmer)
                        distance = len(path) - 1
                        if (ckmer==newKmer and side == "outgoing") or (ckmer!=newKmer and side=="incoming"):
                            c2 = (ckmer,"forward")
                        else:
                            c2 = (ckmer,"backward")
                        #for now, only process match
                        if c2 in self._orientatedCkmers:
                            if side=="outgoing":
                                self._orientatedCkmers[c].set_outgoing(c2,distance,0,False)
                                self._orientatedCkmers[c2].set_incoming(c,distance,0,False)
                            else:
                                self._orientatedCkmers[c].set_incoming(c2,distance,0,False)
                                self._orientatedCkmers[c2].set_outgoing(c,distance,0,False)
                        break
                    else:
                        rightSplitters = list(rightNeighbours.intersection(kmerFound))
                        if len(rightSplitters)>1:
                            for newKmer in rightSplitters:
                                ckmer = haplotyping.General.canonical(newKmer)
                                distance = len(path)
                                if (ckmer==newKmer and side == "outgoing") or (ckmer!=newKmer and side=="incoming"):
                                    c2 = (ckmer,"forward")
                                else:
                                    c2 = (ckmer,"backward")
                                #for now, only process match
                                if c2 in self._orientatedCkmers:
                                    if side=="outgoing":
                                        self._orientatedCkmers[c].set_outgoing(c2,distance,0,False)
                                        self._orientatedCkmers[c2].set_incoming(c,distance,0,False)
                                    else:
                                        self._orientatedCkmers[c].set_incoming(c2,distance,0,False)
                                        self._orientatedCkmers[c2].set_outgoing(c,distance,0,False)
                            break
                        elif len(rightSplitters)==0:
                            break
                        else:
                            newKmer = rightSplitters[0]   
                            path.append(newKmer)
                else:
                    self._logger.error("request to "+self.baseUrl+" didn't succeed")            
        
        for c in ocKmers_disconnected:
            c_data = self._orientatedCkmers[c]
            if len(c_data._incoming)==0:
                glue_missing(c, "incoming")
            if len(c_data._outgoing)==0:
                glue_missing(c, "outgoing")
        
        return
        
        if len(selectedSequencePositions)<0:
            assert self._sequencePositions.issuperset(selectedSequencePositions)
            selectedSequencePositions = sorted(list(selectedSequencePositions))
            for i in range(1,len(selectedSequencePositions)):
                j = allSequencePositions.index(selectedSequencePositions[i])
                assert j>0
                if allSequencePositions[j-1]==selectedSequencePositions[i-1]:
                    kmer1 = self._sequence[selectedSequencePositions[i-1]:selectedSequencePositions[i-1]+self._k]
                    kmer2 = self._sequence[selectedSequencePositions[i]:selectedSequencePositions[i]+self._k]
                    if kmer1==haplotyping.General.canonical(kmer1):
                        c1 = (kmer1,"forward")
                    else:
                        c1 = (haplotyping.General.canonical(kmer1),"backward")
                    assert c1 in ocKmers_openRelevant
                    if kmer2==haplotyping.General.canonical(kmer2):
                        c2 = (kmer2,"forward")
                    else:
                        c2 = (haplotyping.General.canonical(kmer2),"backward")   
                    assert c2 in ocKmers_openRelevant
                    distance = selectedSequencePositions[i] - selectedSequencePositions[i-1]
                    #todo: verify existence and open right side (?)
                    
                    self.findRightConnection(kmer1,distance+10)
                    
                    self.findRightConnection(haplotyping.General.reverse_complement(kmer2),distance+10)
                    
                    self._orientatedCkmers[c1].set_outgoing(c2,distance,0,False) 
                    self._orientatedCkmers[c2].set_incoming(c1,distance,0,False) 
                    
                    
                
    
    class SequenceCkmer(APIGraph.OrientatedCkmer):
        
        def __init__(self, graph: APIGraph, ckmer: str, orientation: str, number: int, split: str):
            """initialise sequence k-mer"""
            #call parent constructor
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).__init__(
                graph, ckmer, orientation, number, split)
            #additional parameters
            self._level : int = None
            self._positions : set = set()            
            self._candidateConnected : dict = {"incoming": False, "outgoing": False} 
        
        def set_incoming(self, orientatedCkmerKey: str, distance: int, number: int, problem: bool):
            #call parent method
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).set_incoming(
                orientatedCkmerKey,distance,number,problem)
            if (self._graph._orientatedCkmers[orientatedCkmerKey]._candidate or 
                self._graph._orientatedCkmers[orientatedCkmerKey]._candidateConnected["incoming"]):
                self._set_connected_via_candidate("incoming")
            
        def set_outgoing(self, orientatedCkmerKey: str, distance: int, number: int, problem: bool):
            #call parent method
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).set_outgoing(
                orientatedCkmerKey,distance,number,problem)
            if (self._graph._orientatedCkmers[orientatedCkmerKey]._candidate or 
                self._graph._orientatedCkmers[orientatedCkmerKey]._candidateConnected["outgoing"]):
                self._set_connected_via_candidate("outgoing")
            
        def add_positions(self, positions: list):
            """add positions to k-mer"""
            self._positions.update([p for p in positions if isinstance(p, int)])
            if len(self._positions)>0:
                self.set_order(min(self._positions))
            
        def set_level(self, level: int):
            """update level with new information"""
            #only update if necessary
            if self._level==None or level<self._level:
                if not self._level==None:
                    self._graph._logger.debug("change level from "+str(self._level)+" to "+str(level))
                self._level=int(level)
                #update levels from connected k-mers
                for incomingOrientatedCkmerKey in self._incoming.keys():
                    self._graph._orientatedCkmers[incomingOrientatedCkmerKey].set_level(level+1)
                for outgoingOrientatedCkmerKey in self._outgoing.keys():
                    self._graph._orientatedCkmers[outgoingOrientatedCkmerKey].set_level(level+1)

        def _set_connected_via_candidate(self, connectionDirection: str):
            assert connectionDirection in ["incoming","outgoing"]
            #only update if necessary
            if not self._candidateConnected[connectionDirection] and not self._candidate:
                if connectionDirection=="incoming":
                    assert len(self._incoming)>0
                elif connectionDirection=="outgoing":
                    assert len(self._outgoing)>0
                self._candidateConnected[connectionDirection] = True
                if connectionDirection=="incoming":
                    if self._candidateConnected["outgoing"]:
                        self.set_candidate()
                    for orientatedCkmerKey in self._outgoing.keys():
                        self._graph._orientatedCkmers[orientatedCkmerKey]._set_connected_via_candidate(connectionDirection)
                elif connectionDirection=="outgoing":
                    if self._candidateConnected["incoming"]:
                        self.set_candidate()
                    for orientatedCkmerKey in self._incoming.keys():
                        self._graph._orientatedCkmers[orientatedCkmerKey]._set_connected_via_candidate(connectionDirection)

        def set_candidate(self):
            #only update if necessary
            if not self._candidate:
                super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).set_candidate()
                                            
            

        
    