import logging,requests
from haplotyping.graph.api import APIGraph
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
        self._sequence_construct_graph()
        
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
        self._bases = {}
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
                self._logger.debug("found %d splitting k-mers in sequence" % (len(self._ckmers)))
                #find first and last splitting k-mer
                self._logger.debug("define k-mer with position %d as start k-mer" % (position_min))
                self._logger.debug("define k-mer with position %d as end k-mer" % (position_max))
                self._start.update([ckmer for ckmer in self._ckmers.keys() if position_min in self._ckmers[ckmer]._positions])
                self._end.update([ckmer for ckmer in self._ckmers.keys() if position_max in self._ckmers[ckmer]._positions])
                #expand k-mer set from the connected set
                candidate_kmer_keys = [ckmerKey for ckmerKey in self._ckmers.keys() if self._ckmers[ckmerKey]._candidate]
                candidate_kmers = set([ckmerKey[0] for ckmerKey in candidate_kmer_keys])
                #update connected to initial candidates
                response_connected = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/connected", 
                                         json = {"kmers": list(candidate_kmers)})
                if response_connected.ok:
                    data_connected = response_connected.json()
                    connected = set(candidate_kmers)
                    for id in data_connected.keys():
                        connected.update(data_connected[id]["ckmers"])
                    self._logger.debug("initially try to connect %d k-mers" % (len(connected)))                 
                    for level in range(10):
                        self._expand_ckmers(level, connected)
                        new_candidate_kmer_keys=[ckmerKey for ckmerKey in self._ckmers.keys() 
                                             if self._ckmers[ckmerKey]._candidate and not ckmerKey[0] in candidate_kmers]
                        new_candidate_kmers = set([ckmerKey[0] for ckmerKey in new_candidate_kmer_keys])
                        if len(new_candidate_kmers)>0:
                            #update connected with new candidates
                            connected.update(new_candidate_kmers)
                            candidate_kmers.update(new_candidate_kmers)
                            response_connected = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/connected", 
                                                 json = {"kmers": list(new_candidate_kmers)})
                            if response_connected.ok:
                                data_connected = response_connected.json()
                                for id in data_connected.keys():
                                    connected.update(data_connected[id]["ckmers"])
                            else:
                                self._logger.error("request to "+self.baseUrl+" didn't succeed") 
                            self._logger.debug("now try to connect %d k-mers" % (len(connected)))   
                else:
                    self._logger.error("request to "+self.baseUrl+" didn't succeed") 
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed")      
            
    def _expand_ckmers(self, level: int, connected: list):
        """
        expand k-mers from defined level by finding direct connections
        """
        assert level>=0
        #define list of k-mers to check
        kmer_keys = [ckmer_key for ckmer_key in self._ckmers.keys() if self._ckmers[ckmer_key]._level==level]
        kmers = set([ckmer_key[0] for ckmer_key in kmer_keys])
        initial_candidates = sum(self._ckmers[ckmer_key]._candidate for ckmer_key in self._ckmers)
        self._logger.debug("from %d k-mers with %d candidates, expand %d at level %d" % 
                           (len(self._ckmers), initial_candidates, len(kmers), level)) 
        #get direct connections
        response = requests.post(self._baseUrl+"split/"+self._uid+"/kmer/direct", 
                                 json = {"kmers": list(kmers)})
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
                    if ckmerKey1 in self._ckmers.keys():
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
                                        if not ckmerKey2 in self._ckmers.keys():
                                            self._create_ckmer(directionItem, orientation2)
                                            set_new_kmers.add(ckmer2)
                                        self._ckmers[ckmerKey1].set_outgoing(ckmerKey2, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)
                                        self._ckmers[ckmerKey2].set_incoming(ckmerKey1, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)
                                        self._ckmers[ckmerKey2].set_level(level+1)
                                    else:
                                        #inward
                                        orientation2 = "forward" if direction2=="right" else "backward"
                                        ckmerKey2 = (ckmer2,orientation2) 
                                        if not ckmerKey2 in self._ckmers.keys():
                                            self._create_ckmer(directionItem, orientation2)
                                            set_new_kmers.add(ckmer2)
                                        self._ckmers[ckmerKey1].set_incoming(ckmerKey2, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)  
                                        self._ckmers[ckmerKey2].set_outgoing(ckmerKey1, 
                                                                             directionItem["connection"]["distance"], 
                                                                             directionItem["connection"]["number"],
                                                                             directionItem["connection"]["problem"]>0)  
                                        self._ckmers[ckmerKey2].set_level(level+1)
            final_candidates = sum(self._ckmers[ckmer]._candidate for ckmer in self._ckmers)
            self._logger.debug("skipped %d and found %d connections to %d k-mers: %d new and %d new candidates" % 
                               (counter_skipped_connections, counter_direct_connections,len(set_direct_kmers),
                                len(set_new_kmers),(final_candidates-initial_candidates))) 
                                    
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed")
            
    def _create_ckmer(self, item: dict, orientation: str):         
        """
        create k-mer and add automatically base(s)
        """
        assert "ckmer" in item.keys() and "number" in item.keys() and "split" in item.keys()
        assert orientation in ["forward","backward"]
        ckmer = item["ckmer"]
        ckmer_key = (ckmer, orientation)
        if not ckmer_key in self._ckmers.keys():
            ckmerEntry = self.SequenceCkmer(self, ckmer, orientation, item["number"], item["split"])
            for direction in item["rightSplitBase"]:
                base = item["rightSplitBase"][direction]["base"]                
                if not base in self._bases.keys():
                    baseEntry = self.SequenceBase(self, base, 
                                            item["rightSplitBase"][direction]["number"], 
                                            item["rightSplitBase"][direction]["ckmers"]) 
        else:
            self._logger.warning("k-mer already exists")
        return self._ckmers[ckmer_key]
            
    class SequenceCkmer(APIGraph.Ckmer):
        
        def __init__(self, graph: APIGraph, ckmer: str, orientation: str, number: int, split: str):
            """initialise sequence k-mer"""
            #call parent constructor
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).__init__(
                graph, ckmer, orientation, number, split)
            #additional parameters
            self._level : int = None
            self._positions : set = set()            
            self._candidateConnected : dict = {"incoming": False, "outgoing": False} 
        
        def set_incoming(self, ckmerKey: str, distance: int, number: int, problem: bool):
            #call parent method
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).set_incoming(
                ckmerKey,distance,number,problem)
            if self._graph._ckmers[ckmerKey]._candidate or self._graph._ckmers[ckmerKey]._candidateConnected["incoming"]:
                self._set_connected_via_candidate("incoming")
            
        def set_outgoing(self, ckmerKey: str, distance: int, number: int, problem: bool):
            #call parent method
            super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).set_outgoing(
                ckmerKey,distance,number,problem)
            if self._graph._ckmers[ckmerKey]._candidate or self._graph._ckmers[ckmerKey]._candidateConnected["outgoing"]:
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
                for incomingCkmerKey in self._incoming.keys():
                    self._graph._ckmers[incomingCkmerKey].set_level(level+1)
                for outgoingCkmerKey in self._outgoing.keys():
                    self._graph._ckmers[outgoingCkmerKey].set_level(level+1)

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
                    for ckmerKey in self._outgoing.keys():
                        self._graph._ckmers[ckmerKey]._set_connected_via_candidate(connectionDirection)
                elif connectionDirection=="outgoing":
                    if self._candidateConnected["incoming"]:
                        self.set_candidate()
                    for ckmerKey in self._incoming.keys():
                        self._graph._ckmers[ckmerKey]._set_connected_via_candidate(connectionDirection)

        def set_candidate(self):
            #only update if necessary
            if not self._candidate:
                super(haplotyping.graph.sequence.SequenceGraph.SequenceCkmer, self).set_candidate()
                                            
            
    class SequenceBase(APIGraph.Base):
        
        def __init__(self, graph: APIGraph, base: str, number: int, ckmers: dict):
            #call parent constructor
            super(haplotyping.graph.sequence.SequenceGraph.SequenceBase, self).__init__(graph, base, number, ckmers) 
        
    