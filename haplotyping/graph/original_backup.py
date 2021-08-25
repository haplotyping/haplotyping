import logging,requests
import graphviz
import haplotyping

class Graph:
    """constructing the De Bruijn graph"""
    
    def __init__(self, url: str, uid: str, sequence: str):
        """construct the De Bruijn graph from a reference sequence using the dataset 
        defined by uid via the API at provided url"""
        #logger
        self._logger = logging.getLogger(__name__)
        #store call parameters
        self.baseUrl = str(url)
        self.uid = str(uid)
        self.sequence = str(sequence).strip()
        self.k = None
        self._logger.info("construct De Bruijn graph for sequence of length "+str(len(self.sequence)))
        self._logger.debug("using dataset "+str(self.uid))
        self._logger.debug("calling API at "+str(self.baseUrl))
        #start initialisation
        self._initial_checks()
        self._construct_graph()
        
    def __repr__(self):
        text = "Graph object for sequence of length %d in dataset" % (len(self.sequence))
        if "variety" in self.dataset and self.dataset["variety"]:
            if "name" in self.dataset["variety"] and self.dataset["variety"]["name"]:
                text = text + " " + str(self.dataset["variety"]["name"])
                additional = []
                if "year" in self.dataset["variety"] and self.dataset["variety"]["year"]:
                    additional.append(str(self.dataset["variety"]["year"]))
                if "origin" in self.dataset["variety"] and self.dataset["variety"]["origin"]:
                    additional.append(str(self.dataset["variety"]["origin"]))
                if len(additional)>0:
                    text = text + "(" + ",".join(additional)+")"
        text = text + " from collection "+str(self.dataset.get("collection","???"))
        return text
            
        
    def _initial_checks(self):
        """
        get dataset and check for k-mer and split
        """
        self._logger.debug("get dataset information from API")        
        #get dataset info
        response = requests.post(self.baseUrl+"dataset/", json = {"uids": [self.uid]})
        if response.ok:
            data = response.json()
            if data.get("total",0)==1:
                self.dataset = data["list"][0]
                if not self.dataset["kmer"]:
                    self._logger.error("dataset "+self.uid+" has no k-mer database")
                if not self.dataset["split"]:
                    self._logger.error("dataset "+self.uid+" has no split k-mer database")    
            else:
                self._logger.error("dataset "+self.uid+" not found")   
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed")
            
    def _construct_graph(self):
        """
        main function to construct the De Bruijn Graph
        """
        self._logger.debug("start constructing the De Bruijn Graph")  
        #initialise
        self.start = set()
        self.end = set()
        self.ckmers = {}
        self.bases = {}
        #get split info
        response = requests.get(self.baseUrl+"split/"+self.uid+"/info")
        if response.ok:
            data = response.json()
            self.k = int(data["k"])
            self._logger.debug("detected k-mer size %d" % (self.k))
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed") 
        #get main k-mers from sequence
        response = requests.post(self.baseUrl+"split/"+self.uid+"/kmer/sequence", json = {"sequence": self.sequence})
        if response.ok:
            data = response.json()
            if len(data)==0:
                self._logger.error("no splitting k-mers found in sequence")
            else:
                #pre-compute positions
                forward_positions = {}
                backward_positions = {}
                for i in range(len(self.sequence)-(self.k-1)):
                    kmer = self.sequence[i:i+self.k]
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
                position_min = len(self.sequence)-self.k
                position_max = 0
                #create initial k-mers
                for item in data:
                    #create k-mer
                    ckmerEntry = self._create_ckmer(item)
                    ckmerEntry.set_level(0)
                    ckmerEntry.set_candidate()
                    #get positions
                    ckmer_forward_positions = set(forward_positions.get(ckmerEntry.ckmer,[]))
                    ckmer_backward_positions = set(backward_positions.get(ckmerEntry.ckmer,[]))
                    #set orientation
                    if len(ckmer_forward_positions)>0 and len(ckmer_backward_positions)>0:
                        ckmerEntry.set_orientation("both")
                        position_min = min(position_min,min(min(ckmer_forward_positions),min(ckmer_backward_positions)))
                        position_max = max(position_max,max(max(ckmer_forward_positions),max(ckmer_backward_positions)))
                    elif len(ckmer_forward_positions)>0:
                        ckmerEntry.set_orientation("forward")
                        position_min = min(position_min,min(ckmer_forward_positions))
                        position_max = max(position_max,max(ckmer_forward_positions))
                    elif len(ckmer_backward_positions)>0:
                        ckmerEntry.set_orientation("backward")
                        position_min = min(position_min,min(ckmer_backward_positions))
                        position_max = max(position_max,max(ckmer_backward_positions))
                    else:
                        ckmerEntry.set_orientation("unknown")
                        self._logger.warning("this should not happen: orientation k-mer unknown at level 0")
                    #set positions
                    ckmerEntry.add_positions(ckmer_forward_positions)
                    ckmerEntry.add_positions(ckmer_backward_positions)
                self._logger.debug("found %d splitting k-mers in sequence" % (len(self.ckmers)))
                #find first and last splitting k-mer
                self._logger.debug("define k-mer with position %d as start k-mer" % (position_min))
                self._logger.debug("define k-mer with position %d as end k-mer" % (position_max))
                self.start.update([ckmer for ckmer in self.ckmers.keys() if position_min in self.ckmers[ckmer].positions])
                self.end.update([ckmer for ckmer in self.ckmers.keys() if position_max in self.ckmers[ckmer].positions])
                #expand k-mer set from the connected set
                candidate_kmers = [ckmer for ckmer in self.ckmers.keys() if self.ckmers[ckmer].candidate]
                #update connected to initial candidates
                response_connected = requests.post(self.baseUrl+"split/"+self.uid+"/kmer/connected", 
                                         json = {"kmers": list(candidate_kmers)})
                if response_connected.ok:
                    data_connected = response_connected.json()
                    connected = set(candidate_kmers)
                    for id in data_connected.keys():
                        connected.update(data_connected[id]["ckmers"])
                    self._logger.debug("initially try to connect %d k-mers" % (len(connected)))                 
                    for level in range(10):
                        self._expand_ckmers(level, connected)
                        new_candidate_kmers=[ckmer for ckmer in self.ckmers.keys() 
                                             if self.ckmers[ckmer].candidate and not ckmer in candidate_kmers]
                        if len(new_candidate_kmers)>0:
                            #update connected with new candidates
                            connected.update(new_candidate_kmers)
                            candidate_kmers.extend(new_candidate_kmers)
                            response_connected = requests.post(self.baseUrl+"split/"+self.uid+"/kmer/connected", 
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
        #define list of k-mers to check
        kmers = [ckmer for ckmer in self.ckmers.keys() if self.ckmers[ckmer].level==level]
        initial_candidates = sum(self.ckmers[ckmer].candidate for ckmer in self.ckmers)
        self._logger.debug("from %d k-mers with %d candidates, expand %d at level %d" % 
                           (len(self.ckmers), initial_candidates, len(kmers), level)) 
        #get direct connections
        response = requests.post(self.baseUrl+"split/"+self.uid+"/kmer/direct", json = {"kmers": kmers})
        counter_direct_connections = 0
        counter_skipped_connections = 0
        set_direct_kmers = set()
        set_new_kmers = set()
        if response.ok:
            data = response.json()
            for item in data:
                fromCkmer = item["ckmer"]
                if not fromCkmer in self.ckmers.keys():
                    self._logger.error("unexpected missing k-mer "+fromCkmer)
                else:                    
                    for direction in item["direct"]:
                        if direction in ["left","right"]:
                            for directionItem in item["direct"][direction]:                                
                                toCkmer = directionItem["ckmer"]
                                if not toCkmer in connected:
                                    counter_skipped_connections+=1
                                else:
                                    #administration
                                    counter_direct_connections+=1
                                    set_direct_kmers.add(toCkmer)
                                    #create and update to-k-mer
                                    if not toCkmer in self.ckmers.keys():
                                        self._create_ckmer(directionItem)
                                        set_new_kmers.add(toCkmer)
                                    self.ckmers[toCkmer].set_level(level+1)
                                    self._direct_connect_ckmers(self.ckmers[fromCkmer],direction,
                                                                self.ckmers[toCkmer],directionItem["connection"]["direction"],
                                                                directionItem["connection"]["distance"],
                                                                directionItem["connection"]["number"],
                                                                directionItem["connection"]["problem"]>0)
                        else:
                            self._logger.error("unexpected direction "+direction)
            final_candidates = sum(self.ckmers[ckmer].candidate for ckmer in self.ckmers)
            self._logger.debug("skipped %d and found %d connections to %d k-mers: %d new and %d new candidates" % 
                               (counter_skipped_connections, counter_direct_connections,len(set_direct_kmers),
                                len(set_new_kmers),(final_candidates-initial_candidates))) 
                                    
        else:
            self._logger.error("request to "+self.baseUrl+" didn't succeed")
            
    def _create_ckmer(self, item: dict):         
        """
        create k-mer and add automatically base(s)
        """
        ckmer = item["ckmer"]
        if not ckmer in self.ckmers.keys():
            ckmerEntry = self.Ckmer(self, ckmer, item["number"], item["split"])
            for direction in item["rightSplitBase"]:
                base = item["rightSplitBase"][direction]["base"]                
                if not base in self.bases.keys():
                    baseEntry = self.Base(self, base, 
                                            item["rightSplitBase"][direction]["number"], 
                                            item["rightSplitBase"][direction]["ckmers"].keys()) 
                ckmerEntry.add_base(direction,base)
                if not (direction==ckmerEntry.split or ckmerEntry.split=="both"):
                    self._logger.error("unexpected base direction for k-mer "+item["ckmer"])
        else:
            self._logger.warning("k-mer already exists")
        return self.ckmers[ckmer]
            
    def _direct_connect_ckmers(self, ckmer1, direction1: str, 
                               ckmer2, direction2: str, 
                               distance: int, number: int, problem: bool):
        """
        symmetrically set direct connection
        """
        entry = {
            "distance": distance,
            "number": number,
            "problem": problem
        }
        ckmer1.set_direct(direction1,ckmer2.ckmer,entry)
        ckmer2.set_direct(direction2,ckmer1.ckmer,entry)
       
    class Ckmer:
        """internal object representing k-mer in the De Bruijn graph"""
        
        def __init__(self, graph, ckmer: str, number: int, split: str):
            #initialise
            self.graph = graph
            self.ckmer = str(ckmer)
            self.split = str(split)
            self.number = int(number)
            self.level = None
            self.candidate = False
            self.orientation = None
            self.positions = set()
            self.bases = {}
            self.left = {}
            self.right = {}
            self.leftCandidateConnected = False
            self.rightCandidateConnected = False
            #checks
            if self.ckmer in graph.ckmers.keys():
                self.graph._logger.error("creating k-mer that already exists")
            #register
            self.graph.ckmers[ckmer] = self
            
        def __repr__(self):
            info = []
            info.append(str(self.number)+"x")
            if not self.split==None:
                info.append("split "+str(self.split))
            if not self.orientation==None:
                info.append("orientation "+str(self.orientation))
            if not self.level==None:
                info.append("level "+str(self.level))
            if len(self.positions)>0:
                info.append("pos "+(",".join([str(p) for p in self.positions])))
            return "Ckmer("+self.ckmer+" ["+(", ".join(info))+"])" 
        
        def set_orientation(self, orientation: str):
            if orientation in ["forward","backward","both"]:
                if not (self.orientation==None or self.orientation==orientation):
                    self.graph._logger.warning("change orientation of k-mer")
            else:
                self.graph._logger.error("incorrect orientation: "+str(orientation))
                
        def set_direct(self, direction, ckmer, entry):
            if direction=="left":
                self.left[ckmer] = entry
                if self.candidate or self.rightCandidateConnected:
                    self.graph.ckmers[ckmer].set_connected_via_candidate(self.ckmer)
                if not self.level==None:
                    self.graph.ckmers[ckmer].set_level(self.level+1)
            elif direction=="right":
                self.right[ckmer] = entry
                if self.candidate or self.leftCandidateConnected:
                    self.graph.ckmers[ckmer].set_connected_via_candidate(self.ckmer)
                if not self.level==None:
                    self.graph.ckmers[ckmer].set_level(self.level+1)    
            else:
                self.graph._logger.error("incorrect direction: "+str(direction))
                
        def add_positions(self, positions: list):
            self.positions.update(positions)
            
        def add_base(self, direction: str, base: str):
            if direction in ["left","right"]:
                if base in self.graph.bases.keys():
                    if self.ckmer in self.graph.bases[base].ckmers:
                        self.bases[direction] = base
                    else:
                        self.graph._logger.error("k-mer not associated with base "+self.ckmer+" - "+base)
                else:
                    self.graph._logger.error("base not defined in graph")
            else:
                self.graph._logger.error("incorrect direction: "+str(direction))
                
        def set_level(self,level:int):
            #only update if necessary
            if self.level==None or level<self.level:
                if not self.level==None:
                    self.graph._logger.debug("change level from "+str(self.level)+" to "+str(level))
                self.level=int(level)
                #update levels from connected k-mers
                for leftCkmer in self.left.keys():
                    self.graph.ckmers[leftCkmer].set_level(level+1)
                for rightCkmer in self.right.keys():
                    self.graph.ckmers[rightCkmer].set_level(level+1)

        def set_connected_via_candidate(self,ckmer):
            #only update if necessary
            if not self.leftCandidateConnected:
                for leftCkmer in self.left.keys():
                    if leftCkmer==ckmer:
                        self.leftCandidateConnected = True
                        if self.rightCandidateConnected and not self.candidate:
                            self.candidate = True
                        for rightCkmer in self.right.keys():
                            self.graph.ckmers[rightCkmer].set_connected_via_candidate(self.ckmer)
                        break
            #only update if necessary            
            if not self.rightCandidateConnected:
                for rightCkmer in self.right.keys():
                    if rightCkmer==ckmer:
                        self.rightCandidateConnected = True
                        if self.leftCandidateConnected and not self.candidate:
                            self.candidate = True
                        for leftCkmer in self.left.keys():
                            self.graph.ckmers[leftCkmer].set_connected_via_candidate(self.ckmer)
                        break
            
        def set_candidate(self):
            #only update if necessary
            if not self.candidate:
                self.candidate = True
                for ckmer in self.left.keys():
                    self.graph.ckmers[ckmer].set_connected_via_candidate(self.ckmer)
                for ckmer in self.right.keys():
                    self.graph.ckmers[ckmer].set_connected_via_candidate(self.ckmer)
                
    class Base:
        """internal object representing base in the De Bruijn graph"""
        
        def __init__(self, graph, base: str, number: int, ckmers: list):
            #initialise 
            self.graph = graph
            self.base = str(base)
            self.number = int(number)
            self.ckmers = ckmers
            #checks
            if self.base in graph.bases.keys():
                self.graph._logger.error("creating base that already exists")
            #register
            self.graph.bases[base] = self
            
        def __repr__(self):
            info = []
            info.append(str(self.number)+"x")
            info.append(str(len(self.ckmers))+" k-mers")
            return "Base("+self.base+" ["+(", ".join(info))+"])"
               
        
    
    