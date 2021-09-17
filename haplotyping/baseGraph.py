from haplotyping.general import General
from graphviz import Digraph

class Graph:
    """basic or minimal version of the De Bruijn graph"""
    def __init__(self):
        self._name : str = None
        self._k : int = None
        self._start : set = set()
        self._end : set = set()
        self._ckmers : dict = {}
        self._bases : dict = {}
            
    def __repr__(self):
        ck = sum(self._ckmers[ck]._candidate for ck in self._ckmers.keys())
        text = "Graph object with %d candidate k-mers" % (ck)
        return text
        
    def set_name(self, name: str):
        self._name = name
        
        
    def visualize(self,*args, **kwargs):
        #initialise configuration
        config = {
            "hideNoneCandidateBase": True,
            "baseFillColorDefault": "white",
            "baseColorDefault": "grey",
            "nodeFillColorDefault": "white",
            "nodeColorDefault": "grey",
            "nodeFillColorCandidate": "yellow",
            "nodeColorCandidate": "black",
        }
        for key, value in kwargs.items():
            if key in config.keys():
                config[key] = value
        #create graph
        g = Digraph("G")
        g.attr(label="test", labelloc="t", nodesep="0", ranksep="0")
        #define base containers and bases
        baseContainers = {}
        ckmerNodes = {}
        def setBaseContainer(base, container):
            baseContainers[base] = container
            for linkedBase in self._bases[base]._linked:
                if not linkedBase in baseContainers.keys():
                    setBaseContainer(linkedBase, container)
        #try to sort bases
        maxBaseOrder = max([max(self._bases[base]._order["forward"] or 0, 
                                self._bases[base]._order["backward"] or 0) for base in self._bases.keys()])
        sortedBaseList = sorted(self._bases.keys(), 
                                key=lambda base: maxBaseOrder+1 if self._bases[base]._order==None 
                                else self._bases[base]._order)
        for base in sortedBaseList:
            #optionally, hide bases without any candidates
            if config["hideNoneCandidateBase"] and not self._bases[base]._candidate:
                continue  
            #create container for baseContainer (multiple bases can be connected by shared k-mers)
            if not base in baseContainers.keys():
                baseContainer=g.subgraph(name="cluster_container_"+base) 
                with baseContainer as bc:
                    bc.attr(style="invis",nodesep="0", ranksep="0")
                setBaseContainer(base,baseContainer) 
            else:
                baseContainer=baseContainers[base]
                with baseContainer as bc:
                    bc.attr(style="filled", color="lightgrey", fillcolor="whitesmoke", label="")
                setBaseContainer(base,baseContainer) 
            #now create the actual base container
            with baseContainer as bc:      
                #add base container
                with bc.subgraph(name="cluster_"+base) as c: 
                    base_style="filled"
                    base_fillcolor=config["baseFillColorDefault"]
                    if not self._bases[base]._order==None:
                        base_fillcolor="lightgrey"
                    base_color=config["baseColorDefault"]
                    base_title = "<" + "<font point-size=\"6\">"
                    if self._bases[base]._orientation=="forward":
                        base_title += base+"*"
                    elif self._bases[base]._orientation=="backward":
                        base_title += "*"+General.reverse_complement(base)
                    else:
                        base_title += "---orientation unknown---"
                    base_title +="<br/><font color=\"grey\">"
                    base_title +=str(self._bases[base]._number)+"x"
                    if not self._bases[base]._order==None:
                        base_title += " "+str(self._bases[base]._order)
                    base_title +="</font></font>" + ">" 
                    c.attr(style=base_style, color=base_color, 
                           fillcolor=base_fillcolor, label=base_title)
                    #add k-mers to base container
                    for ckmer in self._bases[base]._ckmers.keys():  
                        node_style="filled"
                        node_fillcolor=config["nodeFillColorDefault"]
                        node_color=config["nodeColorDefault"]
                        node_penwidth=1
                        if ckmer in self._ckmers.keys():
                            if self._ckmers[ckmer]._candidate:
                                node_fillcolor=config["nodeFillColorCandidate"]
                                node_color=config["nodeColorCandidate"]
                        rightSplit = (ckmer[:-1]==base)
                        leftSplit = (General.reverse_complement(ckmer[1:])==base)
                        if self._bases[base]._orientation=="forward":
                            if rightSplit and not leftSplit:
                                node_letter=ckmer[-1]
                            elif leftSplit and not rightSplit:
                                node_letter=General.reverse_complement(ckmer[0])
                            else:
                                node_letter="?"
                        elif self._bases[base]._orientation=="backward":
                            if rightSplit and not leftSplit:
                                node_letter=General.reverse_complement(ckmer[-1])
                            elif leftSplit and not rightSplit:
                                node_letter=ckmer[0]
                            else:
                                node_letter="?"
                        else:
                            node_letter="?"
                        node_title = "<" + "<font color=\"blue\">"+node_letter+"</font>" 
                        node_title +="<font point-size=\"8\">"
                        node_title +=str(self._bases[base]._ckmers[ckmer])+"x</font>"
                        node_title += ">"     
                        c.node(base+"_"+ckmer,label=node_title, color=node_color, style=node_style, 
                               fillcolor=node_fillcolor, penwidth=str(node_penwidth))
                        #register
                        if ckmer in ckmerNodes.keys():
                            ckmerNodes[ckmer].append(base+"_"+ckmer)
                        else:
                            ckmerNodes[ckmer]=[base+"_"+ckmer]
        #try to sort k-mers
        maxCkmerOrder = max([self._ckmers[ckmer]._order or 0 for ckmer in self._ckmers.keys()])
        sortedCkmerList = sorted(self._ckmers.keys(), 
                                key=lambda ckmer: maxCkmerOrder+1 if self._ckmers[ckmer]._order==None 
                                else self._ckmers[ckmer]._order)
        
        def get_edge_attributes(ckmer1,direction1,ckmer2):
            #parameters
            orientation1 = self._ckmers[ckmer1]._orientation 
            orientation2 = self._ckmers[ckmer2]._orientation 
            left2 = (ckmer1 in self._ckmers[ckmer2]._left.keys())
            right2 = (ckmer1 in self._ckmers[ckmer2]._right.keys())                
            #default                
            edge_direction = "none"
            edge_color = "grey"
            edge_penwidth = "1"
            edge_style = "solid"
            edge_constraint="false"
            if orientation1=="forward":
                if direction1=="left":
                    edge_direction = "back"
                    edge_color = "black"
                    edge_constraint="true"
                elif direction1=="right":                    
                    edge_direction = "forward"
                    edge_color = "black"
                    edge_constraint="true"
            elif orientation1=="backward":
                if direction1=="left":
                    edge_direction = "forward"
                    edge_color = "black"
                    edge_constraint="true"
                elif direction1=="right":                    
                    edge_direction = "back"
                    edge_color = "black"
                    edge_constraint="true"
            elif orientation2=="forward":
                if left2 and not right2:
                    edge_direction = "forward"
                    edge_color = "black"
                    edge_constraint="true"
                elif right2 and not left2:
                    edge_direction = "back"
                    edge_color = "black"
                    edge_constraint="true"
            elif orientation2=="backward":
                if left2 and not right2:
                    edge_direction = "back"
                    edge_color = "black"
                    edge_constraint="true"
                elif right2 and not left2:
                    edge_direction = "forward"
                    edge_color = "black"
                    edge_constraint="true"
            return (edge_direction,edge_color,edge_penwidth,edge_style,edge_constraint)
        
        def get_edge_label(ckmer1,direction1,ckmer2):
            edge_label = "<"
            edge_label += "<font color=\"grey\">"
            if direction1=="left":
                edge_label += str(self._ckmers[ckmer1]._left[ckmer2]["distance"])
            elif direction1=="right":
                edge_label += str(self._ckmers[ckmer1]._right[ckmer2]["distance"])
            edge_label += "</font><br/>"
            edge_label += "<font point-size=\"8\">"
            if direction1=="left":
                edge_label += str(self._ckmers[ckmer1]._left[ckmer2]["number"])+"x"
            elif direction1=="right":
                edge_label += str(self._ckmers[ckmer1]._right[ckmer2]["number"])+"x"
            edge_label += "</font>" 
            edge_label += ">"
            return edge_label
        
        #remember edges
        edges = set()
        for ckmer in sortedCkmerList:
            if ckmer in ckmerNodes.keys():
                for connectedCkmer in self._ckmers[ckmer]._left.keys():
                    if connectedCkmer in ckmerNodes.keys():
                        #only pass once
                        edgeKey = (ckmer+"_"+connectedCkmer) if ckmer<connectedCkmer else (connectedCkmer+"_"+ckmer)
                        if edgeKey in edges:
                            continue
                        else:
                            edges.add(edgeKey)
                        #attributes
                        edge_label = get_edge_label(ckmer,"left",connectedCkmer)
                        (edge_direction,edge_color,edge_penwidth,edge_style,edge_constraint) = get_edge_attributes(
                            ckmer,"left",connectedCkmer)
                        for node in ckmerNodes[ckmer]:
                            for connectedNode in ckmerNodes[connectedCkmer]:
                                g.edge(node,connectedNode, label=edge_label, constraint=edge_constraint,
                                       dir=edge_direction, color=edge_color, penwidth=edge_penwidth)  
                for connectedCkmer in self._ckmers[ckmer]._right.keys():
                    if connectedCkmer in ckmerNodes.keys():
                        #only pass once
                        edgeKey = (ckmer+"_"+connectedCkmer) if ckmer<connectedCkmer else (connectedCkmer+"_"+ckmer)
                        if edgeKey in edges:
                            continue
                        else:
                            edges.add(edgeKey)
                        #attributes
                        edge_label = get_edge_label(ckmer,"right",connectedCkmer)
                        (edge_direction,edge_color,edge_penwidth,edge_style,edge_constraint) = get_edge_attributes(
                            ckmer,"right",connectedCkmer)
                        node = sorted(ckmerNodes[ckmer])[0]
                        connectedNode = sorted(ckmerNodes[connectedCkmer])[0]
                        g.edge(node,connectedNode, label=edge_label, constraint=edge_constraint,
                           dir=edge_direction, color=edge_color, penwidth=edge_penwidth)  
        return g
        
            
    class Ckmer:
        """internal object representing splitting k-mers in the De Bruijn graph"""

        def __init__(self, graph, ckmer: str, orientation: str, number: int, split: str):
            """
            - ckmer: canonical representation of the splitting k-mer
            - orientation: orientation of the k-mer (forward/backward)
            - number: frequency
            - split: location of the split (left/right/both)
            """
            #checks
            assert len(ckmer)==graph._k
            assert orientation in ["forward","backward"]
            assert split in ["left","right","both"]
            assert number>=0
            #initialise
            self._graph = graph
            self._ckmer : str = str(ckmer)
            self._orientation : str = str(orientation)
            self._split : str = str(split)
            self._number : int = int(number)
            #other variables
            self._key = (self._ckmer,self._orientation)
            self._candidate : bool = False
            self._incoming : dict = {}
            self._outgoing : dict = {}
            self._bases : dict = {}
            self._order : int = None
            #check for duplicates
            if self._key in self._graph._ckmers.keys():
                raise Exception("creating k-mer with orientation that already exists")
            #link bases
            if (self._split == "both") or (self._split == "right"):
                base = self._ckmer[:-1]
                if base in self._graph._bases.keys():
                    assert self._ckmer in self._graph._bases[base]._ckmers
                    self._bases["right"] = base
                    self._graph._bases[base]._set_orientation(self._ckmer,"right",self._orientation)
                    if ("left" in self._bases.keys()) and (not self._bases["left"]==base):
                        self._graph._bases[self._bases["right"]]._linked[orientation].add(self._bases["left"])
                        if orientation=="forward":
                            self._graph._bases[self._bases["left"]]._linked["backward"].add(self._bases["right"])
                        else:
                            self._graph._bases[self._bases["left"]]._linked["forward"].add(self._bases["right"])
            if (self._split == "both") or (self._split == "left"):
                base = General.reverse_complement(self._ckmer[1:])
                if base in self._graph._bases.keys():
                    assert self._ckmer in self._graph._bases[base]._ckmers
                    self._bases["left"] = base
                    self._graph._bases[base]._set_orientation(self._ckmer,"left",self._orientation)
                    if ("right" in self._bases.keys()) and (not self._bases["right"]==base):
                        if orientation=="forward":
                            self._graph._bases[self._bases["right"]]._linked["backward"].add(self._bases["left"])
                        else:
                            self._graph._bases[self._bases["right"]]._linked["forward"].add(self._bases["left"])
                        self._graph._bases[self._bases["left"]]._linked[orientation].add(self._bases["right"])
                
            #register
            self._graph._ckmers[self._key] = self
            
        def __repr__(self):
            info = []
            info.append(str(self._number)+"x")
            info.append("split "+str(self._split))
            if self._candidate:
                info.append("candidate")
            return "Ckmer("+self._ckmer+" "+self._orientation+"["+(", ".join(info))+"])" 
        
        def set_incoming(self, ckmerKey: str, distance: int, number: int, problem: bool):
            assert ckmerKey in self._graph._ckmers.keys()
            assert distance>0
            assert number>=0
            self._incoming[ckmerKey] = {"distance": distance, "number": number, "problem": problem} 
            
        def set_outgoing(self, ckmerKey: str, distance: int, number: int, problem: bool):
            assert ckmerKey in self._graph._ckmers.keys()
            assert distance>0
            assert number>=0
            self._outgoing[ckmerKey] = {"distance": distance, "number": number, "problem": problem} 
            
        def set_candidate(self):
            self._candidate = True
            for direction in self._bases.keys():
                base = self._bases[direction]
                self._graph._bases[base].set_candidate()
                
        def set_order(self, order: int):
            self._order = order
            for direction in self._bases.keys():
                base = self._bases[direction]
                self._graph._bases[base].set_order()
            
    class Base:
        """
        internal object representing base in the De Bruijn graph
        """
        
        def __init__(self, graph, base: str, number: int, ckmers: dict):
            assert len(base)==graph._k - 1
            assert number>=0 
            assert len(ckmers)>0
            #initialise 
            self._graph = graph
            self._base : str = str(base)
            self._number : int = int(number)
            self._ckmers : dict = ckmers.copy()
            #other variables
            self._candidate : dict = {"forward": False, "backward": False}
            self._linked : dict = {"forward": set(), "backward": set()}
            self._orientation : dict = {"forward": set(), "backward": set()}
            self._order : dict = {"forward": None, "backward": None}
            #check
            if self._base in self._graph._bases.keys():
                raise Exception("creating base that already exists")
            #link k-mers
            for orientation in ["forward","backward"]:
                for ckmer in self._ckmers:
                    ckmerKey = (ckmer,orientation)
                    #only if k-mer with orientation exists
                    if ckmerKey in self._graph._ckmers.keys():
                        #right
                        if self._base==ckmer[:-1]:
                            self._graph._ckmers[ckmerKey]._bases["right"] = self._base
                            self._set_orientation(ckmer,"right",orientation)
                            if "left" in self._graph._ckmers[ckmerKey]._bases.keys():
                                linkedBase=self._graph._ckmers[ckmerKey]._bases["left"]
                                if not self._base==linkedBase:
                                    self._linked[orientation].add(linkedBase)
                                    if orientation=="forward":
                                        self._graph._bases[linkedBase]._linked["backward"].add(self._base)
                                    else:
                                        self._graph._bases[linkedBase]._linked["forward"].add(self._base)
                        #left
                        if self._base==General.reverse_complement(ckmer[1:]):
                            self._graph._ckmers[ckmerKey]._bases["left"] = self._base
                            self._set_orientation(ckmer,"left",orientation)
                            if "right" in self._graph._ckmers[ckmerKey]._bases.keys():
                                linkedBase=self._graph._ckmers[ckmerKey]._bases["right"]
                                if not self._base==linkedBase:
                                    if orientation=="forward":
                                        self._linked["backward"].add(linkedBase)
                                    else:
                                        self._linked["forward"].add(linkedBase)
                                    self._graph._bases[linkedBase]._linked[orientation].add(self._base)
            self.set_candidate()
            #register
            self._graph._bases[self._base] = self
            
        def __repr__(self):
            info = []
            info.append(str(self._number)+"x")
            info.append(str(len(self._ckmers))+" k-mers")
            if self._candidate:
                info.append("candidate")
            return "Base("+self._base+" ["+(", ".join(info))+"])"
        
        def _set_orientation(self, ckmer: str, split: str, orientation: str):
            assert len(ckmer)==self._graph._k
            assert ckmer in self._ckmers.keys()
            assert split in ["right","left"]
            assert orientation in ["forward","backward"]
            if split=="right":
                self._orientation[orientation].add(ckmer)
            else:
                if orientation=="forward":
                    self._orientation["backward"].add(ckmer)
                else:
                    self._orientation["forward"].add(ckmer)
              
        def set_candidate(self):
            for orientation in ["forward","backward"]:
                self._candidate[orientation]=False
                for ckmer in self._ckmers:
                    #only if k-mer exists
                    ckmerKey = (ckmer,orientation)
                    if ckmerKey in self._graph._ckmers.keys():
                        if self._graph._ckmers[ckmerKey]._candidate:
                            self._candidate[orientation]=True
            
        def set_order(self):
            #try to inherit order from k-mers
            for orientation in ["forward","backward"]:
                self._order[orientation]=None
                for ckmer in self._ckmers:
                    #only if k-mer exists
                    ckmerKey = (ckmer,orientation)
                    if ckmerKey in self._graph._ckmers.keys():
                        if not self._graph._ckmers[ckmerKey]._order==None:
                            if self._order[orientation]==None:
                                self._order[orientation]=self._graph._ckmers[ckmerKey]._order
                            else:
                                self._order[orientation]=min(self._graph._ckmers[ckmerKey]._order,
                                                             self._order[orientation])
                            
