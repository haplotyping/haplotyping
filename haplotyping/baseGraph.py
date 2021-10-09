from haplotyping.general import General
from graphviz import Digraph
import pandas as pd
import networkit as nk
import numpy as np

class Graph:
    """basic or minimal version of the De Bruijn graph"""
    def __init__(self):
        self._name : str = None
        self._k : int = None
        self._start : set = set()
        self._end : set = set()
        self._ckmers : dict = {}
        self._orientatedCkmers : dict = {}
        self._orientatedBases : dict = {}
        self._missingConnections : dict = {}
        self._selected : set = set()
        self._distances = None
            
    def __repr__(self):
        ck = sum(self._orientatedCkmers[ck]._candidate for ck in self._orientatedCkmers.keys())
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
            "nodePenWidthDefault": 1,
            "nodeColorSelected": "blue",
            "nodePenWidthSelected": 3,
            "nodeFillColorCandidate": "yellow",
            "nodeColorCandidate": "black",
            "edgeColorDefault": "grey",
            "edgePenWidthDefault": 1,
            "edgeStyleDefault": "solid",
            "edgeColorAuto": "grey",
            "edgePenWidthAuto": 1,
            "edgeStyleAuto": "dashed",
        }
        for key, value in kwargs.items():
            if key in config.keys():
                config[key] = value
        
        #create graph
        g = Digraph("G")
        g.attr(label="test", labelloc="t", nodesep="0", ranksep="0")
        
        #define base containers and bases
        orientatedBaseContainers = {}
        orientatedCkmerNodes = {}
        
        def setOrientatedBaseContainer(orientatedBase, container):
            orientatedBaseContainers[orientatedBase] = container
            for linkedOrientatedBase in self._orientatedBases[orientatedBase]._linked:
                if not linkedOrientatedBase in orientatedBaseContainers.keys():
                    setOrientatedBaseContainer(linkedOrientatedBase, container)
        
        #try to sort orientated bases
        maxOrientatedBaseOrder = max([self._orientatedBases[orientatedBase]._order or 0 
                                      for orientatedBase in self._orientatedBases.keys()])
        sortedOrientatedBaseList = sorted(self._orientatedBases.keys(), 
                                key=lambda orientatedBase: maxOrientatedBaseOrder+1 
                                          if self._orientatedBases[orientatedBase]._order==None 
                                          else self._orientatedBases[orientatedBase]._order)
        
        #loop over sorted orientated bases
        for orientatedBase in sortedOrientatedBaseList:
            
            #optionally, hide bases without any candidates
            if config["hideNoneCandidateBase"] and not self._orientatedBases[orientatedBase]._candidate:
                continue  
            
            #create container for orientatedBaseContainer (multiple bases can be connected by shared k-mers)
            #and make container visible if it contains multiple bases
            if not orientatedBase in orientatedBaseContainers.keys():
                orientatedBaseContainer=g.subgraph(name="cluster_container_"+
                                                   str(orientatedBase[0])+"_"+str(orientatedBase[1])) 
                with orientatedBaseContainer as obc:
                    obc.attr(style="invis",nodesep="0", ranksep="0")
                setOrientatedBaseContainer(orientatedBase,orientatedBaseContainer) 
            else:
                orientatedBaseContainer=orientatedBaseContainers[orientatedBase]
                with orientatedBaseContainer as obc:
                    obc.attr(style="filled", color="lightgrey", fillcolor="whitesmoke", label="")
                setOrientatedBaseContainer(orientatedBase,orientatedBaseContainer) 
                
            #now create the actual orientated base container in the container
            with orientatedBaseContainer as obc:      
                
                #add orientated base container
                with obc.subgraph(name="cluster_"+str(orientatedBase[0])+"_"+str(orientatedBase[1])) as c: 
                    #define layout
                    base_style="filled"
                    base_fillcolor=config["baseFillColorDefault"]
                    if not self._orientatedBases[orientatedBase]._order==None:
                        base_fillcolor="lightgrey"
                    base_color=config["baseColorDefault"]
                    base_title = "<" + "<font point-size=\"6\">"
                    if orientatedBase[1]=="forward":
                        base_title += orientatedBase[0]+"*"
                    elif orientatedBase[1]=="backward":
                        base_title += "*"+General.reverse_complement(orientatedBase[0])
                    else:
                        base_title += "---orientation unknown---"
                    base_title +="<br/><font color=\"grey\">"
                    base_title +=str(self._orientatedBases[orientatedBase]._number)+"x"
                    if not self._orientatedBases[orientatedBase]._order==None:
                        base_title += " "+str(self._orientatedBases[orientatedBase]._order)
                    base_title +="</font></font>" + ">" 
                    c.attr(style=base_style, color=base_color, 
                           fillcolor=base_fillcolor, label=base_title)
                    #add k-mers to base container
                    for orientatedCkmer in self._orientatedBases[orientatedBase]._orientatedCkmers:  
                        #define layout
                        node_style="filled"
                        node_fillcolor=config["nodeFillColorDefault"]
                        node_color=config["nodeColorDefault"]
                        node_penwidth=config["nodePenWidthDefault"]
                        if orientatedCkmer in self._orientatedCkmers.keys():
                            if self._orientatedCkmers[orientatedCkmer]._candidate:
                                node_fillcolor=config["nodeFillColorCandidate"]
                                node_color=config["nodeColorCandidate"]
                        if orientatedCkmer in self._selected:
                            node_color=config["nodeColorSelected"]
                            node_penwidth=config["nodePenWidthSelected"]
                        node_letter = "?"
                        if orientatedBase[1]=="forward":
                            if orientatedCkmer[1]=="forward":
                                assert orientatedBase[0]==orientatedCkmer[0][:-1]
                                node_letter = orientatedCkmer[0][-1]
                            elif orientatedCkmer[1]=="backward":
                                assert orientatedBase[0]==General.reverse_complement(orientatedCkmer[0][1:])
                                node_letter = General.reverse_complement(orientatedCkmer[0][0])
                        elif orientatedBase[1]=="backward":
                            if orientatedCkmer[1]=="forward":
                                assert General.reverse_complement(orientatedBase[0])==orientatedCkmer[0][1:]
                                node_letter = orientatedCkmer[0][0]
                            elif orientatedCkmer[1]=="backward":
                                assert orientatedBase[0]==orientatedCkmer[0][:-1]
                                node_letter = General.reverse_complement(orientatedCkmer[0][-1])
                        node_title = "<" + "<font color=\"blue\">"+node_letter+"</font>" 
                        node_title +="<font point-size=\"8\">"
                        node_title +=str(self._orientatedCkmers[orientatedCkmer]._number)+"x</font>"
                        node_title += ">"
                        node_key = (str(orientatedBase[0])+"_"+str(orientatedBase[1])+"_"+
                                    str(orientatedCkmer[0])+"_"+str(orientatedCkmer[1]))
                        c.node(node_key,label=node_title, color=node_color, style=node_style, 
                               fillcolor=node_fillcolor, penwidth=str(node_penwidth))
                        #register
                        if orientatedCkmer in orientatedCkmerNodes.keys():
                            assert orientatedBase[1] not in orientatedCkmerNodes[orientatedCkmer]
                            orientatedCkmerNodes[orientatedCkmer][orientatedBase[1]] = node_key
                            for other_node_key in orientatedCkmerNodes[orientatedCkmer].values():
                                if not node_key==other_node_key:
                                    g.edge(other_node_key, node_key, constraint="false",
                                       dir="none", color="grey", style="dashed", penwidth="2")  
                                    
                        else:
                            orientatedCkmerNodes[orientatedCkmer] = {orientatedBase[1]: node_key}
                            
        #try to sort orientated k-mers
        maxOrientatedCkmerOrder = max([self._orientatedCkmers[orientatedCkmer]._order or 0 
                                       for orientatedCkmer in self._orientatedCkmers.keys()])
        sortedOrientatedCkmerList = sorted(self._orientatedCkmers.keys(), 
                                key=lambda orientatedCkmer: maxOrientatedCkmerOrder+1 
                                           if self._orientatedCkmers[orientatedCkmer]._order==None 
                                           else self._orientatedCkmers[orientatedCkmer]._order)
        
        #edge label definition
        def get_edge_label(distance: int, number: int):
            edge_label = "<"
            edge_label += "<font color=\"grey\">"+str(distance)+"</font><br/>"
            edge_label += "<font point-size=\"8\">"+str(number)+"</font>" 
            edge_label += ">"
            return edge_label
        
        #remember edges
        edges = set()

        #loop over the sorted orientated k-mers
        for orientatedCkmer in sortedOrientatedCkmerList:
            if orientatedCkmer in orientatedCkmerNodes.keys():                
                kmerSet = self._orientatedCkmers[orientatedCkmer]._incoming.keys()
                for connectedCkmer in kmerSet:
                    if connectedCkmer in orientatedCkmerNodes.keys():
                        edgeKey = connectedCkmer + orientatedCkmer
                        if edgeKey in edges:
                            #only pass once
                            continue
                        else:
                            edges.add(edgeKey)
                        edge_label = get_edge_label(
                            self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["distance"],
                            self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["number"]
                        )
                        edge_direction = "forward"
                        edge_constraint="true"
                        edge_color = config["edgeColorDefault"]
                        edge_penwidth = config["edgePenWidthDefault"]
                        edge_style = config["edgeStyleDefault"]
                        #detect auto generated edges
                        if self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["number"]==0:
                            edge_color = config["edgeColorAuto"]
                            edge_penwidth = config["edgePenWidthAuto"]
                            edge_style = config["edgeStyleAuto"]
                        #choose preferred nodes for the edges
                        if "backward" in orientatedCkmerNodes[orientatedCkmer]:
                            node = orientatedCkmerNodes[orientatedCkmer]["backward"]
                        else:
                            node = orientatedCkmerNodes[orientatedCkmer]["forward"]
                        if "forward" in orientatedCkmerNodes[connectedCkmer]:
                            connectedNode = orientatedCkmerNodes[connectedCkmer]["forward"]
                        else:
                            connectedNode = orientatedCkmerNodes[connectedCkmer]["backward"]
                        #create the edge
                        g.edge(connectedNode, node, label=edge_label, constraint=edge_constraint, style=edge_style, 
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth))  
                        
        #draw missing connections
        for oc1 in self._missingConnections:
            if oc1 in orientatedCkmerNodes:
                for oc2 in self._missingConnections[oc1]:
                    if oc2 in orientatedCkmerNodes:
                        for d1 in orientatedCkmerNodes[oc1]:
                            n1 = orientatedCkmerNodes[oc1][d1]
                            for d2 in orientatedCkmerNodes[oc2]:
                                n2 = orientatedCkmerNodes[oc2][d2]
                                g.edge(n1,n2, constraint="false",
                                   dir="none", color="grey", penwidth="1", style="dashed")  
                    
        return g
        
    def getDirectDistances(self):
        if self._distances == None:
            self._computeDistances()
        return self._distances["directDistance"]
        
    def getDirectNumbers(self):
        if self._distances == None:
            self._computeDistances()
        return self._distances["directNumber"]    
    
    def getShortestDistances(self):
        if self._distances == None:
            self._computeDistances()
        return self._distances["shortest"]
        
    def getShortestDistance(self, fromOrientatedKmers: set, toOrientatedKmers: set):
        if self._distances == None:
            self._computeDistances()
        distances = []
        for c1 in fromOrientatedKmers:
            for c2 in toOrientatedKmers:
                if self._distances["shortest"].at[c1,c2]!=0:
                    distances.append(self._distances["shortest"].at[c1,c2])
        if len(distances)==0:
            return 0
        else:
            return(sorted(distances)[0])
        
    def _resetDistances(self):
        if not self._distances==None:
            self._logger.debug("reset distances for the De Bruijn Graph")  
            self._distances = None
        
    def _computeDistances(self):
        #initialise
        orientatedCkmerList = [c for c in self._orientatedCkmers]
        self._distances = {
            "directDistance": pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0),
            "directNumber": pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0),
            "shortest": pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0)
        }
        self._logger.debug("compute distances for the De Bruijn Graph")  
        #direct
        for orientatedCkmer in orientatedCkmerList:
            for connectedCkmer in self._orientatedCkmers[orientatedCkmer]._incoming:
                if connectedCkmer in orientatedCkmerList:
                    d = self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["distance"]
                    n = self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["number"]
                    self._distances["directDistance"].at[orientatedCkmer,connectedCkmer] = -d
                    self._distances["directDistance"].at[connectedCkmer,orientatedCkmer] = d
                    self._distances["directNumber"].at[orientatedCkmer,connectedCkmer] = n
                    self._distances["directNumber"].at[connectedCkmer,orientatedCkmer] = n
        self._logger.debug("found 2 x {:d} graph direct connected pairs".
              format(int(np.count_nonzero(self._distances["directDistance"]>0))))   
        #compute distances between all nkmers
        nkg = nk.graph.Graph(n=len(orientatedCkmerList), weighted=True, directed=True)
        #add relations to graph
        for i in range(len(orientatedCkmerList)):
            for j in range(len(orientatedCkmerList)):
                if self._distances["directDistance"].iloc[i,j]>0:
                    nkg.addEdge(i,j,self._distances["directDistance"].iloc[i,j])
        nkd=nk.distance.APSP(nkg)
        nkd.run()
        distance = nkd.getDistances()
        #remove maxint and add negative distances
        for i in range(len(distance)):
            for j in range(len(distance[i])):
                if distance[i][j]>10**100:
                    distance[i][j]=0
                elif distance[i][j]>0:
                    distance[j][i]=-1*distance[i][j]
        self._distances["shortest"] = pd.DataFrame(data=distance, index=orientatedCkmerList, 
                                                   columns=orientatedCkmerList).astype(int)
        self._logger.debug("found 2 x {:d} = {:d} graph disconnected pairs".
              format(int((np.count_nonzero(self._distances["shortest"]==0)-len(self._distances["shortest"].index))/2),
                     np.count_nonzero(self._distances["shortest"]==0)-len(self._distances["shortest"].index)))
        
            
    class OrientatedCkmer:
        """internal object representing splitting orientated k-mers in the De Bruijn graph"""

        def __init__(self, graph, ckmer: str, orientation: str, number: int, split: str):
            """
            - ckmer: canonical representation of the splitting k-mer
            - orientation: orientation from the k-mer (forward/backward)
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
            self._orientatedBases : dict = {}
            self._order : int = None
            #check for duplicates
            if self._key in self._graph._orientatedCkmers.keys():
                raise Exception("creating k-mer with orientation that already exists")
            #register
            self._graph._orientatedCkmers[self._key] = self 
            if self._ckmer in self._graph._ckmers.keys():
                self._graph._ckmers[self._ckmer].add(self._key)
            else:
                self._graph._ckmers[self._ckmer] = set([self._key])
            #link bases
            if (self._split == "both") or (self._split == "right"):
                base = self._ckmer[:-1]
                baseOrientation = self._orientation
                orientatedBaseKey = (base,baseOrientation)
                self._orientatedBases["right"] = orientatedBaseKey
                if not orientatedBaseKey in self._graph._orientatedBases.keys():
                    orientatedBaseEntry = self._graph.OrientatedBase(self._graph, base, baseOrientation)
                self._graph._orientatedBases[orientatedBaseKey].add_ckmer(self._ckmer, self._orientation)      
            if (self._split == "both") or (self._split == "left"):
                base = General.reverse_complement(self._ckmer[1:])
                baseOrientation = "forward" if self._orientation=="backward" else "backward"
                orientatedBaseKey = (base,baseOrientation)
                self._orientatedBases["left"] = orientatedBaseKey
                if not orientatedBaseKey in self._graph._orientatedBases.keys():
                    orientatedBaseEntry = self._graph.OrientatedBase(self._graph, base, baseOrientation)
                self._graph._orientatedBases[orientatedBaseKey].add_ckmer(self._ckmer, self._orientation)     
            
        def __repr__(self):
            info = []
            info.append(str(self._number)+"x")
            info.append("split "+str(self._split))
            if self._candidate:
                info.append("candidate")
            if self._order:
                info.append("order: "+str(self._order))
            return "OrientatedCkmer("+self._ckmer+" "+self._orientation+"["+(", ".join(info))+"])" 
        
        def set_incoming(self, orientatedCkmerKey: str, distance: int, number: int, problem: bool):
            assert orientatedCkmerKey in self._graph._orientatedCkmers.keys()
            assert distance>0
            assert number>=0
            self._incoming[orientatedCkmerKey] = {"distance": distance, "number": number, "problem": problem} 
            
        def set_outgoing(self, orientatedCkmerKey: str, distance: int, number: int, problem: bool):
            assert orientatedCkmerKey in self._graph._orientatedCkmers.keys()
            assert distance>0
            assert number>=0
            self._outgoing[orientatedCkmerKey] = {"distance": distance, "number": number, "problem": problem} 
            
        def set_candidate(self):
            self._candidate = True
            for direction in self._orientatedBases.keys():
                orientatedBaseKey = self._orientatedBases[direction]
                self._graph._orientatedBases[orientatedBaseKey].set_candidate()
                
        def set_order(self, order: int):
            self._order = order
            for direction in self._orientatedBases.keys():
                orientatedBaseKey = self._orientatedBases[direction]
                self._graph._orientatedBases[orientatedBaseKey].set_order()
            
    class OrientatedBase:
        """
        internal object representing orientated base in the De Bruijn graph
        """
        
        def __init__(self, graph, base: str, orientation: str):
            """
            - base: base for right splitting k-mer
            - orientation: orientation from the right splitting k-mer base (forward/backward)
            """
            #checks
            assert len(base)==graph._k - 1
            assert orientation in ["forward","backward"]
            #initialise 
            self._graph = graph
            self._base : str = str(base)
            self._orientation : str = str(orientation)            
            #other variables
            self._key = (self._base,self._orientation)
            self._candidate : bool = False
            self._order : int = None
            self._number : int = 0
            self._orientatedCkmers : set = set()
            self._linked : set = set()
            #check
            if self._key in self._graph._orientatedBases.keys():
                raise Exception("creating base that already exists")
            #register
            self._graph._orientatedBases[self._key] = self
            
        def __repr__(self):
            info = []
            info.append(str(self._number)+"x")
            info.append(str(len(self._orientatedCkmers))+" k-mers")
            if self._candidate:
                info.append("candidate")
            if self._order:
                info.append("order: "+str(self._order))
            return "OrientatedBase("+self._base+" "+self._orientation+"["+(", ".join(info))+"])"
        
        def set_candidate(self):
            self._candidate = True
            
        def set_order(self):
            self._order = None
            for orientatedCkmerKey in self._orientatedCkmers:
                if not self._graph._orientatedCkmers[orientatedCkmerKey]._order == None:
                    if self._order == None:
                        self._order = self._graph._orientatedCkmers[orientatedCkmerKey]._order
                    else:
                        self._order = min(self._order, self._graph._orientatedCkmers[orientatedCkmerKey]._order)
                                           
        
        def add_ckmer(self, ckmer, orientation):
            orientatedCkmerKey = (ckmer, orientation)            
            if orientatedCkmerKey in self._graph._orientatedCkmers.keys():
                assert orientatedCkmerKey not in self._orientatedCkmers
                self._orientatedCkmers.add(orientatedCkmerKey)
                self._number += self._graph._orientatedCkmers[orientatedCkmerKey]._number
                #check candidate
                self._candidate = (self._candidate or self._graph._orientatedCkmers[orientatedCkmerKey]._candidate)
                #check order
                if not self._graph._orientatedCkmers[orientatedCkmerKey]._order == None:
                    if self._order == None:
                        self._order = self._graph._orientatedCkmers[orientatedCkmerKey]._order
                    else:
                        self._order = min(self._order, self._graph._orientatedCkmers[orientatedCkmerKey]._order)
                #check linked
                for baseKey in self._graph._orientatedCkmers[orientatedCkmerKey]._orientatedBases.values():
                    if not (baseKey == self._key):
                        self._linked.add(baseKey)
                        self._graph._orientatedBases[baseKey]._linked.add(self._key)
            else:
                raise Exception("can't add non-existing orientated canonical k-mer")
                            
