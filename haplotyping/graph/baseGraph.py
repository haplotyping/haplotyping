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
        
        
    def visualize(self, *args, **kwargs):
        
        #initialise configuration
        config = {
            "showAllBases": False,
            "showPrePostBases": True,
            "showPrePostBasesMaxSteps": 2,
            "showDeadEnds": True,
            "baseStyleDefault": "filled",
            "baseFillColorDefault": "white",
            "basePenWidthDefault": 1,
            "baseColorDefault": "grey",
            "baseStyleOrder": "filled",
            "baseFillColorOrder": "grey95",
            "basePenWidthOrder": 1.5,
            "baseColorOrder": "black",
            "baseStyleCandidate": "filled",
            "baseFillColorCandidate": "grey85",
            "basePenWidthCandidate": 1.5,
            "baseColorCandidate": "black",
            "baseStylePrePost": "filled,dashed",
            "baseFillColorPrePost": "white",
            "basePenWidthPrePost": 1,
            "baseColorPrePost": "grey",
            "nodeFillColorDefault": "white",
            "nodeColorDefault": "grey",
            "nodePenWidthDefault": 1,
            "nodeColorSelected": "blue",
            "nodePenWidthSelected": 3,
            "nodeFillColorCandidate": "lightyellow1",
            "nodeColorCandidate": "black",
            "nodeFillColorStartEnd": "yellow",
            "nodeColorStartEnd": "black",
            "edgeColorDefault": "grey",
            "edgePenWidthDefault": 1,
            "edgeStyleDefault": "solid",
            "edgeColorAuto": "grey",
            "edgePenWidthAuto": 1,
            "edgeStyleAuto": "dashed",
            "edgeColorMissingPath": "red", 
            "edgeStyleMissingPath": "dashed",
        }
        for key, value in kwargs.items():
            if key in config.keys():
                config[key] = value
        
        #create graph
        g = Digraph("G")
        graph_label = "Graph" if self._name==None else str(self._name)
        g.attr(label=graph_label, labelloc="t", nodesep="0", ranksep="0")
        
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
        test = False
        for orientatedBase in sortedOrientatedBaseList:
            
            #optionally, hide bases
            if not config["showAllBases"]:
                if (config["showPrePostBases"] and not (self._orientatedBases[orientatedBase]._preStart==None and 
                                                        self._orientatedBases[orientatedBase]._postEnd==None)):
                    if config["showPrePostBasesMaxSteps"] and isinstance(config["showPrePostBasesMaxSteps"],int):
                        if not (self._orientatedBases[orientatedBase]._preStart==None or 
                                self._orientatedBases[orientatedBase]._preStart<=config["showPrePostBasesMaxSteps"]):
                            continue
                        elif not (self._orientatedBases[orientatedBase]._postEnd==None or 
                                self._orientatedBases[orientatedBase]._postEnd<=config["showPrePostBasesMaxSteps"]):
                            continue
                elif self._orientatedBases[orientatedBase]._candidate:
                    pass
                else:
                    continue
            
            #create container for orientatedBaseContainer (multiple bases can be connected by shared k-mers)
            #and make container visible if it contains multiple bases
            if not orientatedBase in orientatedBaseContainers.keys():
                orientatedBaseContainerName = "cluster_container_{}_{}".format(orientatedBase[0],orientatedBase[1])
                orientatedBaseContainer=g.subgraph(name=orientatedBaseContainerName)
                with orientatedBaseContainer as obc:
                    obc.attr(style="invis",nodesep="0", ranksep="0")
                setOrientatedBaseContainer(orientatedBase,orientatedBaseContainerName) 
            else:
                orientatedBaseContainerName=orientatedBaseContainers[orientatedBase]
                orientatedBaseContainer=g.subgraph(name=orientatedBaseContainerName)
                with orientatedBaseContainer as obc:
                    obc.attr(style="filled", color="lightgrey", fillcolor="whitesmoke", label="")
                setOrientatedBaseContainer(orientatedBase,orientatedBaseContainerName) 
              
            #now create the actual orientated base container in the container
            orientatedBaseContainer=g.subgraph(name=orientatedBaseContainerName)  
            with orientatedBaseContainer as obc:      
                
                #add orientated base container
                with obc.subgraph(name="cluster_"+str(orientatedBase[0])+"_"+str(orientatedBase[1])) as c: 
                    #define layout
                    base_style=config["baseStyleDefault"]
                    base_fillcolor=config["baseFillColorDefault"]
                    base_color=config["baseColorDefault"]
                    base_penwidth=config["basePenWidthDefault"]
                    if self._orientatedBases[orientatedBase]._candidate:
                        base_style=config["baseStyleCandidate"]
                        base_fillcolor=config["baseFillColorCandidate"]
                        base_color=config["baseColorCandidate"]
                        base_penwidth=config["basePenWidthCandidate"]
                    elif not (self._orientatedBases[orientatedBase]._preStart==None and 
                          self._orientatedBases[orientatedBase]._postEnd==None):
                        base_style=config["baseStylePrePost"]
                        base_fillcolor=config["baseFillColorPrePost"]
                        base_color=config["baseColorPrePost"]
                        base_penwidth=config["basePenWidthPrePost"]
                    else:
                        base_style=config["baseStyleOrder"]
                        base_fillcolor=config["baseFillColorOrder"]
                        base_color=config["baseColorOrder"]
                        base_penwidth=config["basePenWidthOrder"]
                    base_label = self._visualize_base_label(orientatedBase) 
                    c.attr(style=base_style, color=base_color, 
                           fillcolor=base_fillcolor, label=base_label, penwidth=str(base_penwidth))
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
                        if orientatedCkmer in self._start or orientatedCkmer in self._end:
                            node_fillcolor=config["nodeFillColorStartEnd"]
                            node_color=config["nodeColorStartEnd"]
                        if orientatedCkmer in self._selected:
                            node_color=config["nodeColorSelected"]
                            node_penwidth=config["nodePenWidthSelected"]
                        node_label = self._visualize_node_label(orientatedBase, orientatedCkmer)
                        node_key = (str(orientatedBase[0])+"_"+str(orientatedBase[1])+"_"+
                                    str(orientatedCkmer[0])+"_"+str(orientatedCkmer[1]))
                        c.node(node_key,label=node_label, color=node_color, style=node_style, 
                               fillcolor=node_fillcolor, penwidth=str(node_penwidth))
                        if ((config["showAllBases"] or config["showDeadEnds"]) 
                            and not self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd==None):
                            dead_key = (str(orientatedBase[0])+"_"+str(orientatedBase[1])+"_"+
                                        str(self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[0]))
                            g.node(dead_key, shape="point")
                            g.edge(dead_key, node_key, constraint="true",
                              label=self._visualize_dead_end_label(
                                  self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[1]),
                                   dir="forward", color=config["edgeColorAuto"], weight="1",
                                   style=config["edgeStyleAuto"], penwidth=str(config["edgePenWidthAuto"]))
                        if ((config["showAllBases"] or config["showDeadEnds"])
                            and not self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd==None):
                            dead_key = (str(orientatedBase[0])+"_"+str(orientatedBase[1])+"_"+
                                        str(self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[0]))
                            g.node(dead_key, shape="point")
                            g.edge(node_key, dead_key, constraint="true", 
                              label=self._visualize_dead_end_label(
                                  self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[1]),
                                   dir="forward", color=config["edgeColorAuto"], weight="1",
                                   style=config["edgeStyleAuto"], penwidth=str(config["edgePenWidthAuto"]))
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
                        edge_label = self._visualize_edge_label(
                            self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]
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
                        if self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["path"]==None:
                            edge_color = config["edgeColorMissingPath"]
                            edge_style = config["edgeStyleMissingPath"]
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
    
    #base label definition
    def _visualize_base_label(self, orientatedBase: str):
        base_label = "<" + "<font point-size=\"6\">"
        if orientatedBase[1]=="forward":
            base_label += orientatedBase[0]+"*"
        elif orientatedBase[1]=="backward":
            base_label += "*"+General.reverse_complement(orientatedBase[0])
        else:
            base_label += "---orientation unknown---"
        base_label +="<br/><font color=\"grey\">"
        base_label +=str(self._orientatedBases[orientatedBase]._number)+"x"
        base_label +="</font></font>" + ">"
        return base_label
    
    #edge label definition
    def _visualize_edge_label(self, info: dict):
        edge_label = "<"
        edge_label += "<font color=\"grey\">"+str(info["distance"])+"</font><br/>"
        edge_label += "<font point-size=\"8\">"+str(info["number"])+"x</font>" 
        edge_label += ">"
        return edge_label
    
    def _visualize_dead_end_label(self, distance: int):
        edge_label = "<"
        edge_label += "<font color=\"grey\">"+str(distance)+"</font>" 
        edge_label += ">"
        return edge_label
    
    #node label definition
    def _visualize_node_label(self, orientatedBase: str, orientatedCkmer: str):
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
        node_label = "<" + "<font color=\"blue\">"+node_letter+"</font>" 
        node_label +="<font point-size=\"8\">"
        node_label +=str(self._orientatedCkmers[orientatedCkmer]._number)+"x</font>"
        node_label += ">"
        return node_label
    
    def unset_start(self, orientatedCkmerKey: str):
        if orientatedCkmerKey in self._start:
            self._start.remove(orientatedCkmerKey)
            #note: don't undo preStart (?)
                
    def set_start(self, orientatedCkmerKey: str):
        if not orientatedCkmerKey in self._start:
            self._start.add(orientatedCkmerKey)
            orientatedCkmer = self._orientatedCkmers[orientatedCkmerKey]
            for incomingOrientatedCkmerKey in orientatedCkmer._incoming.keys():
                self._orientatedCkmers[incomingOrientatedCkmerKey].set_preStart(1)
            
    def unset_end(self, orientatedCkmerKey: str):
        if orientatedCkmerKey in self._end:
            self._end.remove(orientatedCkmerKey)
            #note: don't undo postEnd (?)
            
    def set_end(self, orientatedCkmerKey: str):
        if not orientatedCkmerKey in self._end:
            self._end.add(orientatedCkmerKey)
            orientatedCkmer = self._orientatedCkmers[orientatedCkmerKey]
            for outgoingOrientatedCkmerKey in orientatedCkmer._outgoing.keys():
                self._orientatedCkmers[outgoingOrientatedCkmerKey].set_postEnd(1)
    
    def getCandidates(self):
        candidates = []
        for orientatedCkmer in self._orientatedCkmers.keys():
            if self._orientatedCkmers[orientatedCkmer]._candidate:
                candidates.append(orientatedCkmer)
        return set(candidates)
    
    def getConnectedCandidates(self):
        d = self.getShortestDistances()
        connectedSets = []
        candidates = self.getCandidates()
        #compute connected candidates
        for k in candidates:    
            c = set([m for m in candidates if m==k or not d[k][m]==0])
            newConnectedSets = [c]
            for connectedSet in connectedSets:
                if len(c.intersection(connectedSet))>0:
                    newConnectedSets[0].update(connectedSet)
                else:
                    newConnectedSets.append(connectedSet)
            connectedSets=newConnectedSets
        #construct result
        result = []
        for connectedSet in connectedSets:
            entry = {"connected": connectedSet, "start": set(), "end": set()}
            connectedList = list(connectedSet)
            for k in connectedList:
                if max(d[k][connectedList].values)==0:
                    entry["start"].add(k)
                if min(d[k][connectedList].values)==0:
                    entry["end"].add(k)
            result.append(entry)
        return result
        
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
                if self._distances["directDistance"].values[i,j]>0:
                    nkg.addEdge(i,j,self._distances["directDistance"].values[i,j])
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
            self._incomingDeadEnd : tuple = None
            self._outgoingDeadEnd : tuple = None
            self._orientatedBases : dict = {}
            self._order : int = None
            self._preStart : int = None
            self._postEnd : int = None
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
        
        def set_incoming_dead_end(self, orientatedCkmerKey, distance: int, path: str):
            assert len(self._incoming)==0
            self._incomingDeadEnd = (orientatedCkmerKey, distance, path)
            #consistency check
            if not path==None:
                assert len(orientatedCkmerKey[0])==self._graph._k
                assert len(path)==self._graph._k+distance
                if self._orientation=="forward":
                    assert path[-self._graph._k:]==self._ckmer
                else:
                    assert path[-self._graph._k:]==General.reverse_complement(self._ckmer)
                if orientatedCkmerKey[1]=="forward":
                    assert path[0:self._graph._k]==orientatedCkmerKey[0]
                else:
                    assert path[0:self._graph._k]==General.reverse_complement(orientatedCkmerKey[0])
        
        def set_outgoing_dead_end(self, orientatedCkmerKey, distance: int, path: str):
            assert len(self._outgoing)==0
            self._outgoingDeadEnd = (orientatedCkmerKey, distance, path)
            #consistency check
            if not path==None:
                assert len(orientatedCkmerKey[0])==self._graph._k
                assert len(path)==self._graph._k+distance
                if self._orientation=="forward":
                    assert path[0:self._graph._k]==self._ckmer
                else:
                    assert path[0:self._graph._k]==General.reverse_complement(self._ckmer)
                if orientatedCkmerKey[1]=="forward":
                    assert path[-self._graph._k:]==orientatedCkmerKey[0]
                else:
                    assert path[-self._graph._k:]==General.reverse_complement(orientatedCkmerKey[0])
        
        def set_incoming(self, orientatedCkmerKey, distance: int, number: int, problem: bool = None, path: str = None):
            assert orientatedCkmerKey in self._graph._orientatedCkmers.keys()
            assert distance>0
            assert number>=0
            assert path == None or len(path)==distance+self._graph._k
            assert self._incomingDeadEnd == None
            if distance<=self._graph._k:
                #check shared
                if distance<self._graph._k:
                    if self._orientation=="forward":
                        shared = self._ckmer[0:self._graph._k-distance]
                    else:
                        shared = General.reverse_complement(self._ckmer)[0:self._graph._k-distance]
                    if orientatedCkmerKey[1]=="forward":
                        assert shared==orientatedCkmerKey[0][-len(shared):]
                    else:
                        assert shared==General.reverse_complement(orientatedCkmerKey[0])[-len(shared):]
                #compute path and compare with provided
                if orientatedCkmerKey[1]=="forward":
                    computedPath = orientatedCkmerKey[0]
                else:
                    computedPath = General.reverse_complement(orientatedCkmerKey[0])
                if self._orientation=="forward":
                    computedPath = computedPath + self._ckmer[-distance:]
                else:
                    computedPath = computedPath + General.reverse_complement(self._ckmer)[-distance:]
                if not path == None:
                    assert path == computedPath
                else:
                    path = computedPath
            #set or update connection
            if orientatedCkmerKey in self._outgoing.keys():
                assert self._incoming[orientatedCkmerKey]["distance"] == distance
                assert self._incoming[orientatedCkmerKey]["number"] == number
            else:    
                self._incoming[orientatedCkmerKey] = {
                    "distance": distance, "number": number, 
                    "problem": problem, "path": None
                }
            if not problem==None:
                self._incoming[orientatedCkmerKey]["problem"] = problem
            if not path==None:
                self._incoming[orientatedCkmerKey]["path"] = path    
            #update start and end based on connections
            if orientatedCkmerKey in self._graph._end:
                self.set_postEnd(1)
            elif not self._graph._orientatedCkmers[orientatedCkmerKey]._postEnd==None:
                self.set_postEnd(self._graph._orientatedCkmers[orientatedCkmerKey]._postEnd+1)
            if self._key in self._graph._start:
                self._graph._orientatedCkmers[orientatedCkmerKey].set_preStart(1)
            elif not self._preStart == None:
                self._graph._orientatedCkmers[orientatedCkmerKey].set_preStart(self._preStart+1)
            
        def set_outgoing(self, orientatedCkmerKey, distance: int, number: int, problem: bool = None, path: str = None):
            assert orientatedCkmerKey in self._graph._orientatedCkmers.keys()
            assert distance>0
            assert number>=0
            assert path == None or len(path)==distance+self._graph._k
            assert self._outgoingDeadEnd == None
            if distance<=self._graph._k:
                #check shared
                if distance<self._graph._k:
                    if self._orientation=="forward":
                        shared = self._ckmer[-(self._graph._k-distance):]
                    else:
                        shared = General.reverse_complement(self._ckmer)[-(self._graph._k-distance):]
                    if orientatedCkmerKey[1]=="forward":
                        assert shared==orientatedCkmerKey[0][0:len(shared)]
                    else:
                        assert shared==General.reverse_complement(orientatedCkmerKey[0])[0:len(shared)]
                #compute path and compare with provided
                if self._orientation=="forward":
                    computedPath = self._ckmer
                else:
                    computedPath = General.reverse_complement(self._ckmer)
                if orientatedCkmerKey[1]=="forward":
                    computedPath = computedPath + orientatedCkmerKey[0][-distance:]
                else:
                    computedPath = computedPath + General.reverse_complement(orientatedCkmerKey[0])[-distance:]
                if not path == None:
                    assert path == computedPath
                else:
                    path = computedPath
            if orientatedCkmerKey in self._outgoing.keys():
                assert self._outgoing[orientatedCkmerKey]["distance"] == distance
                assert self._outgoing[orientatedCkmerKey]["number"] == number
            else:    
                self._outgoing[orientatedCkmerKey] = {
                    "distance": distance, "number": number, 
                    "problem": problem, "path": None
                }
            if not problem==None:
                self._outgoing[orientatedCkmerKey]["problem"] = problem
            if not path==None:
                self._outgoing[orientatedCkmerKey]["path"] = path    
            #update start and end based on connections
            if orientatedCkmerKey in self._graph._start:
                self.set_preStart(1)
            elif not self._graph._orientatedCkmers[orientatedCkmerKey]._preStart==None:
                self.set_preStart(self._graph._orientatedCkmers[orientatedCkmerKey]._preStart+1)
            if self._key in self._graph._end:
                self._graph._orientatedCkmers[orientatedCkmerKey].set_postEnd(1)
            elif not self._postEnd == None:
                self._graph._orientatedCkmers[orientatedCkmerKey].set_postEnd(self._postEnd+1)
            
        def unset_candidate(self):
            self._candidate = False
            #note: don't undo base setting(?)
                
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
                
        def set_preStart(self, steps: int):
            if (self._preStart == None) or (self._preStart>steps):
                self._preStart = steps
                for direction in self._orientatedBases.keys():
                    orientatedBaseKey = self._orientatedBases[direction]
                    self._graph._orientatedBases[orientatedBaseKey].set_preStart(steps)
                for incomingOrientatedCkmerKey in self._incoming.keys():
                    self._graph._orientatedCkmers[incomingOrientatedCkmerKey].set_preStart(steps+1)
                    
        def set_postEnd(self, steps: int):
            if (self._postEnd == None) or (self._postEnd>steps):
                self._postEnd = steps
                for direction in self._orientatedBases.keys():
                    orientatedBaseKey = self._orientatedBases[direction]
                    self._graph._orientatedBases[orientatedBaseKey].set_postEnd(steps)
                for outgoingOrientatedCkmerKey in self._outgoing.keys():
                    self._graph._orientatedCkmers[outgoingOrientatedCkmerKey].set_postEnd(steps+1)
             
            
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
            self._number : int = 0
            self._order : int = None
            self._orientatedCkmers : set = set()
            self._linked : set = set()
            self._preStart : int = None
            self._postEnd : int = None
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
            
        def set_preStart(self,steps):
            if (self._preStart == None) or (self._preStart>steps):
                self._preStart = steps
                for orientatedCkmerKey in self._orientatedCkmers:
                    self._graph._orientatedCkmers[orientatedCkmerKey].set_preStart(steps)

        def set_postEnd(self,steps):
            if (self._postEnd == None) or (self._postEnd>steps):
                self._postEnd = steps   
                for orientatedCkmerKey in self._orientatedCkmers:
                    self._graph._orientatedCkmers[orientatedCkmerKey].set_postEnd(steps)
            
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
                            
