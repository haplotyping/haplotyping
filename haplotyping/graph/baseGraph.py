from haplotyping.general import General
from graphviz import Digraph
import html
import pandas as pd
import networkit as nk
import numpy as np
from typing import Literal

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
        self._arms = []
            
    def __repr__(self):
        text = "Graph object"
        if self._name:
            text = "{} {}".format(text,self._name)
        text = "{} with {} candidate k-mers".format(text,len(self.getCandidates()))
        return text
        
    def _setName(self, name: str):
        self._name = name
    
    def name(self):
        return self._name
        
    def visualize(self, *args, **kwargs):
        
        #initialise configuration
        config = {
            "showAllBases": False,
            "showPrePostBases": True,
            "showPrePostBasesMaxSteps": 2,
            "showDeadEnds": True,
            "hideDeadEndBefore": [],
            "hideDeadEndAfter": [],
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
            "nodeFillColorIncomingArm": "lightblue",
            "nodeColorIncomingArm": "black",
            "nodeFillColorOutgoingArm": "lightcoral",
            "nodeColorOutgoingArm": "black",
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
            "restrictedListOfOrientatedCkmers": False,
            "containerGraph": False,
            "prefix": "graph"
        }
        for key, value in kwargs.items():
            if key in config.keys():
                config[key] = value
        
        #create graph
        if config["containerGraph"]:
            g = config["containerGraph"]
        else:
            g = Digraph()
            graph_label = "Graph"
            if self._name:
                graph_label = "{} {}".format(graph_label,html.escape(self._name))
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
        maxOrientatedBaseOrder = (0 if len(self._orientatedBases)==0 
                                  else max([self._orientatedBases[orientatedBase]._order or 0 
                                      for orientatedBase in self._orientatedBases.keys()]))
        sortedOrientatedBaseList = sorted(self._orientatedBases.keys(), 
                                key=lambda orientatedBase: maxOrientatedBaseOrder+1 
                                          if self._orientatedBases[orientatedBase]._order==None 
                                          else self._orientatedBases[orientatedBase]._order)
        
        #restrict bases
        if config["restrictedListOfOrientatedCkmers"]:
            restrictedListOfBases = set()
            for orientatedCkmer in config["restrictedListOfOrientatedCkmers"]:
                restrictedListOfBases.update(self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
            sortedOrientatedBaseList = [base for base in sortedOrientatedBaseList if base in restrictedListOfBases]
            
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
                elif self._orientatedBases[orientatedBase].candidate():
                    pass
                else:
                    continue
            
            #create container for orientatedBaseContainer (multiple bases can be connected by shared k-mers)
            #and make container visible if it contains multiple bases
            if not orientatedBase in orientatedBaseContainers.keys():
                orientatedBaseContainerName = "cluster_container_{}_{}_{}".format(config["prefix"],
                                                          orientatedBase[0],orientatedBase[1])
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
                with obc.subgraph(name="cluster_{}_{}_{}".format(config["prefix"],
                                                         orientatedBase[0],orientatedBase[1])) as c: 
                    #define layout
                    base_style=config["baseStyleDefault"]
                    base_fillcolor=config["baseFillColorDefault"]
                    base_color=config["baseColorDefault"]
                    base_penwidth=config["basePenWidthDefault"]
                    if self._orientatedBases[orientatedBase].candidate():
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
                            if self._orientatedCkmers[orientatedCkmer].candidate():
                                node_fillcolor=config["nodeFillColorCandidate"]
                                node_color=config["nodeColorCandidate"]
                            elif self._orientatedCkmers[orientatedCkmer].incomingArmType():
                                node_fillcolor=config["nodeFillColorIncomingArm"]
                                node_color=config["nodeColorIncomingArm"]
                            elif self._orientatedCkmers[orientatedCkmer].outgoingArmType():
                                node_fillcolor=config["nodeFillColorOutgoingArm"]
                                node_color=config["nodeColorOutgoingArm"]
                        if orientatedCkmer in self._start or orientatedCkmer in self._end:
                            node_fillcolor=config["nodeFillColorStartEnd"]
                            node_color=config["nodeColorStartEnd"]
                        if orientatedCkmer in self._selected:
                            node_color=config["nodeColorSelected"]
                            node_penwidth=config["nodePenWidthSelected"]
                        node_label = self._visualize_node_label(orientatedBase, orientatedCkmer)
                        node_key = "{}_{}_{}_{}_{}".format(config["prefix"],
                                orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                        c.node(node_key, label=node_label, color=node_color, style=node_style, 
                               fillcolor=node_fillcolor, penwidth=str(node_penwidth))
                        if ((config["showAllBases"] or config["showDeadEnds"]) 
                            and not orientatedCkmer in config["hideDeadEndBefore"]
                            and not self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd==None):
                            dead_key = "{}_{}_{}_{}".format(config["prefix"],
                                orientatedBase[0],orientatedBase[1],
                                self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[0])
                            g.node(dead_key, shape="point")
                            g.edge(dead_key, node_key, constraint="true",
                              label=self._visualize_dead_end_label(
                                  self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[1]),
                                   dir="forward", color=config["edgeColorAuto"], weight="1",
                                   style=config["edgeStyleAuto"], penwidth=str(config["edgePenWidthAuto"]))
                        if ((config["showAllBases"] or config["showDeadEnds"])
                            and not orientatedCkmer in config["hideDeadEndAfter"]
                            and not self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd==None):
                            dead_key = "{}_{}_{}_{}".format(config["prefix"],
                                orientatedBase[0],orientatedBase[1],
                                self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[0])
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
        maxOrientatedCkmerOrder = (0 if len(self._orientatedBases)==0 
                                  else max([self._orientatedBases[orientatedBase]._order or 0 
                                      for orientatedBase in self._orientatedBases.keys()]))
        sortedOrientatedCkmerList = sorted(self._orientatedCkmers.keys(), 
                                key=lambda orientatedCkmer: maxOrientatedCkmerOrder+1 
                                           if self._orientatedCkmers[orientatedCkmer]._order==None 
                                           else self._orientatedCkmers[orientatedCkmer]._order)

        
        #restrict k-mers
        if config["restrictedListOfOrientatedCkmers"]:
            sortedOrientatedCkmerList = [ckmer for ckmer in sortedOrientatedCkmerList 
                                         if ckmer in config["restrictedListOfOrientatedCkmers"]]
        
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
    
    def _unsetStart(self, orientatedCkmer: str):
        if orientatedCkmer in self._start:
            self._start.remove(orientatedCkmer)
            #note: don't undo preStart (?)
                
    def _setStart(self, orientatedCkmer: str):
        self._orientatedCkmers[orientatedCkmer]._setCandidate()
        if not orientatedCkmer in self._start:
            self._start.add(orientatedCkmer)
            for incomingOrientatedCkmer in self._orientatedCkmers[orientatedCkmer]._incoming.keys():
                self._orientatedCkmers[incomingOrientatedCkmer]._setPreStart(1)
            
    def _unsetEnd(self, orientatedCkmer: str):
        if orientatedCkmer in self._end:
            self._end.remove(orientatedCkmer)
            #note: don't undo postEnd (?)
            
    def _setEnd(self, orientatedCkmer: str):
        self._orientatedCkmers[orientatedCkmer]._setCandidate()
        if not orientatedCkmer in self._end:
            self._end.add(orientatedCkmer)
            for outgoingOrientatedCkmer in self._orientatedCkmers[orientatedCkmer]._outgoing.keys():
                self._orientatedCkmers[outgoingOrientatedCkmer]._setPostEnd(1)
    
    def getOrientatedCkmers(self):
        return set(self._orientatedCkmers.keys())
    
    def getStartCandidates(self):
        candidates = []
        for orientatedCkmer in self._orientatedCkmers.keys():
            if orientatedCkmer in self._start and self._orientatedCkmers[orientatedCkmer].candidate():
                candidates.append(orientatedCkmer)
        return set(candidates)
    
    def getEndCandidates(self):
        candidates = []
        for orientatedCkmer in self._orientatedCkmers.keys():
            if orientatedCkmer in self._end and self._orientatedCkmers[orientatedCkmer].candidate():
                candidates.append(orientatedCkmer)
        return set(candidates)
    
    def getCandidates(self):
        candidates = []
        for orientatedCkmer in self._orientatedCkmers.keys():
            if self._orientatedCkmers[orientatedCkmer].candidate():
                candidates.append(orientatedCkmer)
        return set(candidates)
    
    """get orientated candidate bases"""
    def getCandidateBases(self):
        candidateBases = set()
        for orientatedCkmer in self._orientatedCkmers.keys():
            if self._orientatedCkmers[orientatedCkmer].candidate():
                candidateBases.update(self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
        return candidateBases
    
    """get connected sets of orientated candidate k-mers"""
    def getConnectedCandidates(self, baseConnected=False):
        d = self.getShortestDistances()
        connectedSets = []
        candidates = self.getCandidates()
        #compute connected candidates
        for k in candidates:    
            if baseConnected:
                kmerSet = set([k])
                for orientatedBase in self._orientatedCkmers[k]._orientatedBases.values():
                    kmerSet.update(self._orientatedBases[orientatedBase]._orientatedCkmers)
                c = set()
                for ks in kmerSet:
                    c.update([m for m in candidates if m==ks or not d[ks][m]==0])
            else:
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
                if k in self._start:
                    entry["start"].add(k)
                elif max(d[k][connectedList].values)==0:
                    self._logger.debug("unexpected start entry {}".format(k))
                    entry["start"].add(k)
                if k in self._end:
                    entry["end"].add(k)
                elif min(d[k][connectedList].values)==0:
                    self._logger.debug("unexpected end entry {}".format(k))
                    entry["end"].add(k)
            result.append(entry)
        return result
        
    """get direct distances between orientated k-mers"""
    def getDirectDistances(self):
        if self._distances == None:
            self._computeDistances()
        return self._distances["directDistance"]
        
    """get frequencies read evidence of direct connections between orientated k-mers"""
    def getDirectNumbers(self):
        if self._distances == None:
            self._computeDistances()
        return self._distances["directNumber"]    
    
    """get shortest distances between orientated k-mers"""
    def getShortestDistances(self):
        if self._distances == None:
            self._computeDistances()
        return self._distances["shortest"]
        
    def getShortestDistance(self, fromOrientatedKmers: set, toOrientatedKmers: set):
        """
        get the shortest distance between two sets of orientated k-mers
        """
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
        
    def getArms(self):
        return self._arms
    
    def detectArms(self):
        """
        detect arms
        """
        if len(self._arms)==0:
            candidates = list(self.getCandidates())
            candidateBases = list(self.getCandidateBases())
            dist = self.getShortestDistances()
            other = set(self._orientatedCkmers.keys()).difference(candidates)
            incomingArms = {}
            outgoingArms = {}
            for orientatedCkmer in other:
                entryDistances = dist[orientatedCkmer][candidates]
                incomingDistances = [d for d in entryDistances if d<0]
                outgoingDistances = [d for d in entryDistances if d>0]
                if len(incomingDistances)>0:
                    d = max(incomingDistances)
                    for candidateCkmer,candidateDistance in entryDistances.items():
                        if candidateDistance==d:
                            if candidateCkmer in incomingArms:
                                incomingArms[candidateCkmer]["size"] = max(incomingArms[candidateCkmer]["size"],d)
                                incomingArms[candidateCkmer]["number"] = max(incomingArms[candidateCkmer]["number"],
                                                                              self._orientatedCkmers[orientatedCkmer]._number)
                                incomingArms[candidateCkmer]["n"]+=1
                                incomingArms[candidateCkmer]["orientatedBases"].update(
                                    self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
                                incomingArms[candidateCkmer]["orientatedCkmers"].add(orientatedCkmer)
                            else:
                                incomingArms[candidateCkmer] = {"size": d, "n": 1, 
                                        "number": self._orientatedCkmers[orientatedCkmer]._number,
                                        "orientatedBases": set(
                                            self._orientatedCkmers[orientatedCkmer]._orientatedBases.values()),
                                        "orientatedCkmers": set([orientatedCkmer])}
                if len(outgoingDistances)>0:
                    d = min(outgoingDistances)
                    for candidateCkmer,candidateDistance in entryDistances.items():
                        if candidateDistance==d:
                            if candidateCkmer in outgoingArms:
                                outgoingArms[candidateCkmer]["size"] = max(outgoingArms[candidateCkmer]["size"],d)
                                outgoingArms[candidateCkmer]["number"] = max(outgoingArms[candidateCkmer]["number"],
                                                                              self._orientatedCkmers[orientatedCkmer]._number)
                                outgoingArms[candidateCkmer]["n"]+=1
                                outgoingArms[candidateCkmer]["orientatedBases"].update(
                                    self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
                                outgoingArms[candidateCkmer]["orientatedCkmers"].add(orientatedCkmer)
                            else:
                                outgoingArms[candidateCkmer] = {"size": d, "n": 1, 
                                        "number": self._orientatedCkmers[orientatedCkmer]._number,
                                        "orientatedBases": set(
                                            self._orientatedCkmers[orientatedCkmer]._orientatedBases.values()),
                                        "orientatedCkmers": set([orientatedCkmer])}
            for candidateCkmer in incomingArms:
                if candidateCkmer in self._start:
                    continue
                elif len(incomingArms[candidateCkmer]["orientatedBases"].difference(candidateBases))>0:
                    arm = self.Arm(self, candidateCkmer, "incoming")
                    for orientatedCkmer in incomingArms[candidateCkmer]["orientatedCkmers"]:
                        arm._add(orientatedCkmer)
            for candidateCkmer in outgoingArms:
                if candidateCkmer in self._end:
                    continue
                elif len(outgoingArms[candidateCkmer]["orientatedBases"].difference(candidateBases))>0:
                    arm = self.Arm(self, candidateCkmer, "outgoing")
                    for orientatedCkmer in outgoingArms[candidateCkmer]["orientatedCkmers"]:
                        arm._add(orientatedCkmer)

    
    def resetArms(self):
        for orientatedCkmer in self._orientatedCkmers:
            self._orientatedCkmers[orientatedCkmer]._incomingArm = None
            self._orientatedCkmers[orientatedCkmer]._outgoingArm = None
        for orientatedCkmer in self._orientatedCkmers:
            self._orientatedCkmers[orientatedCkmer]._unsetArm()          
        self._arms = []
                
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
        
        TYPE = Literal["candidate","incomingArm","outgoingArm"]

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
            self._type : TYPE = None
            self._arm = []
            self._incoming : dict = {}
            self._outgoing : dict = {}
            self._incomingDeadEnd : tuple = None
            self._outgoingDeadEnd : tuple = None
            self._incomingArm = None
            self._outgoingArm = None
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
                self._graph._orientatedBases[orientatedBaseKey]._addOrientatedCkmer(self._ckmer, self._orientation)      
            if (self._split == "both") or (self._split == "left"):
                base = General.reverse_complement(self._ckmer[1:])
                baseOrientation = "forward" if self._orientation=="backward" else "backward"
                orientatedBaseKey = (base,baseOrientation)
                self._orientatedBases["left"] = orientatedBaseKey
                if not orientatedBaseKey in self._graph._orientatedBases.keys():
                    orientatedBaseEntry = self._graph.OrientatedBase(self._graph, base, baseOrientation)
                self._graph._orientatedBases[orientatedBaseKey]._addOrientatedCkmer(self._ckmer, self._orientation)     
            
        def __repr__(self):
            info = []
            info.append("{}x".format(self._number))
            info.append("split {}".format(self._split))
            if self._type:
                info.append("type: {}".format(self._type))
            if self._order:
                info.append("order: {}".format(self._order))
            return "OrientatedCkmer({} {}[{}])".format(self._ckmer,self._orientation,", ".join(info))
        
        def _setIncomingArm(self, arm):
            assert self._incomingArm == None
            self._incomingArm = arm
        
        def _setOutgoingArm(self, arm):
            assert self._outgoingArm == None
            self._outgoingArm = arm
            
        def _setIncomingDeadEnd(self, orientatedCkmer, distance: int, path: str):
            assert len(self._incoming)==0
            self._incomingDeadEnd = (orientatedCkmer, distance, path)
            #consistency check
            if not path==None:
                assert len(orientatedCkmer[0])==self._graph._k
                assert len(path)==self._graph._k+distance
                if self._orientation=="forward":
                    assert path[-self._graph._k:]==self._ckmer
                else:
                    assert path[-self._graph._k:]==General.reverse_complement(self._ckmer)
                if orientatedCkmer[1]=="forward":
                    assert path[0:self._graph._k]==orientatedCkmer[0]
                else:
                    assert path[0:self._graph._k]==General.reverse_complement(orientatedCkmer[0])
        
        def _setOutgoingDeadEnd(self, orientatedCkmer, distance: int, path: str):
            assert len(self._outgoing)==0
            self._outgoingDeadEnd = (orientatedCkmer, distance, path)
            #consistency check
            if not path==None:
                assert len(orientatedCkmer[0])==self._graph._k
                assert len(path)==self._graph._k+distance
                if self._orientation=="forward":
                    assert path[0:self._graph._k]==self._ckmer
                else:
                    assert path[0:self._graph._k]==General.reverse_complement(self._ckmer)
                if orientatedCkmer[1]=="forward":
                    assert path[-self._graph._k:]==orientatedCkmer[0]
                else:
                    assert path[-self._graph._k:]==General.reverse_complement(orientatedCkmer[0])
        
        def _setIncoming(self, orientatedCkmer, distance: int, number: int, problem: bool = None, path: str = None):
            assert orientatedCkmer in self._graph._orientatedCkmers
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
                    if orientatedCkmer[1]=="forward":
                        assert shared==orientatedCkmer[0][-len(shared):]
                    else:
                        assert shared==General.reverse_complement(orientatedCkmer[0])[-len(shared):]
                #compute path and compare with provided
                if orientatedCkmer[1]=="forward":
                    computedPath = orientatedCkmer[0]
                else:
                    computedPath = General.reverse_complement(orientatedCkmer[0])
                if self._orientation=="forward":
                    computedPath = computedPath + self._ckmer[-distance:]
                else:
                    computedPath = computedPath + General.reverse_complement(self._ckmer)[-distance:]
                if not path == None:
                    assert path == computedPath
                else:
                    path = computedPath
            #set or update connection
            if orientatedCkmer in self._outgoing.keys():
                if not self._incoming[orientatedCkmer]["problem"]:
                    assert self._incoming[orientatedCkmer]["distance"] == distance
                    assert self._incoming[orientatedCkmer]["number"] == number
            else:    
                self._incoming[orientatedCkmer] = {
                    "distance": distance, "number": number, 
                    "problem": problem, "path": None
                }
            if not problem==None:
                self._incoming[orientatedCkmer]["problem"] = problem
            if not path==None:
                self._incoming[orientatedCkmer]["path"] = path    
            #update start and end based on connections
            if orientatedCkmer in self._graph._end:
                self._setPostEnd(1)
            elif not self._graph._orientatedCkmers[orientatedCkmer]._postEnd==None:
                self._setPostEnd(self._graph._orientatedCkmers[orientatedCkmer]._postEnd+1)
            if self._key in self._graph._start:
                self._graph._orientatedCkmers[orientatedCkmer]._setPreStart(1)
            elif not self._preStart == None:
                self._graph._orientatedCkmers[orientatedCkmer]._setPreStart(self._preStart+1)
            
        def _setOutgoing(self, orientatedCkmer, distance: int, number: int, problem: bool = None, path: str = None):
            assert orientatedCkmer in self._graph._orientatedCkmers.keys()
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
                    if orientatedCkmer[1]=="forward":
                        assert shared==orientatedCkmer[0][0:len(shared)]
                    else:
                        assert shared==General.reverse_complement(orientatedCkmer[0])[0:len(shared)]
                #compute path and compare with provided
                if self._orientation=="forward":
                    computedPath = self._ckmer
                else:
                    computedPath = General.reverse_complement(self._ckmer)
                if orientatedCkmer[1]=="forward":
                    computedPath = computedPath + orientatedCkmer[0][-distance:]
                else:
                    computedPath = computedPath + General.reverse_complement(orientatedCkmer[0])[-distance:]
                if not path == None:
                    assert path == computedPath
                else:
                    path = computedPath
            if orientatedCkmer in self._outgoing.keys():
                if not self._outgoing[orientatedCkmer]["problem"]:
                    assert self._outgoing[orientatedCkmer]["distance"] == distance
                    assert self._outgoing[orientatedCkmer]["number"] == number
            else:    
                self._outgoing[orientatedCkmer] = {
                    "distance": distance, "number": number, 
                    "problem": problem, "path": None
                }
            if not problem==None:
                self._outgoing[orientatedCkmer]["problem"] = problem
            if not path==None:
                self._outgoing[orientatedCkmer]["path"] = path    
            #update start and end based on connections
            if orientatedCkmer in self._graph._start:
                self._setPreStart(1)
            elif not self._graph._orientatedCkmers[orientatedCkmer]._preStart==None:
                self._setPreStart(self._graph._orientatedCkmers[orientatedCkmer]._preStart+1)
            if self._key in self._graph._end:
                self._graph._orientatedCkmers[orientatedCkmer]._setPostEnd(1)
            elif not self._postEnd == None:
                self._graph._orientatedCkmers[orientatedCkmer]._setPostEnd(self._postEnd+1)
            
        def _unsetCandidate(self):
            if self.candidate():
                self._type = None
                for orientatedBase in self._orientatedBases.values():
                    self._graph._orientatedBases[orientatedBase]._setType()                                
                
        def _setCandidate(self):
            self._type = "candidate"
            for orientatedBase in self._orientatedBases.values():
                self._graph._orientatedBases[orientatedBase]._setType()                                                
                
        def candidate(self):
            return self._type=="candidate"
                
        def _unsetArm(self):
            if self.arm():
                self._type = None
                for orientatedBase in self._orientatedBases.values():
                    self._graph._orientatedBases[orientatedBase]._setType() 
                self._arm = []
                
        def _setArm(self, arm):
            if arm.armType()=="incoming":
                self._type = "incomingArm"
                if len(self._arm)>0:
                    for otherArm in self._arm:
                        assert otherArm.armType()=="incoming"
                self._arm.append(arm)
            elif arm.armType()=="outgoing":
                self._type = "outgoingArm"
                if len(self._arm)>0:
                    for otherArm in self._arm:
                        assert otherArm.armType()=="outgoing"
                self._arm.append(arm)
            for orientatedBase in self._orientatedBases.values():
                self._graph._orientatedBases[orientatedBase]._setType()
                
        def arm(self):
            return self._type=="incomingArm" or self._type=="outgoingArm"
                
        def armKey(self):
            if self.arm():
                return [arm.key() for arm in self._arm]
            else:
                return None
                
        def incomingArmType(self):
            return self._type=="incomingArm"
                
        def outgoingArmType(self):
            return self._type=="outgoingArm"
                
        def incomingArm(self):
            return self._incomingArm
                
        def outgoingArm(self):
            return self._outgoingArm
                
        def _setOrder(self, order: int):
            self._order = order
            for direction in self._orientatedBases.keys():
                orientatedBaseKey = self._orientatedBases[direction]
                self._graph._orientatedBases[orientatedBaseKey]._setOrder()
                
        def order(self):
            return self._order
                
        def _setPreStart(self, steps: int):
            if (self._preStart == None) or (self._preStart>steps):
                self._preStart = steps
                for direction in self._orientatedBases.keys():
                    orientatedBaseKey = self._orientatedBases[direction]
                    self._graph._orientatedBases[orientatedBaseKey]._setPreStart(steps)
                for incomingOrientatedCkmer in self._incoming.keys():
                    self._graph._orientatedCkmers[incomingOrientatedCkmer]._setPreStart(steps+1)
                    
        def _setPostEnd(self, steps: int):
            if (self._postEnd == None) or (self._postEnd>steps):
                self._postEnd = steps
                for direction in self._orientatedBases.keys():
                    orientatedBaseKey = self._orientatedBases[direction]
                    self._graph._orientatedBases[orientatedBaseKey]._setPostEnd(steps)
                for outgoingOrientatedCkmer in self._outgoing.keys():
                    self._graph._orientatedCkmers[outgoingOrientatedCkmer]._setPostEnd(steps+1)
             
            
    class OrientatedBase:
        """
        internal object representing orientated base in the De Bruijn graph
        """
        
        TYPE = Literal["candidate","arm"]
        
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
            self._type : TYPE = None
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
            info.append("{}x".format(self._number))
            info.append("{} k-mers".format(len(self._orientatedCkmers)))
            if self._type:
                info.append("type: {}".format(self._type))
            if self._order:
                info.append("order: {}".format(self._order))
            return "OrientatedBase({} {} [{}])".format(self._base,self._orientation,", ".join(info))
        
        def _setType(self):
            self._type = None
            for orientatedCkmer in self._orientatedCkmers:
                if self._graph._orientatedCkmers[orientatedCkmer].candidate():
                    self._type="candidate"
                    self._preStart = None
                    break
                elif self._graph._orientatedCkmers[orientatedCkmer].arm():
                    self._type="arm"
            
        def candidate(self):
            return self._type == "candidate"
            
        def arm(self):
            return self._type == "arm"
            
        def _setPreStart(self,steps):
            if self.candidate():
                pass
            elif (self._preStart == None) or (self._preStart>steps):
                self._preStart = steps
                for orientatedCkmer in self._orientatedCkmers:
                    self._graph._orientatedCkmers[orientatedCkmer]._setPreStart(steps)

        def _setPostEnd(self,steps):
            if self.candidate():
                pass
            elif (self._postEnd == None) or (self._postEnd>steps):
                self._postEnd = steps   
                for orientatedCkmer in self._orientatedCkmers:
                    self._graph._orientatedCkmers[orientatedCkmer]._setPostEnd(steps)
            
        def _setOrder(self):
            self._order = None
            for orientatedCkmer in self._orientatedCkmers:
                if not self._graph._orientatedCkmers[orientatedCkmer]._order == None:
                    if self._order == None:
                        self._order = self._graph._orientatedCkmers[orientatedCkmer]._order
                    else:
                        self._order = min(self._order, self._graph._orientatedCkmers[orientatedCkmer]._order)
                                           
        
        def _addOrientatedCkmer(self, ckmer, orientation):
            orientatedCkmer = (ckmer, orientation)            
            if orientatedCkmer in self._graph._orientatedCkmers.keys():
                assert orientatedCkmer not in self._orientatedCkmers
                self._orientatedCkmers.add(orientatedCkmer)
                self._number += self._graph._orientatedCkmers[orientatedCkmer]._number
                #check type
                self._setType()
                #check order
                if not self._graph._orientatedCkmers[orientatedCkmer]._order == None:
                    if self._order == None:
                        self._order = self._graph._orientatedCkmers[orientatedCkmer]._order
                    else:
                        self._order = min(self._order, self._graph._orientatedCkmers[orientatedCkmer]._order)
                #check linked
                for baseKey in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                    if not (baseKey == self._key):
                        self._linked.add(baseKey)
                        self._graph._orientatedBases[baseKey]._linked.add(self._key)
            else:
                raise Exception("can't add non-existing orientated canonical k-mer")
                            
    class Arm:
        """
        internal object representing arm
        """
        
        TYPE = Literal["incoming","outward"]
        
        def __init__(self, graph, orientatedCkmer, type: TYPE):
            self._type = type
            self._graph = graph
            self._orientatedCkmers = set()
            self._connection = None
            self._maxFreq = 0
            self._size = 0
            self._key = None
            #mark incoming/outgoing arm
            assert orientatedCkmer in self._graph._orientatedCkmers
            if type=="incoming":
                assert self._graph._orientatedCkmers[orientatedCkmer].candidate()
                self._graph._orientatedCkmers[orientatedCkmer]._setIncomingArm(self)
                self._connection = orientatedCkmer                
            elif type=="outgoing":
                assert self._graph._orientatedCkmers[orientatedCkmer].candidate()
                self._graph._orientatedCkmers[orientatedCkmer]._setOutgoingArm(self)
                self._connection = orientatedCkmer                      
            #register arm
            self._graph._arms.append(self)
            self._key = len(self._graph._arms)
            
        def __repr__(self):
            if self._type=="incoming":
                text = "Incoming Arm"
            elif self._type=="outgoing":
                text = "Outgoing Arm"
            else:
                text = "Arm"
            text = "{}, {} entries, maximum frequency {}".format(text,self._size,self._maxFreq)
            return text
            
        def _add(self, orientatedCkmer):
            assert orientatedCkmer in self._graph._orientatedCkmers
            if self._type == "incoming":
                self._graph._orientatedCkmers[orientatedCkmer]._setArm(self)
            elif self._type == "outgoing":
                self._graph._orientatedCkmers[orientatedCkmer]._setArm(self)
            self._orientatedCkmers.add(orientatedCkmer)
            self._maxFreq = max(self._maxFreq,self._graph._orientatedCkmers[orientatedCkmer]._number)
            distances = self._graph.getShortestDistances()
            self._size = max(self._size, distances[orientatedCkmer][self._connection], 
                             distances[self._connection][orientatedCkmer])
            
        def armType(self):
            return self._type
        
        def maxFreq(self):
            return self._maxFreq
        
        def key(self):
            return self._key
        
        def n(self):
            return len(self._orientatedCkmers)
        
        def size(self):
            return self._size
        
        def connection(self):
            return self._connection
            
            
            