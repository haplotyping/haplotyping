from haplotyping.general import General
from graphviz import Digraph
import html
import pandas as pd
import networkit as nk
import networkx as nx
import numpy as np
from typing import Literal

class Graph():
    """basic or minimal version of the De Bruijn graph"""
    def __init__(self):
        self._name : str = None
        self._k : int = None
        self._start : set = set()
        self._end : set = set()
        self._ckmers : dict = {}
        self._orientatedCkmers : dict = {}
        self._orientatedBases : dict = {}
        self._selected : set = set()
        self._connected = None
        self._directDistances = None
        self._directNumbers = None
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
        """visualize the De Bruijn graph"""
        #initialise configuration
        initConfig = {
            "type": "basic",
            "showDeadEnds": True,
            "showArms": False,
            "showPotentialConnectedArms": False,
            "hideDeadEndBefore": [],
            "hideDeadEndAfter": [],
            "nodeFillColorDefault": "white",
            "nodeFillColorSelection": "orange",
            "nodeColorDefault": "grey",
            "nodePenWidthDefault": 1,
            "nodeColorSelected": "blue",
            "nodePenWidthSelected": 3,
            "nodeFillColorIncomingArm": "lightblue",
            "nodeColorIncomingArm": "black",
            "nodeFillColorOutgoingArm": "lightcoral",
            "nodeColorOutgoingArm": "black",
            "nodeFillColorCandidate": "lightyellow1",
            "nodeFillColorSelectionCandidate": "red",
            "nodeColorCandidate": "black",
            "nodeFillColorStartEnd": "yellow",
            "nodeFillColorSelectionStartEnd": "darkred",
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
            "selection": False,
            "containerGraph": False,
            "prefix": "graph"
        }
        #construct config
        config = kwargs.copy()
        for key, value in initConfig.items():
            if not key in config.keys():
                config[key] = initConfig[key]

        #translate selection to orientated k-mers
        config["selectionOrientatedCkmers"] = False
        selection = set()
        if config["selection"]:
            if isinstance(config["selection"],str):
                for i in range(len(config["selection"])-(self._k-1)):
                    kmer = config["selection"][i:i+self._k]
                    ckmer = General.canonical(kmer)
                    if ckmer in self._ckmers:
                        selection.update(self._ckmers[ckmer]._orientated)
            elif hasattr(config["selection"], "__iter__"):
                for kmer in config["selection"]:
                    if isinstance(kmer,str) and (len(kmer)==self._k):
                        ckmer = General.canonical(kmer)
                        if ckmer in self._ckmers.keys():
                            selection.update(self._ckmers[ckmer]._orientated)
                    elif isinstance(kmer,tuple) and (len(kmer)==2) and (len(kmer[0])==self._k):
                        if kmer in self._orientatedCkmers.keys():
                            selection.add(kmer)
            if len(selection)>0:
                config["selectionOrientatedCkmers"] = list(selection)

        #create graph
        if config["containerGraph"]:
            g = config["containerGraph"]
        else:
            g = Digraph()
            graph_label = "Graph"
            if self._name:
                graph_label = "{} {}".format(graph_label,html.escape(self._name))
            g.attr(label=graph_label, labelloc="t", nodesep="0", ranksep="0")

        if config["type"] == "basic":
            return self._visualizeBasic(g, **config)
        elif config["type"] == "advanced":
            return self._visualizeAdvanced(g, **config)
        else:
            if config["containerGraph"]:
                return ([],{})
            else:
                return g
    
    def _visualizeBasic(self, g, *args, **kwargs):

        def connectArm(arm_key,arm, **kwargs):
            node_key = self._visualize_node_key(config["prefix"],arm.connection())
            if arm.armType()=="incoming":
                g.edge(arm_key,node_key,style="dashed", 
                                       color="grey", rankdir="lr", constraint="true")
            elif arm.armType()=="outgoing":
                g.edge(node_key,arm_key,style="dashed", 
                                       color="grey", rankdir="lr", constraint="true")   
         
        #initialise configuration
        basicConfig = {
            "showAllNodes": False,
            "nodeNumberFontSize": 10,
            "nodePenWidthSelected": 2,
            "edgeDistanceFontSize": 10,
            "edgeNumberFontSize": 8,
            "edgeArrowSize": 0.5
        }
        #construct config
        config = kwargs.copy()
        for key, value in basicConfig.items():
            if not key in config.keys():
                config[key] = basicConfig[key]

        g.attr(rankdir="TB", nodesep="0.3")

        #remember edges and nodes
        edges = set()
        nodes = {}

        if config["restrictedListOfOrientatedCkmers"]:
            for orientatedCkmer in self._orientatedCkmers:
                if not orientatedCkmer in config["restrictedListOfOrientatedCkmers"]:
                    continue
                elif not (config["showAllNodes"] or self._orientatedCkmers[orientatedCkmer].candidate()):
                    continue
                else:
                    #define layout
                    node_style="filled"
                    node_color=config["nodeColorDefault"]
                    node_penwidth=config["nodePenWidthDefault"]
                    if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                        node_fillcolor=config["nodeFillColorSelection"]
                    else:
                        node_fillcolor=config["nodeFillColorDefault"]
                    if self._orientatedCkmers[orientatedCkmer].candidate():
                        if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                            node_fillcolor=config["nodeFillColorSelectionCandidate"]
                        else:
                            node_fillcolor=config["nodeFillColorCandidate"]
                        node_color=config["nodeColorCandidate"]
                    if orientatedCkmer in self._start or orientatedCkmer in self._end:
                        if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                            node_fillcolor=config["nodeFillColorSelectionStartEnd"]
                        else:
                            node_fillcolor=config["nodeFillColorStartEnd"]
                        node_color=config["nodeColorStartEnd"]
                    if orientatedCkmer in self._selected:
                        node_color=config["nodeColorSelected"]
                        node_penwidth=config["nodePenWidthSelected"]
                    node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                    node_label = self._visualize_basic_node_label(orientatedCkmer,
                                                                    {"number": config["nodeNumberFontSize"]}, config)
                    nodes[orientatedCkmer] = {"node": node_key}
                    g.node(node_key, shape="circle",
                           label=node_label, color=node_color, style=node_style, 
                           fillcolor=node_fillcolor, penwidth=str(node_penwidth),
                           width="0.2", margin="0.01")

            for orientatedCkmer1 in self._orientatedCkmers:
                if not orientatedCkmer1 in config["restrictedListOfOrientatedCkmers"]:
                    continue
                elif not (config["showAllNodes"] or self._orientatedCkmers[orientatedCkmer1].candidate()):
                    continue
                else:
                    node_key1 = self._visualize_node_key(config["prefix"],orientatedCkmer1)
                    for orientatedCkmer2 in self._orientatedCkmers[orientatedCkmer1]._outgoing:
                        if not orientatedCkmer2 in config["restrictedListOfOrientatedCkmers"]:
                            continue
                        elif not (config["showAllNodes"] or self._orientatedCkmers[orientatedCkmer2].candidate()):
                            continue
                        else:
                            node_key2 = self._visualize_node_key(config["prefix"],orientatedCkmer2)
                            #only pass once
                            edgeKey = (node_key1, node_key2)
                            if edgeKey in edges:
                                continue
                            else:
                                edges.add(edgeKey)
                            props = self._orientatedCkmers[orientatedCkmer1]._outgoing[orientatedCkmer2]
                            edge_label = self._visualize_edge_label(props, {"distance": config["edgeDistanceFontSize"],
                                                                            "number": config["edgeNumberFontSize"]}, config)
                            # edge_label = "{}".format(props["distance"])
                            edge_direction = "forward"
                            edge_constraint="true"
                            edge_color = config["edgeColorDefault"]
                            edge_penwidth = config["edgePenWidthDefault"]
                            edge_style = config["edgeStyleDefault"]
                            edge_arrowsize = config["edgeArrowSize"]
                            #detect auto generated edges
                            if props["number"]==0:
                                edge_color = config["edgeColorAuto"]
                                edge_penwidth = config["edgePenWidthAuto"]
                                edge_style = config["edgeStyleAuto"]
                            if props["path"] is None:
                                edge_color = config["edgeColorMissingPath"]
                                edge_style = config["edgeStyleMissingPath"]
                            #create the edge
                            g.edge(node_key1, node_key2, label=edge_label, 
                                   constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth)) 
        else:            
            #start with connected candidates
            connectedCandidates = self.getConnectedCandidates()
            for i in range(len(connectedCandidates)):
                #add nodes
                g.node("{}_connected_{}_start".format(config["prefix"],i), shape="point", style="invis")
                for orientatedCkmer in connectedCandidates[i]["connected"]:
                    #define layout
                    node_style="filled"
                    node_fillcolor=config["nodeFillColorCandidate"]
                    node_color=config["nodeColorCandidate"]
                    node_penwidth=config["nodePenWidthDefault"]
                    if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                        node_fillcolor=config["nodeFillColorSelectionCandidate"]
                    else:
                        node_fillcolor=config["nodeFillColorCandidate"]
                    if orientatedCkmer in self._start or orientatedCkmer in self._end:
                        if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                            node_fillcolor=config["nodeFillColorSelectionStartEnd"]
                        else:
                            node_fillcolor=config["nodeFillColorStartEnd"]
                        node_color=config["nodeColorStartEnd"]
                    if orientatedCkmer in self._selected:
                        node_color=config["nodeColorSelected"]
                        node_penwidth=config["nodePenWidthSelected"]
                    node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                    node_label = self._visualize_basic_node_label(orientatedCkmer,
                                                                    {"number": config["nodeNumberFontSize"]}, config)
                    nodes[orientatedCkmer] = {"node": node_key}
                    g.node(node_key, shape="circle",
                           label=node_label, color=node_color, style=node_style, 
                           fillcolor=node_fillcolor, penwidth=str(node_penwidth),
                           width="0.2", margin="0.01")
                g.node("{}_connected_{}_end".format(config["prefix"],i), shape="point", style="invis")    
    
                #add connections
                for orientatedCkmer in connectedCandidates[i]["start"]:
                    if orientatedCkmer in connectedCandidates[i]["connected"]:
                        node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                        g.edge("{}_connected_{}_start".format(config["prefix"],i), node_key,
                               weight="1", style="invis")
                for orientatedCkmer1 in connectedCandidates[i]["connected"]:
                    node_key1 = self._visualize_node_key(config["prefix"],orientatedCkmer1)
                    for orientatedCkmer2 in self._orientatedCkmers[orientatedCkmer1]._outgoing:
                        if orientatedCkmer2 in connectedCandidates[i]["connected"]:
                            node_key2 = self._visualize_node_key(config["prefix"],orientatedCkmer2)
                            props = self._orientatedCkmers[orientatedCkmer1]._outgoing[orientatedCkmer2]
                            #only pass once
                            edgeKey = (node_key1, node_key2)
                            if edgeKey in edges:
                                continue
                            else:
                                edges.add(edgeKey)
                            props = self._orientatedCkmers[orientatedCkmer1]._outgoing[orientatedCkmer2]
                            edge_label = self._visualize_edge_label(props, {"distance": config["edgeDistanceFontSize"],
                                                                            "number": config["edgeNumberFontSize"]}, config)
                            # edge_label = "{}".format(props["distance"])
                            edge_direction = "forward"
                            edge_constraint="true"
                            edge_color = config["edgeColorDefault"]
                            edge_penwidth = config["edgePenWidthDefault"]
                            edge_style = config["edgeStyleDefault"]
                            edge_arrowsize = config["edgeArrowSize"]
                            #detect auto generated edges
                            if props["number"]==0:
                                edge_color = config["edgeColorAuto"]
                                edge_penwidth = config["edgePenWidthAuto"]
                                edge_style = config["edgeStyleAuto"]
                            if props["path"] is None:
                                edge_color = config["edgeColorMissingPath"]
                                edge_style = config["edgeStyleMissingPath"]
                            #create the edge
                            g.edge(node_key1, node_key2, label=edge_label, 
                                   constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth)) 
                for orientatedCkmer in connectedCandidates[i]["end"]:
                    if orientatedCkmer in connectedCandidates[i]["connected"]:
                        node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                        g.edge(node_key,"{}_connected_{}_end".format(config["prefix"],i),weight="1",style="invis")
    
            if config["showAllNodes"]:
                #add all other nodes
                for orientatedCkmer in self._orientatedCkmers:
                    if not orientatedCkmer in nodes:
                        #define layout
                        node_style="filled"
                        if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                            node_fillcolor=config["nodeFillColorSelection"]
                        else:
                            node_fillcolor=config["nodeFillColorDefault"]
                        node_color=config["nodeColorDefault"]
                        node_penwidth=config["nodePenWidthDefault"]
                        node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                        node_label = self._visualize_basic_node_label(orientatedCkmer,
                                                                        {"number": config["nodeNumberFontSize"]}, config)
                        nodes[orientatedCkmer] = {"node": node_key}
                        g.node(node_key, shape="circle",
                               label=node_label, color=node_color, style=node_style, 
                               fillcolor=node_fillcolor, penwidth=str(node_penwidth),
                               width="0.2", margin="0.01")
                #add all missing connections
                for orientatedCkmer1 in self._orientatedCkmers:
                    node_key1 = self._visualize_node_key(config["prefix"],orientatedCkmer1)
                    for orientatedCkmer2 in self._orientatedCkmers[orientatedCkmer1]._outgoing:
                        node_key2 = self._visualize_node_key(config["prefix"],orientatedCkmer2)
                        #only pass once
                        edgeKey = (node_key1, node_key2)
                        if edgeKey in edges:
                            continue
                        else:
                            edges.add(edgeKey)
                        props = self._orientatedCkmers[orientatedCkmer1]._outgoing[orientatedCkmer2]
                        edge_label = self._visualize_edge_label(props, {"distance": config["edgeDistanceFontSize"],
                                                                        "number": config["edgeNumberFontSize"]}, config)
                        edge_direction  = "forward"
                        edge_constraint = "true"
                        edge_color = config["edgeColorDefault"]
                        edge_penwidth = config["edgePenWidthDefault"]
                        edge_style = config["edgeStyleDefault"]
                        edge_arrowsize = config["edgeArrowSize"]
                        #detect auto generated edges
                        if props["number"]==0:
                            edge_color = config["edgeColorAuto"]
                            edge_penwidth = config["edgePenWidthAuto"]
                            edge_style = config["edgeStyleAuto"]
                        if props["path"] is None:
                            edge_color = config["edgeColorMissingPath"]
                            edge_style = config["edgeStyleMissingPath"]
                        #create the edge
                        g.edge(node_key1, node_key2, label=edge_label, 
                               constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                               dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth)) 
        
        #dead ends
        if ((config["showAllNodes"] or config["showDeadEnds"])):
            edge_arrowsize = config["edgeArrowSize"]
            for orientatedCkmer in nodes:
                if (not orientatedCkmer in config["hideDeadEndBefore"]
                    and not self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd is None):
                    node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                    dead_key = "{}_{}_{}_{}".format(config["prefix"],
                                orientatedCkmer[0],orientatedCkmer[1],
                                self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[0])
                    g.node(dead_key, shape="point")
                    edge_label = self._visualize_dead_end_label(
                          self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[1], {"distance": config["edgeDistanceFontSize"],
                                                                    "number": config["edgeNumberFontSize"]}, config)
                    edge_direction  = "forward"
                    edge_constraint = "true"
                    edge_color = config["edgeColorAuto"]
                    edge_penwidth = config["edgePenWidthAuto"]
                    edge_style = config["edgeStyleAuto"]
                    edge_arrowsize = config["edgeArrowSize"]
                    g.edge(dead_key, node_key, label=edge_label, 
                           constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                           dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth))
                if (not orientatedCkmer in config["hideDeadEndAfter"]
                    and not self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd is None):
                    node_key = self._visualize_node_key(config["prefix"],orientatedCkmer)
                    dead_key = "{}_{}_{}_{}".format(config["prefix"],
                        orientatedCkmer[0],orientatedCkmer[1],
                        self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[0])
                    g.node(dead_key, shape="point")
                    edge_label = self._visualize_dead_end_label(
                          self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[1], {"distance": config["edgeDistanceFontSize"],
                                                                    "number": config["edgeNumberFontSize"]}, config)
                    edge_direction  = "forward"
                    edge_constraint = "true"
                    edge_color = config["edgeColorAuto"]
                    edge_penwidth = config["edgePenWidthAuto"]
                    edge_style = config["edgeStyleAuto"]
                    edge_arrowsize = config["edgeArrowSize"]
                    g.edge(node_key, dead_key, label=edge_label, 
                           constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                           dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth))

        arms = self.getArms()
        if config["showAllNodes"]:
            #warn if relevant
            if config["showArms"] and len(arms)>0:
                self._logger.warning("can't show {} detected arms if all nodes are shown".format(len(arms)))
            if config["showPotentialConnectedArms"]:
                connectedCandidateArms = self._detectConnectedArmsCandidates()
                if len(connectedCandidateArms)>0:
                    self._logger.warning("can't show {} detected connected arms candidates if all nodes are shown".format(
                        len(connectedCandidateArms)))
        else:
            processedArms = set()
            #show arms
            if config["showArms"]:
                #show potential connected arms
                if config["showPotentialConnectedArms"]:
                    connectedCandidateArms = self._detectConnectedArmsCandidates()
                    for i in range(len(connectedCandidateArms)):
                        arm1 = self.getArm(connectedCandidateArms[i][0])
                        arm2 = self.getArm(connectedCandidateArms[i][1])
                        if arm1.connection() in nodes and arm2.connection() in nodes:
                            processedArms.add(arm1.id())
                            processedArms.add(arm2.id())
                            arm_key1, arm_key2 = self._visualizeConnectedArms(g,i,arm1,arm2,**config)
                            connectArm(arm_key1,arm1,**config)
                            connectArm(arm_key2,arm2,**config)
                #loop over arms
                for i in range(len(arms)):
                    arm = arms[i]
                    #only if no potential connected arms
                    if arm.id() in processedArms:
                        continue
                    else:
                        processedArms.add(arm.id())
                    #only if visible
                    if arm.connection() in nodes:
                        #create arm
                        arm_key = self._visualizeArm(g,arm,**config)
                        connectArm(arm_key,arm,**config)
            #show potential connected arms
            elif config["showPotentialConnectedArms"]:
                connectedCandidateArms = self._detectConnectedArmsCandidates()
                for i in range(len(connectedCandidateArms)):
                    arm1 = self.getArm(connectedCandidateArms[i][0])
                    arm2 = self.getArm(connectedCandidateArms[i][1])
                    if arm1.connection() in nodes and arm2.connection() in nodes:
                        node_key1 = self._visualize_node_key(config["prefix"],arm1.connection())
                        node_key2 = self._visualize_node_key(config["prefix"],arm2.connection())
                        self._visualizeConnectedArmsConnection(g, node_key1,node_key2, **config)
        if config["containerGraph"]:
            return (list(edges), nodes)
        else:
            return g


    def _visualizeAdvanced(self, g, *args, **kwargs):
        def connectArm(arm_key, arm, **kwargs):
            orientatedCkmersArm = arm._orientatedCkmers.intersection(orientatedCkmerNodes.keys())
            if len(orientatedCkmersArm)==0:
                for node_key in orientatedCkmerNodes[arm.connection()].values():
                    if arm.armType()=="incoming":
                        g.edge(arm_key,node_key,style="dashed", 
                                               color="grey", rankdir="lr", constraint="true")
                    elif arm.armType()=="outgoing":
                        g.edge(node_key,arm_key,style="dashed", 
                                               color="grey", rankdir="lr", constraint="true")
            else:
                for orientatedCkmer in orientatedCkmersArm:
                    for node_key in orientatedCkmerNodes[orientatedCkmer].values():
                        if arm.armType()=="incoming":
                            if(len(set(self._orientatedCkmers[orientatedCkmer]._incoming.keys())
                                   .intersection(orientatedCkmerNodes))>0):
                                break
                            g.edge(arm_key,node_key,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                        elif arm.armType()=="outgoing":
                            if(len(set(self._orientatedCkmers[orientatedCkmer]._outgoing.keys())
                                   .intersection(orientatedCkmerNodes))>0):
                                break
                            g.edge(node_key,arm_key,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                        
        #initialise configuration
        advancedConfig = {
            "showAllBases": False,
            "showPrePostBases": True,
            "showPrePostBasesMaxSteps": 2,
            "baseFontSize": 6,
            "baseNumberFontSize": 6,
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
            "nodeNumberFontSize": 8,
            "nodeLetterFontSize": 14,
            "edgeDistanceFontSize": "14",
            "edgeNumberFontSize": "8"
        }
        #construct config
        config = kwargs.copy()
        for key, value in advancedConfig.items():
            if not key in config.keys():
                config[key] = advancedConfig[key]

        #define base containers and bases
        orientatedBaseContainers = {}
        orientatedCkmerNodes = {}
        edges = set()
        
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
                                          if self._orientatedBases[orientatedBase]._order is None 
                                          else self._orientatedBases[orientatedBase]._order)
        
        #restrict bases
        if config["restrictedListOfOrientatedCkmers"]:
            restrictedListOfBases = set()
            for orientatedCkmer in config["restrictedListOfOrientatedCkmers"]:
                restrictedListOfBases.update(self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
            sortedOrientatedBaseList = [base for base in sortedOrientatedBaseList if base in restrictedListOfBases]
        
        #loop over sorted orientated bases
        for orientatedBase in sortedOrientatedBaseList:
            
            #optionally, hide bases
            if not config["showAllBases"]:
                if (config["showPrePostBases"] and not (self._orientatedBases[orientatedBase]._preStart is None and 
                                                        self._orientatedBases[orientatedBase]._postEnd is None)):
                    if config["showPrePostBasesMaxSteps"] and isinstance(config["showPrePostBasesMaxSteps"],int):
                        if not (self._orientatedBases[orientatedBase]._preStart is None or 
                                self._orientatedBases[orientatedBase]._preStart<=config["showPrePostBasesMaxSteps"]):
                            continue
                        elif not (self._orientatedBases[orientatedBase]._postEnd is None or 
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
                    obc.attr(style="invis", nodesep="0", ranksep="0")
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
                    elif not (self._orientatedBases[orientatedBase]._preStart is None and 
                          self._orientatedBases[orientatedBase]._postEnd is None):
                        base_style=config["baseStylePrePost"]
                        base_fillcolor=config["baseFillColorPrePost"]
                        base_color=config["baseColorPrePost"]
                        base_penwidth=config["basePenWidthPrePost"]
                    else:
                        base_style=config["baseStyleOrder"]
                        base_fillcolor=config["baseFillColorOrder"]
                        base_color=config["baseColorOrder"]
                        base_penwidth=config["basePenWidthOrder"]
                    base_label = self._visualize_base_label(orientatedBase, {"base": config["baseFontSize"],
                                                                             "number": config["baseNumberFontSize"]}, config) 
                    c.attr(style=base_style, color=base_color, 
                           fillcolor=base_fillcolor, label=base_label, penwidth=str(base_penwidth))
                    #add k-mers to base container
                    for orientatedCkmer in self._orientatedBases[orientatedBase]._orientatedCkmers:  
                        #define layout
                        node_style="filled"
                        if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                            node_fillcolor=config["nodeFillColorSelection"]
                        else:
                            node_fillcolor=config["nodeFillColorDefault"]
                        node_color=config["nodeColorDefault"]
                        node_penwidth=config["nodePenWidthDefault"]
                        if orientatedCkmer in self._orientatedCkmers.keys():
                            if self._orientatedCkmers[orientatedCkmer].candidate():
                                if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                                    node_fillcolor=config["nodeFillColorSelectionCandidate"]
                                else:
                                    node_fillcolor=config["nodeFillColorCandidate"]
                                node_color=config["nodeColorCandidate"]
                            elif (not config["showAllBases"]) and config["showArms"]:
                                if self._orientatedCkmers[orientatedCkmer].incomingArmType():
                                    node_fillcolor=config["nodeFillColorIncomingArm"]
                                    node_color=config["nodeColorIncomingArm"]
                                elif self._orientatedCkmers[orientatedCkmer].outgoingArmType():
                                    node_fillcolor=config["nodeFillColorOutgoingArm"]
                                    node_color=config["nodeColorOutgoingArm"]
                        if orientatedCkmer in self._start or orientatedCkmer in self._end:
                            if config["selectionOrientatedCkmers"] and (orientatedCkmer in config["selectionOrientatedCkmers"]):
                                node_fillcolor=config["nodeFillColorSelectionStartEnd"]
                            else:
                                node_fillcolor=config["nodeFillColorStartEnd"]
                            node_color=config["nodeColorStartEnd"]
                        if orientatedCkmer in self._selected:
                            node_color=config["nodeColorSelected"]
                            node_penwidth=config["nodePenWidthSelected"]
                        node_label = self._visualize_advanced_node_label(orientatedBase, orientatedCkmer,
                                                                {"letter": config["nodeLetterFontSize"], 
                                                                 "number": config["nodeNumberFontSize"]}, config)
                        node_key = self._visualize_node_key(config["prefix"],orientatedCkmer, orientatedBase)
                        c.node(node_key, label=node_label, color=node_color, style=node_style, 
                               fillcolor=node_fillcolor, penwidth=str(node_penwidth))
                        if ((config["showAllBases"] or config["showDeadEnds"]) 
                            and not orientatedCkmer in config["hideDeadEndBefore"]
                            and not self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd is None):
                            dead_key = "{}_{}_{}_{}".format(config["prefix"],
                                orientatedBase[0],orientatedBase[1],
                                self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[0])
                            g.node(dead_key, shape="point")
                            g.edge(dead_key, node_key, constraint="true",
                              label=self._visualize_dead_end_label(
                                  self._orientatedCkmers[orientatedCkmer]._incomingDeadEnd[1], {}, config),
                                   dir="forward", color=config["edgeColorAuto"], weight="1",
                                   style=config["edgeStyleAuto"], penwidth=str(config["edgePenWidthAuto"]))
                        if ((config["showAllBases"] or config["showDeadEnds"])
                            and not orientatedCkmer in config["hideDeadEndAfter"]
                            and not self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd is None):
                            dead_key = "{}_{}_{}_{}".format(config["prefix"],
                                orientatedBase[0],orientatedBase[1],
                                self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[0])
                            g.node(dead_key, shape="point")
                            g.edge(node_key, dead_key, constraint="true", 
                              label=self._visualize_dead_end_label(
                                  self._orientatedCkmers[orientatedCkmer]._outgoingDeadEnd[1], {}, config),
                                   dir="forward", color=config["edgeColorAuto"], weight="1",
                                   style=config["edgeStyleAuto"], penwidth=str(config["edgePenWidthAuto"]))
                        #register
                        if orientatedCkmer in orientatedCkmerNodes.keys():
                            assert orientatedBase[1] not in orientatedCkmerNodes[orientatedCkmer]
                            orientatedCkmerNodes[orientatedCkmer]["{}_{}".format(config["prefix"],orientatedBase[1])] = node_key
                            for other_node_key in orientatedCkmerNodes[orientatedCkmer].values():
                                if not node_key==other_node_key:
                                    g.edge(other_node_key, node_key, constraint="false",
                                       dir="none", color="grey", style="dashed", penwidth="2")  
                                    
                        else:
                            orientatedCkmerNodes[orientatedCkmer] = {"{}_{}".format(config["prefix"],orientatedBase[1]): node_key}
                            
        #try to sort orientated k-mers
        maxOrientatedCkmerOrder = (0 if len(self._orientatedBases)==0 
                                  else max([self._orientatedBases[orientatedBase]._order or 0 
                                      for orientatedBase in self._orientatedBases.keys()]))
        sortedOrientatedCkmerList = sorted(self._orientatedCkmers.keys(), 
                                key=lambda orientatedCkmer: maxOrientatedCkmerOrder+1 
                                           if self._orientatedCkmers[orientatedCkmer]._order is None 
                                           else self._orientatedCkmers[orientatedCkmer]._order)

        
        #restrict k-mers
        if config["restrictedListOfOrientatedCkmers"]:
            sortedOrientatedCkmerList = [ckmer for ckmer in sortedOrientatedCkmerList 
                                         if ckmer in config["restrictedListOfOrientatedCkmers"]]
        
        #loop over the sorted orientated k-mers
        forward_key = "{}_{}".format(config["prefix"],"forward")
        backward_key = "{}_{}".format(config["prefix"],"backward")
        for orientatedCkmer in sortedOrientatedCkmerList:
            if orientatedCkmer in orientatedCkmerNodes.keys():                
                kmerSet = self._orientatedCkmers[orientatedCkmer]._incoming.keys()
                for connectedCkmer in kmerSet:
                    if connectedCkmer in orientatedCkmerNodes.keys():                        
                        edge_label = self._visualize_edge_label(
                            self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer], 
                                                    {"distance": config["edgeDistanceFontSize"],
                                                     "number": config["edgeNumberFontSize"]}, config
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
                        if backward_key in orientatedCkmerNodes[orientatedCkmer]:
                            node = orientatedCkmerNodes[orientatedCkmer][backward_key]
                        else:
                            node = orientatedCkmerNodes[orientatedCkmer][forward_key]
                        if forward_key in orientatedCkmerNodes[connectedCkmer]:
                            connectedNode = orientatedCkmerNodes[connectedCkmer][forward_key]
                        else:
                            connectedNode = orientatedCkmerNodes[connectedCkmer][backward_key]
                        if self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["path"] is None:
                            edge_color = config["edgeColorMissingPath"]
                            edge_style = config["edgeStyleMissingPath"]
                        #create the edge only once
                        edgeKey = (connectedNode, node)
                        if not edgeKey in edges:
                            edges.add(edgeKey)
                            g.edge(connectedNode, node, label=edge_label, constraint=edge_constraint, style=edge_style, 
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth))  
                        
        arms = self.getArms()
        if config["showAllBases"]:
            #warn if relevant
            if config["showArms"] and len(arms)>0:
                self._logger.warning("can't show {} detected arms if all bases are shown".format(len(arms)))
            if config["showPotentialConnectedArms"]:
                connectedCandidateArms = self._detectConnectedArmsCandidates()
                if len(connectedCandidateArms)>0:
                    self._logger.warning("can't show {} detected connected arms candidates if all bases are shown".format(
                        len(connectedCandidateArms)))
        else:
            processedArms = set()
            #show arms
            if config["showArms"]:
                #show potential connected arms
                if config["showPotentialConnectedArms"]:
                    connectedCandidateArms = self._detectConnectedArmsCandidates()
                    for i in range(len(connectedCandidateArms)):
                        arm1 = self.getArm(connectedCandidateArms[i][0])
                        arm2 = self.getArm(connectedCandidateArms[i][1])
                        if arm1.connection() in orientatedCkmerNodes and arm2.connection() in orientatedCkmerNodes:
                            processedArms.add(arm1.id())
                            processedArms.add(arm2.id())
                            arm_key1, arm_key2 = self._visualizeConnectedArms(g,i,arm1,arm2,**config)
                            connectArm(arm_key1,arm1,**config)
                            connectArm(arm_key2,arm2,**config)
                #loop over arms
                for i in range(len(arms)):
                    arm = arms[i]
                    #only if no potential connected arms
                    if arm.id() in processedArms:
                        continue
                    else:
                        processedArms.add(arm.id())
                    #only if visible
                    if arm.connection() in orientatedCkmerNodes:
                        #create arm
                        arm_key = self._visualizeArm(g,arm,**config)
                        connectArm(arm_key,arm,**config)
            #show potential connected arms
            elif config["showPotentialConnectedArms"]:
                connectedCandidateArms = self._detectConnectedArmsCandidates()
                for i in range(len(connectedCandidateArms)):
                    arm1 = self.getArm(connectedCandidateArms[i][0])
                    arm2 = self.getArm(connectedCandidateArms[i][1])
                    if arm1.connection() in orientatedCkmerNodes and arm2.connection() in orientatedCkmerNodes:
                        orientatedCkmersArm1 = arm1._orientatedCkmers.intersection(orientatedCkmerNodes.keys())
                        orientatedCkmersArm2 = arm2._orientatedCkmers.intersection(orientatedCkmerNodes.keys())
                        if len(orientatedCkmersArm1)==0:
                            orientatedCkmersArm1.add(arm1.connection())
                        if len(orientatedCkmersArm2)==0:
                            orientatedCkmersArm1.add(arm2.connection())
                        for orientatedCkmer1 in orientatedCkmersArm1:
                            for node_key1 in orientatedCkmerNodes[orientatedCkmer1].values():
                                for orientatedCkmer2 in orientatedCkmersArm2:
                                    for node_key2 in orientatedCkmerNodes[orientatedCkmer2].values():
                                        self._visualizeConnectedArmsConnection(g,node_key1,node_key2, **config)
        if config["containerGraph"]:
            return (list(edges), orientatedCkmerNodes)
        else:
            return g
    
    #arm
    def _visualizeArm(self, g, arm, *args, **kwargs):
        initConfig = {
                "type": "basic",
                "prefix": "arm",
                "armIncomingStyle": "filled",
                "armIncomingFillColor": "lightblue",
                "armIncomingPenWidth": 1,
                "armIncomingColor": "black",
                "armOutgoingStyle": "filled",
                "armOutgoingFillColor": "lightcoral",
                "armOutgoingPenWidth": 1,
                "armOutgoingColor": "black"
            }
        config = kwargs.copy()
        for key, value in initConfig.items():
            if not key in config.keys():
                config[key] = value

        armPrefix = "arm_{}_{}".format(config["prefix"],arm.id())
                    
        arm_name = "Arm {}".format(arm.id())
        if arm.armType()=="incoming":
            arm_name = "Incoming {}".format(arm_name)
            arm_style = config["armIncomingStyle"]
            arm_fillcolor = config["armIncomingFillColor"]
            arm_penwidth = config["armIncomingPenWidth"]
            arm_color = config["armIncomingColor"]
        elif arm.armType()=="outgoing":
            arm_name = "Outgoing {}".format(arm_name)
            arm_style = config["armOutgoingStyle"]
            arm_fillcolor = config["armOutgoingFillColor"]
            arm_penwidth = config["armOutgoingPenWidth"]
            arm_color = config["armOutgoingColor"]
        else:
            return
        armGraph=g.subgraph(name="cluster_graph_{}".format(armPrefix))
        with armGraph as ag:
            if config["type"]=="basic":
                arm_label = ("<" + 
                             "<font point-size=\"12\" color=\"grey\">{}</font><br/>".format(arm_name) + 
                             "<font point-size=\"10\">{} nodes</font><br/>".format(arm.n()) + 
                             "<font point-size=\"8\">max {}x</font>".format(arm.maxFreq()) + 
                             ">")
            
            else:
                arm_label = ("<" + 
                         "<font point-size=\"12\" color=\"grey\">{}</font><br/>".format(arm_name) + 
                         "<font point-size=\"10\">size {} with {} nodes</font><br/>".format(
                             arm.size(),arm.n()) + 
                         "<font point-size=\"8\">maximum frequency: {}x</font><br/>".format(
                             arm.maxFreq()) + 
                         ">")
            ag.attr(label=arm_label, style=arm_style, fillcolor=arm_fillcolor, 
                color=arm_color, penwidth=str(arm_penwidth), labelloc="t", nodesep="0", ranksep="0")

            arm_key = "arm_{}_{}".format(armPrefix,arm.armType())
            ag.node(arm_key, shape="point")
        return arm_key

    def _visualizeConnectedArms(self, g, id, arm1, arm2, **kwargs):
        initConfig = {
            "prefix": "graph",
            "connectedArmsStyle": "dashed,filled",
            "connectedArmsFillColor": "orange",
            "connectedArmsPenWidth": 1,
            "connectedArmsColor": "black"
        }
        config = kwargs.copy()
        for key,value in initConfig.items():
            if not key in config.keys():
                config[key] = value
        connectedArmsGraph=g.subgraph(name="cluster_{}_connectedarms_{}_{}_{}".format(config["prefix"],id,arm1.id(),arm2.id()))
        with connectedArmsGraph as tg:
            connectedarms_label = "Potentially Connected"
            connectedarms_style = config["connectedArmsStyle"]
            connectedarms_fillcolor = config["connectedArmsFillColor"]
            connectedarms_penwidth = config["connectedArmsPenWidth"]
            connectedarms_color = config["connectedArmsColor"]
            tg.attr(label=connectedarms_label, style=connectedarms_style, fillcolor=connectedarms_fillcolor, 
                    color=connectedarms_color, penwidth=str(connectedarms_penwidth), 
                    constraint="false", rankdir="TB", labelloc="t", nodesep="0", ranksep="0")
            #arms
            arm_key1 = self._visualizeArm(tg,arm1,**config)
            arm_key2 = self._visualizeArm(tg,arm2,**config)
            tg.edge(arm_key1,arm_key2,weight="1",style="invis")
        return arm_key1,arm_key2

    def _visualizeConnectedArmsConnection(self, g, arm_key1, arm_key2, **kwargs):
        initConfig = {
            "connectedArmsEdgeFontColor": "blue",
            "connectedArmsEdgeStyle": "dotted",
            "connectedArmsEdgePenwidth": 2,
            "connectedArmsEdgeFontSize": 12,
            "connectedArmsEdgeColor": "blue"
        }
        config = kwargs.copy()
        for key,value in initConfig.items():
            if not key in config.keys():
                config[key] = value
        edge_label = ("<" + 
                     "<font point-size=\"{}\" color=\"{}\">possible<br/> connected arms</font>".format(
                         config["connectedArmsEdgeFontSize"], config["connectedArmsEdgeColor"]
                     ) + 
                     ">")
        edge_style = config["connectedArmsEdgeStyle"]
        edge_color = config["connectedArmsEdgeColor"]
        edge_penwidth = config["connectedArmsEdgePenwidth"]
        g.edge(arm_key1, arm_key2,style=edge_style, 
                   label=edge_label, color=edge_color, rankdir="lr", 
                   constraint="true", penwidth=str(edge_penwidth))
        
    #base label definition
    def _visualize_base_label(self, orientatedBase: tuple, fontSize:dict = {}, config:dict = {}):
        base_label = "<"
        base_label += "<font point-size=\"{}\">".format(fontSize.get("base","0"))
        if orientatedBase[1]=="forward":
            base_label += orientatedBase[0]+"*"
        elif orientatedBase[1]=="backward":
            base_label += "*"+General.reverse_complement(orientatedBase[0])
        else:
            base_label += "---orientation unknown---"
        base_label +="<br/><font point-size=\"{}\" color=\"grey\">{}x</font></font>".format(
            fontSize.get("number","0"),self._orientatedBases[orientatedBase]._number)
        base_label +=">"
        return base_label
    
    #edge label definition
    def _visualize_edge_label(self, info: dict, fontSize:dict = {}, config:dict = {}):
        edge_label = "<"
        edge_label += "<font point-size=\"{}\" color=\"grey\">{}</font><br/>".format(fontSize.get("distance","0"), info["distance"])
        edge_label += "<font point-size=\"{}\">{}x</font>".format(fontSize.get("number","0"), info["number"])
        edge_label += ">"
        return edge_label
    
    def _visualize_dead_end_label(self, distance: int, fontSize:dict = {}, config:dict = {}):
        edge_label = "<"
        edge_label += "<font point-size=\"{}\" color=\"grey\">{}</font>".format(fontSize.get("distance","0"), distance)
        edge_label += ">"
        return edge_label
    
    #node key definition
    def _visualize_node_key(self, prefix, orientatedCkmer: tuple, orientatedBase: tuple = None):
        if orientatedBase is None:
            return "{}_{}_{}".format(prefix,
                                    orientatedCkmer[0],orientatedCkmer[1])
        else:
            return "{}_{}_{}_{}_{}".format(prefix,
                                    orientatedBase[0],orientatedBase[1],
                                    orientatedCkmer[0],orientatedCkmer[1])
    
    #node label definition
    def _visualize_basic_node_label(self, orientatedCkmer: tuple, fontSize:dict = {}, config:dict = {}):
        node_label = "<" + "<font point-size=\"{}\">{}x</font>".format(
            fontSize.get("number","0"),self._orientatedCkmers[orientatedCkmer]._number)
        node_label += ">"
        return node_label
    
    def _visualize_advanced_node_label(self, orientatedBase: tuple, orientatedCkmer: tuple, fontSize:dict = {}, config:dict = {}):
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
        node_label = "<" + "<font point-size=\"{}\" color=\"blue\">{}</font>".format(fontSize.get("letter","0"),node_letter)
        node_label +="<font point-size=\"{}\">{}x</font>".format(
            fontSize.get("number","0"),self._orientatedCkmers[orientatedCkmer]._number)
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
        connected = self.getConnected()
        connectedSets = []
        #compute connected candidates        
        candidates = self.getCandidates()
        for k in candidates:    
            if baseConnected:
                kmerSet = set([k])
                for orientatedBase in self._orientatedCkmers[k]._orientatedBases.values():
                    kmerSet.update(self._orientatedBases[orientatedBase]._orientatedCkmers)
                c = set()
                for ks in kmerSet:
                    c.update([m for m in candidates if m==ks or connected[ks][m] or connected[m][ks]])
            else:
                c = set([m for m in candidates if m==k or connected[k][m] or connected[m][k]]) 
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
                elif sum(connected[connectedList][k])==0:
                    if k in self._end:
                        pass
                    else:
                        self._logger.debug("possible new end entry")
                    entry["end"].add(k)
                if k in self._end:
                    entry["end"].add(k)
                elif sum(connected[k][connectedList])==0:
                    if k in self._start:
                        pass
                    else:
                        self._logger.debug("possible new start entry")
                    entry["start"].add(k)
            result.append(entry)
        return result
        
    """get direct distances between orientated k-mers"""
    def getDirectDistances(self):
        if self._directDistances is None:
            self._computeDirectDistances()
        else:
            list1 = set([c for c in self._orientatedCkmers])
            list2 = set(self._directDistances.index)
            if not list1==list2:
                self._computeDirectDistances()    
        return self._directDistances.copy()
    
    def _computeDirectDistances(self):
        #initialise
        orientatedCkmerList = [c for c in self._orientatedCkmers]
        self._directDistances = pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0)
        self._logger.debug("compute direct distances for the De Bruijn Graph")  
        for orientatedCkmer in orientatedCkmerList:
            for connectedCkmer in self._orientatedCkmers[orientatedCkmer]._incoming:
                if connectedCkmer in orientatedCkmerList:
                    d = self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["distance"]
                    self._directDistances.at[orientatedCkmer,connectedCkmer] = -d
                    self._directDistances.at[connectedCkmer,orientatedCkmer] = d
        
    """get frequencies read evidence of direct connections between orientated k-mers"""
    def getDirectNumbers(self):
        if self._directNumbers is None:
            self._computeDirectNumbers()
        else:
            list1 = set([c for c in self._orientatedCkmers])
            list2 = set(self._directNumbers.index)
            if not list1==list2:
                self._computeDirectNumbers()    
        return self._directNumbers.copy()
                           
    def _computeDirectNumbers(self):
        #initialise
        orientatedCkmerList = [c for c in self._orientatedCkmers]
        self._directNumbers = pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0)
        self._logger.debug("compute direct numbers for the De Bruijn Graph")  
        for orientatedCkmer in orientatedCkmerList:
            for connectedCkmer in self._orientatedCkmers[orientatedCkmer]._incoming:
                if connectedCkmer in orientatedCkmerList:
                    n = self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["number"]
                    self._directNumbers.at[orientatedCkmer,connectedCkmer] = n
                    self._directNumbers.at[connectedCkmer,orientatedCkmer] = n

    """get connected orientated k-mers"""
    def getConnected(self):
        if self._connected is None:
            self._computeConnected()
        else:
            list1 = set([c for c in self._orientatedCkmers])
            list2 = set(self._connected.index)
            if not list1==list2:
                self._computeConnected()
        return self._connected.copy()
    
    """compute connected orientated k-mers"""
    def _computeConnected(self):
        orientatedCkmerList = [c for c in self._orientatedCkmers]
        connected = np.identity(len(orientatedCkmerList), dtype="bool")
        outgoing = np.zeros((len(orientatedCkmerList),len(orientatedCkmerList)), dtype = "bool") 
        for i in range(len(orientatedCkmerList)):
            cKmer = self._orientatedCkmers[orientatedCkmerList[i]]
            if len(cKmer._outgoing)>0:
                outgoing[i,[orientatedCkmerList.index(c) 
                            for c in cKmer._outgoing 
                            if c in orientatedCkmerList]] = 1            
        G = nx.from_numpy_array(outgoing, create_using=nx.DiGraph)
        for i in range(len(outgoing)):
            connected[i,list(nx.dfs_tree(G, source=i))] = True
        self._connected = pd.DataFrame(connected, index=orientatedCkmerList, columns=orientatedCkmerList)
        
    
#     """get shortest distances between orientated k-mers"""
#     def getShortestDistances(self):
#         if self._distances is None:
#             self._computeDistances()
#         return self._distances["shortest"]
        
#     def getShortestDistance(self, fromOrientatedKmers: set, toOrientatedKmers: set):
#         """
#         get the shortest distance between two sets of orientated k-mers
#         """
#         if self._distances is None:
#             self._computeDistances()
#         distances = []
#         for c1 in fromOrientatedKmers:
#             for c2 in toOrientatedKmers:
#                 if self._distances["shortest"].at[c1,c2]!=0:
#                     distances.append(self._distances["shortest"].at[c1,c2])
#         if len(distances)==0:
#             return 0
#         else:
#             return(sorted(distances)[0])
        
    def getShortestDirectDistance(self, fromOrientatedKmers: set, toOrientatedKmers: set):
        """
        get the shortest direct distance between two sets of orientated k-mers
        """
        distances = []
        for c1 in fromOrientatedKmers:
            if c1 in self._orientatedCkmers:
                for c2 in toOrientatedKmers.intersection(self._orientatedCkmers[c1]._outgoing.keys()):
                    distances.append(self._orientatedCkmers[c1]._outgoing[c2]["distance"])
        if len(distances)==0:
            return 0
        else:
            return(sorted(distances)[0])
        
    def _detectArms(self):
        """
        detect arms
        """
        if len(self._arms)==0:
            candidates = list(self.getCandidates())
            candidateBases = list(self.getCandidateBases())
            dist = self.getDirectDistances()
            connected = self.getConnected()
            other = set(self._orientatedCkmers.keys()).difference(candidates)
            incomingArms = {}
            outgoingArms = {}
            for orientatedCkmer in other:
                entryDistances = dist[orientatedCkmer][candidates]
                incomingDistances = [d for d in entryDistances if d<0]
                outgoingDistances = [d for d in entryDistances if d>0]
                if len(incomingDistances)>0:
                    d = max(incomingDistances)
                    #compute the candidate k-mers where the arm is incoming
                    for candidateCkmer,candidateDistance in entryDistances.items():
                        assert self._orientatedCkmers[candidateCkmer].candidate()
                        if candidateDistance==d:
                            if candidateCkmer in incomingArms:
                                if not orientatedCkmer in incomingArms[candidateCkmer]["orientatedCkmers"]:
                                    incomingArms[candidateCkmer]["n"]+=1
                                    incomingArms[candidateCkmer]["number"]=max(incomingArms[candidateCkmer]["number"],
                                        self._orientatedCkmers[orientatedCkmer]._number)
                                    incomingArms[candidateCkmer]["orientatedCkmers"].add(orientatedCkmer)
                                    incomingArms[candidateCkmer]["orientatedBases"].update(
                                        self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
                                incomingArms[candidateCkmer]["size"] = max(incomingArms[candidateCkmer]["size"],d)
                            else:
                                incomingArms[candidateCkmer] = {"size": d, "n": 0, "number": 0,
                                        "orientatedBases": set(), "orientatedCkmers": set()}
                                for entry in other.intersection(connected.index[connected[orientatedCkmer]]):
                                    incomingArms[candidateCkmer]["n"]+=1
                                    incomingArms[candidateCkmer]["number"]=max(incomingArms[candidateCkmer]["number"],
                                        self._orientatedCkmers[entry]._number)
                                    incomingArms[candidateCkmer]["orientatedCkmers"].add(entry)
                                    incomingArms[candidateCkmer]["orientatedBases"].update(
                                        self._orientatedCkmers[entry]._orientatedBases.values())
                if len(outgoingDistances)>0:
                    d = min(outgoingDistances)
                    for candidateCkmer,candidateDistance in entryDistances.items():                        
                        assert self._orientatedCkmers[candidateCkmer].candidate()
                        if candidateDistance==d:
                            if candidateCkmer in outgoingArms:
                                if not orientatedCkmer in outgoingArms[candidateCkmer]["orientatedCkmers"]:
                                    outgoingArms[candidateCkmer]["n"]+=1
                                    outgoingArms[candidateCkmer]["number"]=max(outgoingArms[candidateCkmer]["number"],
                                        self._orientatedCkmers[orientatedCkmer]._number)
                                    outgoingArms[candidateCkmer]["orientatedCkmers"].add(orientatedCkmer)
                                    outgoingArms[candidateCkmer]["orientatedBases"].update(
                                        self._orientatedCkmers[orientatedCkmer]._orientatedBases.values())
                                outgoingArms[candidateCkmer]["size"] = min(outgoingArms[candidateCkmer]["size"],d)
                            else:
                                outgoingArms[candidateCkmer] = {"size": d, "n": 0, "number": 0,
                                        "orientatedBases": set(), "orientatedCkmers": set()}
                                for entry in other.intersection(
                                        connected.columns[connected.loc[[orientatedCkmer]].values[0]]):
                                    outgoingArms[candidateCkmer]["n"]+=1
                                    outgoingArms[candidateCkmer]["number"]=max(outgoingArms[candidateCkmer]["number"],
                                        self._orientatedCkmers[entry]._number)
                                    outgoingArms[candidateCkmer]["orientatedCkmers"].add(entry)
                                    outgoingArms[candidateCkmer]["orientatedBases"].update(
                                        self._orientatedCkmers[entry]._orientatedBases.values())
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

    def _resetArms(self):
        for orientatedCkmer in self._orientatedCkmers:
            self._orientatedCkmers[orientatedCkmer]._incomingArm = None
            self._orientatedCkmers[orientatedCkmer]._outgoingArm = None
        for orientatedCkmer in self._orientatedCkmers:
            self._orientatedCkmers[orientatedCkmer]._unsetArm()          
        self._arms = []

    def getArms(self):
        return self._arms

    def connectArms(self, **args):

        numberOfCandidates = args.get("numberOfCandidates",3)
        numberOfKmers = args.get("numberOfKmers",10)
        maximumLength = args.get("maximumLength",1000)
        maximumSteps = args.get("maximumSteps",100)
        maximumPaths = args.get("maximumPaths",0)
        minimumReadDepth = args.get("minimumReadDepth",2)
        minimumGlueSize = args.get("minimumGlueSize",10)
        
        def processArmReads(reads,arm,newArmKmers={}):
            #first raw filtering
            filteredReads = []
            candidateKmers = set([entry[0] for entry in self.getCandidates()])
            armKmers = set([item[0] for item in arm._orientatedCkmers])
            for read in reads:
                #compute type
                kmerTypes = []        
                for kmer in read["kmers"]:
                    if (kmer["ckmer"] in armKmers) or (kmer["ckmer"] in newArmKmers):
                        kmerTypes.append("a")
                    elif kmer["ckmer"] in candidateKmers:
                        kmerTypes.append("c")
                    else:
                        kmerTypes.append("u")
                #filter
                if not "u" in kmerTypes : 
                    #only if no new arm k-mers (initial)
                    if (len(newArmKmers)==0) and not ("a" in kmerTypes and "c" in kmerTypes):
                        continue #nothing new or relevant
                elif not ("a" in kmerTypes or "c" in kmerTypes):
                    continue #no arm or candidate
                elif kmerTypes[0]=="u" and kmerTypes[-1]=="u":
                    continue #doesn't belong to the graph
                elif not (kmerTypes[-1]=="u" or kmerTypes[0]=="u"):
                    continue #problematic
                filteredReads.append(([item.copy() for item in read["kmers"]],kmerTypes,len(read["kmers"]) * [read["number"]],))
            
            #orientate or filter
            orientatedReads = []
            for entry in filteredReads:
                orientation = None
                for i in range(len(entry[0])):
                    if entry[1][i]=="u":
                        continue
                    else:
                        if entry[0][i]["orientation"]=="forward":
                            ockmer = (entry[0][i]["ckmer"],"forward")
                            rckmer = (entry[0][i]["ckmer"],"backward")
                        elif entry[0][i]["orientation"]=="backward":
                            ockmer = (entry[0][i]["ckmer"],"backward")
                            rckmer = (entry[0][i]["ckmer"],"forward")
                        else:
                            orientation = "conflict"
                            continue
                        if ockmer in self.getCandidates() or ockmer in arm._orientatedCkmers:
                            if orientation is None:
                                orientation="forward"
                            elif not orientation=="forward":
                                orientation = "conflict"
                                continue
                        elif rckmer in self.getCandidates() or rckmer in arm._orientatedCkmers:
                            if orientation is None:
                                orientation="backward"
                            elif not orientation=="backward":
                                orientation = "conflict"
                                continue
                        elif ockmer[0] in newArmKmers:
                            if newArmKmers[ockmer[0]]==ockmer[1]:
                                if orientation is None:
                                    orientation="forward"
                                elif not orientation=="forward":
                                    orientation = "conflict"
                                    continue
                            elif newArmKmers[rckmer[0]]==rckmer[1]:
                                if orientation is None:
                                    orientation="backward"
                                elif not orientation=="backward":
                                    orientation = "conflict"
                                    continue
                            else:
                                orientation = "conflict"
                                continue
                        else:
                            orientation = "conflict"
                            continue
                if orientation=="forward":
                    if arm.armType()=="incoming":
                        if entry[1][0]=="u" and not entry[1][-1]=="u":
                            if "a" in entry[1] and "u" in entry[1][entry[1].index("a"):]:
                                pass #unknown after first arm k-mer
                            elif "c" in entry[1] and "a" in entry[1][entry[1].index("c"):]:
                                pass #arm k-mer after first candidate k-mer
                            elif "c" in entry[1] and "u" in entry[1][entry[1].index("c"):]:
                                pass #unknown after first candidate k-mer
                            elif "a" in entry[1]: #should at least contain arm-kmer
                                orientatedReads.append([[item.copy() for item in entry[0]],entry[1].copy(),entry[2]])
                        elif "a" in entry[1] and "c" in entry[1] and not "u" in entry[1]: 
                            if "a" in entry[1][entry[1].index("c"):]:
                                pass #arm k-mer after candidate k-mer
                            else:
                                orientatedReads.append([[item.copy() for item in entry[0]],entry[1].copy(),entry[2]])
                    elif arm.armType()=="outgoing":
                        if entry[1][-1]=="u" and not entry[1][0]=="u":
                            if "a" in entry[1] and "u" in entry[1][:entry[1].index("a")]:
                                pass #unknown before first arm k-mer
                            elif "c" in entry[1] and "a" in entry[1][:entry[1].index("c")]:
                                pass #arm k-mer before first candidate k-mer
                            elif "c" in entry[1] and "u" in entry[1][:entry[1].index("c")]:
                                pass #unknown before first candidate k-mer
                            elif "a" in entry[1]: #should at least contain arm-kmer
                                orientatedReads.append([[item.copy() for item in entry[0]],entry[1].copy(),entry[2]])
                        elif "a" in entry[1] and "c" in entry[1] and not "u" in entry[1]:
                            if "a" in entry[1][:entry[1].index("c")]:
                                pass #arm k-mer before candidate k-mer
                            else:
                                orientatedReads.append([[item.copy() for item in entry[0]],entry[1].copy(),entry[2]])
                elif orientation=="backward":
                    reverseReads = [item.copy() for item in entry[0][::-1]]
                    reverseTypes = [item for item in entry[1][::-1]]
                    maxPosition = entry[0][-1]["position"]
                    for i in range(len(reverseReads)):
                        reverseReads[i]["position"] = maxPosition - reverseReads[i]["position"]
                        reverseReads[i]["orientation"] = ("forward" if reverseReads[i]["orientation"]=="backward" else (
                                                          "backward" if reverseReads[i]["orientation"]=="forward" else "unknown"))
                    if arm.armType()=="incoming":
                        if reverseTypes[0]=="u" and not reverseTypes[-1]=="u":
                            if "a" in reverseTypes and "u" in reverseTypes[reverseTypes.index("a"):]:
                                pass #unknown after first arm k-mer
                            elif "c" in reverseTypes and "a" in reverseTypes[reverseTypes.index("c"):]:
                                pass #arm k-mer after first candidate k-mer
                            elif "c" in reverseTypes and "u" in reverseTypes[reverseTypes.index("c"):]:
                                pass #unknown after first candidate k-mer
                            elif "a" in reverseTypes: #should at least contain arm-kmer
                                orientatedReads.append([reverseReads,reverseTypes,entry[2]])
                        elif "a" in reverseTypes and "c" in reverseTypes and not "u" in reverseTypes:
                            if "a" in reverseTypes[:reverseTypes.index("c")]:
                                pass #arm k-mer before candidate k-mer
                            else:
                                orientatedReads.append([reverseReads,reverseTypes,entry[2]])
                    elif arm.armType()=="outgoing":
                        if reverseTypes[-1]=="u" and not reverseTypes[0]=="u":
                            if "a" in reverseTypes and "u" in reverseTypes[:reverseTypes.index("a")]:
                                pass #unknown before first arm k-mer
                            elif "c" in reverseTypes and "a" in reverseTypes[:reverseTypes.index("c")]:
                                pass #arm k-mer before first candidate k-mer
                            elif "c" in reverseTypes and "u" in reverseTypes[:reverseTypes.index("c")]:
                                pass #unknown before first candidate k-mer
                            elif "a" in reverseTypes: #should at least contain arm-kmer
                                orientatedReads.append([reverseReads,reverseTypes,entry[2]])
                        elif "a" in reverseTypes and "c" in reverseTypes and not "u" in reverseTypes:
                            if "a" in reverseTypes[:reverseTypes.index("c")]:
                                pass #arm k-mer before candidate k-mer
                            else:
                                orientatedReads.append([reverseReads,reverseTypes,entry[2]])
            
            #build response
            response = []
            for entry in orientatedReads:
                kmerList = []
                typeList = []
                positionList = []
                readNumberList = []
                for kmerInfo,kmerType,frequency in zip(entry[0],entry[1],entry[2]):
                    kmerList.append((kmerInfo["ckmer"],kmerInfo["orientation"],kmerInfo["split"],kmerInfo["number"]))
                    typeList.append(kmerType)
                    positionList.append(kmerInfo["position"])
                    readNumberList.append(frequency)
                response.append([kmerList,typeList,positionList,readNumberList])
                    
            return response
    
        def expandGlueResponses(arm,glueResponse):
            for item in glueResponse:
                if not "u" in item[1]:
                    if arm.armType()=="outgoing":
                        expanded = True
                        while expanded:
                            expanded = False
                            kmer = item[0][-1]
                            direct = self._api.getSplitDirect(self._uid,kmer[0])
                            if direct:
                                if kmer[1]=="forward":
                                    options = direct.get("direct",{}).get("right",{})
                                else:
                                    options = direct.get("direct",{}).get("left",{})
                                if len(options)==1:
                                    expanded = True
                                    if options[0]["connection"]["direction"]=="left":
                                        item[0].append((options[0]["ckmer"],"forward",options[0]["split"],options[0]["number"]))
                                    else:
                                        item[0].append((options[0]["ckmer"],"backward",options[0]["split"],options[0]["number"]))
                                    item[1].append("u")
                                    item[2].append(item[2][-1]+options[0]["connection"]["distance"])
                                    item[3].append(options[0]["connection"]["number"])
                    elif arm.armType()=="incoming":
                        expanded = True
                        while expanded:
                            expanded = False
                            kmer = item[0][0]
                            direct = self._api.getSplitDirect(self._uid,kmer[0])
                            if direct:
                                if kmer[1]=="forward":
                                    options = direct.get("direct",{}).get("left",{})
                                else:
                                    options = direct.get("direct",{}).get("right",{})
                                if len(options)==1:
                                    expanded = True
                                    if options[0]["connection"]["direction"]=="right":
                                        item[0] = [(options[0]["ckmer"],"forward",options[0]["split"],options[0]["number"])] + item[0]
                                    else:
                                        item[0] = [(options[0]["ckmer"],"backward",options[0]["split"],options[0]["number"])] + item[0]
                                    item[1] = ["u"] + item[1]
                                    item[2] = [0] + [x+options[0]["connection"]["distance"] for x in item[2]]
                                    item[3] = [options[0]["connection"]["number"]] + item[3]
    
        def glue(reads,arm,initialList=None):
            #initiate
            starters = []
            others = [item for item in reads]
            #helper function
            def _getNewStarters(entries):
                hasNeighbour = set()
                for entry in entries:
                    for i in range(1,len(entry[0])):
                        if arm.armType()=="incoming":
                            hasNeighbour.add(entry[0][i-1])
                        elif arm.armType()=="outgoing":
                            hasNeighbour.add(entry[0][i])
                if arm.armType()=="incoming":
                    newStarters = [entry for entry in entries if not entry[0][-1] in hasNeighbour]
                    newOthers = [entry for entry in entries if entry[0][-1] in hasNeighbour]
                elif arm.armType()=="outgoing":
                    newStarters = [entry for entry in entries if not entry[0][0] in hasNeighbour]
                    newOthers = [entry for entry in entries if entry[0][0] in hasNeighbour]
                return newStarters,newOthers    
            #glue reads
            glueList = []
            if not initialList is None:
                for entry in initialList:
                    glueList.append(entry)
            #process all reads
            while True:
                if len(others)==0:
                    break
                else:
                    starters,others = _getNewStarters(others)
                    if len(starters)==0:
                        break
                    else:
                        #first check if read contained, update counters
                        uncontainedStarters = []
                        for entry in starters:
                            contained = False
                            if arm.armType()=="incoming":
                                for glueEntry in glueList:
                                    if entry[0][-1] in glueEntry[0]:
                                        id = glueEntry[0].index(entry[0][-1])
                                        glueCandidate = glueEntry[0][:id+1]
                                        if entry[0] == glueCandidate[-len(entry[0]):]:
                                            contained = True
                                            for i in range(len(entry[0])):
                                                glueEntry[3][i+1+id-len(entry[0])] += entry[3][i]
                            elif arm.armType()=="outgoing":
                                for glueEntry in glueList:
                                    if entry[0][0] in glueEntry[0]:
                                        id = glueEntry[0].index(entry[0][0])
                                        glueCandidate = glueEntry[0][id:]
                                        if entry[0] == glueCandidate[0:len(entry[0])]:
                                            contained = True
                                            for i in range(len(entry[0])):
                                                glueEntry[3][i+id] += entry[3][i]
                            if not contained:
                                uncontainedStarters.append(entry)
                        #then, try to glue
                        newGlueEntries = []
                        for entry in starters:
                            glued = False
                            if arm.armType()=="incoming":
                                for glueEntry in glueList:
                                    if not "u" in glueEntry[1] and not initialList is None:
                                        pass
                                    elif entry[0][-1] in glueEntry[0]:
                                        id = glueEntry[0].index(entry[0][-1])
                                        glueCandidate = glueEntry[0][:id+1]
                                        if glueCandidate == entry[0][-len(glueCandidate):]:
                                            glued = True
                                            newGlueEntry = [[],[],[],[]]
                                            for i in range(len(entry[0])):
                                                newGlueEntry[0].append(tuple(entry[0][i]))
                                                newGlueEntry[1].append(entry[1][i])
                                                newGlueEntry[2].append(entry[2][i])
                                                newGlueEntry[3].append(entry[3][i])
                                            for i in range(id+1):
                                                newGlueEntry[3][len(entry[0])+i-id-1] += glueEntry[3][i]
                                            positionStart = entry[2][-1] - glueEntry[2][id]
                                            for i in range(id+1,len(glueEntry[0])):
                                                newGlueEntry[0].append(tuple(glueEntry[0][i]))
                                                newGlueEntry[1].append(glueEntry[1][i])
                                                newGlueEntry[2].append(glueEntry[2][i]+positionStart)
                                                newGlueEntry[3].append(glueEntry[3][i])
                                            if len([True for item in newGlueEntries if item[0]==newGlueEntry[0]])==0:
                                                newGlueEntries.append(newGlueEntry)
                            elif arm.armType()=="outgoing":
                                for glueEntry in glueList:
                                    if not "u" in glueEntry[1] and not initialList is None:
                                        pass
                                    elif entry[0][0] in glueEntry[0]:
                                        id = glueEntry[0].index(entry[0][0])
                                        glueCandidate = glueEntry[0][id:]
                                        if glueCandidate == entry[0][0:len(glueCandidate)]:
                                            glued = True
                                            newGlueEntry = [
                                                [tuple(item) for item in glueEntry[0][0:id]],
                                                glueEntry[1][0:id].copy(),
                                                glueEntry[2][0:id].copy(),
                                                glueEntry[3][0:id].copy()]
                                            positionStart = glueEntry[2][id] - entry[2][0]
                                            for i in range(len(entry[0])):
                                                newGlueEntry[0].append(tuple(entry[0][i]))
                                                newGlueEntry[1].append(entry[1][i])
                                                newGlueEntry[2].append(entry[2][i]+positionStart)
                                                newGlueEntry[3].append(entry[3][i])
                                            for i in range(id,len(glueEntry[0])):
                                                newGlueEntry[3][i]+=glueEntry[3][i]
                                            if len([True for item in newGlueEntries if item[0]==newGlueEntry[0]])==0:
                                                newGlueEntries.append(newGlueEntry)
                            if initialList is None and not glued:
                                newGlueEntries.append(entry)
                        #now combine with glueEntries
                        glueList.extend(newGlueEntries)
                        #and filter this by removin unnecessary subsets
                        filteredGlueList = []
                        glueListReport = []
                        for i in range(len(glueList)):
                            reportEntry = {"isRealSubset": set(), "unconfirmed": False, 
                                           "ignore": False, "duplicate": set(), "subsets": set()}
                            for k in range(len(glueList[i][0])):
                                if glueList[i][3][k]<minimumReadDepth:
                                    if glueList[i][1][k]=="u":
                                        reportEntry["unconfirmed"] = True
                                    elif glueList[i][1][k]=="a":
                                        if arm.armType()=="incoming":
                                            if (not "u" in glueList[i][1]) or (k>glueList[i][1].index("a")+numberOfKmers-1):
                                                reportEntry["ignore"] = True
                                        elif arm.armType()=="outgoing":
                                            if (not "u" in glueList[i][1]) or (k<glueList[i][1].index("u")-numberOfKmers):
                                                reportEntry["ignore"] = True
                            for j in range(len(glueList)):
                                if not i==j:
                                    if glueList[i][0]==glueList[j][0]:
                                        if i>j:
                                            reportEntry["duplicate"].add(j)
                                    else:
                                        if glueList[i][0][0] in glueList[j][0]:
                                            #todo: loop over all matches
                                            p = glueList[j][0].index(glueList[i][0][0])
                                            if glueList[i][0]==glueList[j][0][p:p+len(glueList[i][0])]:
                                                reportEntry["isRealSubset"].add(j)
                            glueListReport.append(reportEntry)
                        #get confirmed subsets
                        for i in range(len(glueList)):
                            if not (glueListReport[i]["unconfirmed"] or glueListReport[i]["ignore"]):
                                for j in glueListReport[i]["isRealSubset"]:
                                    glueListReport[j]["subsets"].add(i)
        
                        includedSubsets = set()
                        for i in range(len(glueList)):
                            if len(glueListReport[i]["duplicate"])>0:
                                pass
                            elif glueListReport[i]["ignore"]:
                                pass
                            elif len(glueListReport[i]["isRealSubset"])==0 and len(glueListReport[i]["duplicate"])==0:
                                filteredGlueList.append(glueList[i])
                                if glueListReport[i]["unconfirmed"] and len(glueListReport[i]["subsets"])>0:
                                    #compute largest confirmed subset(s)
                                    subsets = glueListReport[i]["subsets"]
                                    for j in glueListReport[i]["subsets"]:
                                        subsets = subsets.difference(glueListReport[j]["subsets"])
                                    #todo: select
                                    for j in subsets:
                                        if len(glueListReport[j]["duplicate"])==0:
                                            includedSubsets.add(j)
                        for i in range(len(includedSubsets)):
                            filteredGlueList.append(glueList[i])
                        glueList = filteredGlueList
        
            expandGlueResponses(arm,glueList)                
            
            return glueList
    
        def expandArm(arm,glueResponse,allNewKmersArm,n,processedReads):
            newKmersArm = {}
            #get new k-mers for arm
            while len(newKmersArm)<n:
                newKmers = 0
                for item in glueResponse:
                    if "u" in item[1]:
                        id = (item[1].index("u") if arm.armType()=="outgoing" 
                              else (len(item[1]) - item[1][::-1].index("u") - 1))
                        newKmersArm[item[0][id][0]] = item[0][id][1]
                        item[1][id] = "a"
                        newKmers+=1
                if newKmers==0:
                    break
            #check glue responses for missing k-mers in selection
            while True:
                newKmers = 0
                for item in glueResponse:
                    for i in range(len(item[0])):
                        if item[0][i][0] in newKmersArm:
                            item[1][i] = "a"
                    if arm.armType()=="outgoing":
                        while "u" in item[1] and "a" in item[1] and item[1].index("u")<(len(item[1])-item[1][::-1].index("a")-1):
                            id = item[1].index("u")
                            if not item[0][id][0] in newKmersArm:
                                newKmersArm[item[0][id][0]] = item[0][id][1]
                                newKmers+=1
                            item[1][id] = "a"
                    elif arm.armType()=="incoming":
                        while "u" in item[1] and "a" in item[1] and (len(item[1])-item[1][::-1].index("u")-1)>item[1].index("a"):
                            id = (len(item[1])-item[1][::-1].index("u")-1)
                            if not item[0][id][0] in newKmersArm:
                                newKmersArm[item[0][id][0]] = item[0][id][1]
                                newKmers+=1
                            item[1][id] = "a"
                if newKmers==0:
                    break
            #get all reads
            newReadsArm = self._api.getSplitReads(self._uid,newKmersArm.keys())
            allNewKmersArm.update(newKmersArm)
            responsenew = processArmReads(newReadsArm,arm,allNewKmersArm)
            responsefiltered = []
            for processedRead in responsenew:
                entry = tuple([item[0] for item in processedRead[0]])
                if not entry in processedReads:
                    processedReads.add(entry)
                    responsefiltered.append(processedRead)
            glueResponse = glue(responsefiltered,arm,glueResponse)
            return glueResponse,allNewKmersArm,len(newKmersArm),len(responsefiltered)
    
        def getPathLengths(glueResponse):
            if len(glueResponse)>0:
                minLength = min([item[2][-1] for item in glueResponse])
                maxLength = max([item[2][-1] for item in glueResponse])
            else:
                minLength = 0
                maxLength = 0
            length = "{}".format(maxLength) if minLength==maxLength else "{}-{}".format(minLength,maxLength)
            return minLength,maxLength,length
    
        def getSolutions(glueResponse1,glueResponse2):
            solutions = []
            for path1 in glueResponse1:
                kmerset1 = set([item[0] for item in path1[0]])
                for path2 in glueResponse2:
                    kmerset2 = set([item[0] for item in path2[0]])
                    if len(kmerset1.intersection(kmerset2))>0:
                        for i in range(len(path1[0])):
                            if path1[0][i] in path2[0]:
                                j = path2[0].index(path1[0][i])
                                p1 = path1[0][i:]
                                p2 = path2[0][j:j+len(p1)]
                                if p1==p2 and (minimumGlueSize==0 or len(p1)>=minimumGlueSize):
                                    solution = path1[0] + path2[0][j+len(p1):]
                                    numbers = path1[3][0:i]
                                    for k in range(i,i+len(p1)):
                                        numbers.append(max(path1[3][k],path2[3][j-i+k]))
                                    numbers = numbers + path2[3][j+len(p1):]
                                    if min(numbers)>=minimumReadDepth:
                                        solutions.append(solution)
                                break
            return solutions
    
        def processSolution(solution):
            #get direct connections
            kmers = set()
            for i in range(len(solution)):
                kmers.add(solution[i][0])
            directConnections = {item["ckmer"]: item for item in self._api.getSplitDirect(self._uid,kmers)}
            #process
            previousCkmer = None
            for i in range(len(solution)):
                item = solution[i]
                ckmer = self._createOrientatedCkmer(item[0], item[1], item[3], item[2])
                if previousCkmer:
                    if previousCkmer._key in ckmer._incoming and ckmer._key in previousCkmer._outgoing:
                        pass
                    else:
                        direct = directConnections.get(item[0],{})
                        if item[1]=="forward":
                            entries = direct.get("direct",{}).get("left",[])
                        else:
                            entries = direct.get("direct",{}).get("right",[])
                        connection = None
                        for entry in entries:
                            if entry["ckmer"]==previousCkmer._ckmer:
                                if (previousCkmer._orientation=="forward") and (entry.get("connection",{}).get("direction","")=="right"):
                                    connection = entry.get("connection",{})
                                    break
                                elif (previousCkmer._orientation=="backward") and (entry.get("connection",{}).get("direction","")=="left"):
                                    connection = entry.get("connection",{})
                                    break
                        if connection:
                            if not previousCkmer._key in ckmer._incoming:
                                ckmer._setIncoming(previousCkmer._key,connection["distance"], connection["number"], connection["problem"]>0)
                            if not ckmer._key in previousCkmer._outgoing:
                                previousCkmer._setOutgoing(ckmer._key,connection["distance"], connection["number"], connection["problem"]>0)
                        else:
                            #should not happen
                            pass
                previousCkmer = ckmer
    
        #main function
        self._detectArms()
        armPairs = self._detectConnectedArmsCandidates()
        for pairCounter in range(len(armPairs)):
            solutions = []
            self._logger.debug("trying to connect {} of {} pair(s) of arms".format(pairCounter+1,len(armPairs)))
            arm1 = self.getArm(armPairs[pairCounter][0])
            arm2 = self.getArm(armPairs[pairCounter][1])
            arm1Kmers = set([item[0] for item in arm1._orientatedCkmers])
            arm2Kmers = set([item[0] for item in arm2._orientatedCkmers])
            assert arm1.armType()=="outgoing"
            assert arm2.armType()=="incoming"
            #get orientated k-mers
            candidatesArm1 = set([entry for entry in [arm1.connection()] 
                                      if entry in self.getCandidates()])
            candidatesArm2 = set([entry for entry in [arm2.connection()] 
                                      if entry in self.getCandidates()])
            for i in range(numberOfCandidates):
                newCandidatesArm1 = set()
                for entry in candidatesArm1:
                    newCandidatesArm1.update([newEntry for newEntry in self._orientatedCkmers[entry]._incoming
                                              if newEntry in self.getCandidates()])
                candidatesArm1.update(newCandidatesArm1)
                newCandidatesArm2 = set()
                for entry in candidatesArm2:
                    newCandidatesArm2.update([newEntry for newEntry in self._orientatedCkmers[entry]._outgoing
                                              if newEntry in self.getCandidates()])
                candidatesArm2.update(newCandidatesArm2)
            #get plain k-mers
            candidateKmersArm1 = set([entry[0] for entry in candidatesArm1])
            candidateKmersArm2 = set([entry[0] for entry in candidatesArm2])
    
            readsArm1 = self._api.getSplitReads(self._uid,arm1Kmers.union(candidateKmersArm1))
            response1 = processArmReads(readsArm1,arm1)
            processedReads1 = set([tuple([item[0] for item in processedRead[0]]) for processedRead in response1])
            glueResponse1 = glue(response1,arm1)
            glueResponse1 = [entry for entry in glueResponse1 if "c" in entry[1]]
            for entry in glueResponse1:
                for i in range(len(entry[0])):
                    if entry[0][i][0] in arm1Kmers or entry[0][i][0] in candidateKmersArm1:
                        entry[3][i]=max(entry[3][i],minimumReadDepth)
            minLength1,maxLength1,length1 = getPathLengths(glueResponse1)
            self._logger.debug("starting with {} outgoing path(s) of length {}".format(len(glueResponse1),len(response1),length1))
            
            readsArm2 = self._api.getSplitReads(self._uid,arm2Kmers.union(candidateKmersArm2))
            response2 = processArmReads(readsArm2,arm2)
            processedReads2 = set([tuple([item[0] for item in processedRead[0]]) for processedRead in response2])
            glueResponse2 = glue(response2,arm2)
            glueResponse2 = [entry for entry in glueResponse2 if "c" in entry[1]]
            for entry in glueResponse2:
                for i in range(len(entry[0])):
                    if entry[0][i][0] in arm2Kmers or entry[0][i][0] in candidateKmersArm2:
                        entry[3][i]=max(entry[3][i],minimumReadDepth)
            minLength2,maxLength2,length2 = getPathLengths(glueResponse2)
            self._logger.debug("starting with {} incoming path(s) of length {}".format(len(glueResponse2),len(response2),length2))
            
            allNewKmersArm1 = {}
            allNewKmersArm2 = {}
            counter = 0
            solutions = []
            while True:
                counter+=1
                glueResponse1,allNewKmersArm1,newKmersArm1,newReads1 = expandArm(
                    arm1,glueResponse1,allNewKmersArm1,numberOfKmers,processedReads1)
                minLength1,maxLength1,length1 = getPathLengths(glueResponse1)
                self._logger.debug("step {}: expand {} k-mers with {} reads to {} outgoing path(s) of length {}".format(
                    counter,newKmersArm1,newReads1,len(glueResponse1),length1))
                glueResponse2,allNewKmersArm2,newKmersArm2,newReads2 = expandArm(
                    arm2,glueResponse2,allNewKmersArm2,numberOfKmers,processedReads2)
                minLength2,maxLength2,length2 = getPathLengths(glueResponse2)
                self._logger.debug("step {}: expand {} k-mers with {} reads to {} incoming path(s) of length {}".format(
                    counter,newKmersArm2,newReads2,len(glueResponse2),length2))
                #no solution possible (for now)
                if newKmersArm1==0 and newKmersArm2==0:
                    break
                #check
                kmerset1 = set()
                kmerset2 = set()
                for entry in glueResponse1:
                    kmerset1.update([item[0] for item in entry[0]])
                for entry in glueResponse2:
                    kmerset2.update([item[0] for item in entry[0]])
                solutions = getSolutions(glueResponse1,glueResponse2)
                #check for solutions
                if len(solutions)>0:
                    break
                #check limits
                if maximumSteps>0 and counter>=maximumSteps:
                    break
                elif maximumLength>0 and maxLength1+maxLength2>maximumLength:
                    break
                elif maximumPaths>0 and len(glueResponse1)+len(glueResponse2)>maximumPaths:
                    break
            if len(solutions)>0:
                self._logger.debug("{} solution(s) for {} of {} pair(s) of arms".format(len(solutions),pairCounter+1,len(armPairs)))
                for solution in solutions:
                    self._logger.debug("add solution with {} nodes".format(len(solution)))
                    processSolution(solution)
            else:
                self._logger.debug("no solution for {} of {} pair(s) of arms".format(pairCounter+1,len(armPairs)))
    
        #reset arms
        self._resetArms()
        self._detectArms()


    def getArm(self, armId):
        if armId>=0 and armId<=len(self._arms):
            return self._arms[armId]
        else:
            return None

    def _detectConnectedArmsCandidates(self):
        return []
        
    def _resetDistances(self):
        if not self._connected is None:
            self._logger.debug("reset computed connections for the De Bruijn Graph")  
            self._connected = None
        
#     def _computeDistances(self):
#         #initialise
#         orientatedCkmerList = [c for c in self._orientatedCkmers]
#         self._distances = {
#             "directDistance": pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0),
#             "directNumber": pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0),
#             "shortest": pd.DataFrame(index = orientatedCkmerList, columns=orientatedCkmerList).fillna(0)
#         }
#         self._logger.debug("compute distances for the De Bruijn Graph")  
#         #direct
#         for orientatedCkmer in orientatedCkmerList:
#             for connectedCkmer in self._orientatedCkmers[orientatedCkmer]._incoming:
#                 if connectedCkmer in orientatedCkmerList:
#                     d = self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["distance"]
#                     n = self._orientatedCkmers[orientatedCkmer]._incoming[connectedCkmer]["number"]
#                     self._distances["directDistance"].at[orientatedCkmer,connectedCkmer] = -d
#                     self._distances["directDistance"].at[connectedCkmer,orientatedCkmer] = d
#                     self._distances["directNumber"].at[orientatedCkmer,connectedCkmer] = n
#                     self._distances["directNumber"].at[connectedCkmer,orientatedCkmer] = n
#         self._logger.debug("found 2 x {:d} graph direct connected pairs".
#               format(int(np.count_nonzero(self._distances["directDistance"]>0))))   
#         #compute distances between all nkmers
#         nkg = nk.graph.Graph(n=len(orientatedCkmerList), weighted=True, directed=True)
#         #add relations to graph
#         for i in range(len(orientatedCkmerList)):
#             for j in range(len(orientatedCkmerList)):
#                 if self._distances["directDistance"].values[i,j]>0:
#                     nkg.addEdge(i,j,self._distances["directDistance"].values[i,j])
#         nkd=nk.distance.APSP(nkg)
#         nkd.run()
#         distance = nkd.getDistances()
#         #remove maxint and add negative distances
#         for i in range(len(distance)):
#             for j in range(len(distance[i])):
#                 if distance[i][j]>10**100:
#                     distance[i][j]=0
#                 elif distance[i][j]>0:
#                     distance[j][i]=-1*distance[i][j]
#         self._distances["shortest"] = pd.DataFrame(data=distance, index=orientatedCkmerList, 
#                                                    columns=orientatedCkmerList).astype(int)
#         self._logger.debug("found 2 x {:d} = {:d} graph disconnected pairs".
#               format(int((np.count_nonzero(self._distances["shortest"]==0)-len(self._distances["shortest"].index))/2),
#                      np.count_nonzero(self._distances["shortest"]==0)-len(self._distances["shortest"].index)))
        
        
    def _createOrientatedCkmer(self, ckmer: str, orientation: str, number: int, split: str):
        ckmerKey = (ckmer,orientation)
        if ckmerKey in self._orientatedCkmers:
            return self._orientatedCkmers[ckmerKey]
        else:
            self._resetDistances()
            return OrientatedCkmer(self,ckmer,orientation,number,split)
        
    class CanonicalKmer:
        """internal object representing splitting canonical k-mers in the De Bruijn graph"""
        
        SIDE = Literal["left","right"]
        
        def __init__(self, graph, ckmer: str, number: int, split: str):
            """
            - ckmer: canonical representation of the splitting k-mer
            - number: frequency
            - split: location of the split (left/right/both)
            """
            #checks
            assert len(ckmer)==graph._k
            assert split in ["left","right","both"]
            assert number>=0
            #initialise
            self._graph = graph
            self._ckmer : str = str(ckmer)
            self._split : str = str(split)
            self._number : int = int(number)
            #other variables
            self._left : dict = {}
            self._right : dict = {}
            self._orientated : set = set([])
            #check for duplicates
            if self._ckmer in self._graph._ckmers.keys():
                raise Exception("creating canonical k-mer that already exists")
            #register
            self._graph._ckmers[self._ckmer] = self
            
        def __repr__(self):
            info = []
            info.append("{}x".format(self._number))
            info.append("split {}".format(self._split))
            return "CanonicalCkmer({}[{}])".format(self._ckmer,", ".join(info))
        
        def _setLeft(self, ckmer, side: SIDE, distance: int, number: int, problem: bool = None):
            assert ckmer in self._graph._ckmers.keys()
            assert distance>0
            assert number>=0
            #set or update connection
            if (ckmer in self._left.keys()) and (side in self._left[ckmer].keys()):
                if not self._left[ckmer][side]["problem"]:
                    assert self._left[ckmer][side]["distance"] == distance
                    assert self._left[ckmer][side]["number"] == number
                elif not problem is None:
                    self._left[ckmer][side]["problem"] = problem
            else:    
                if not ckmer in self._left.keys():
                    self._left[ckmer] = {}
                self._left[ckmer][side] = {
                    "distance": distance, "number": number, 
                    "problem": problem
                }  
            #check existence orientated versions
            for keyFrom in self._orientated:
                if keyFrom[1]=="forward":
                    if side=="left":
                        keyTo = (ckmer, "backward")
                    else:
                        keyTo = (ckmer, "forward")
                else:
                    if side=="left":
                        keyTo = (ckmer, "forward")
                    else:
                        keyTo = (ckmer, "backward")
                if not keyTo in self._graph._ckmers[ckmer]._orientated:
                    self._graph._ckmers[ckmer]._orientate(keyTo[1])
            for keyTo in self._graph._ckmers[ckmer]._orientated:
                if keyTo[1]=="forward":
                    if side=="left":
                        keyFrom = (ckmer, "backward")
                    else:
                        keyFrom = (ckmer, "forward")
                else:
                    if side=="left":
                        keyFrom = (ckmer, "forward")
                    else:
                        keyFrom = (ckmer, "backward")
                if not keyFrom in self._orientated:
                    self._graph._ckmers[ckmer]._orientate(keyTo[1])
            
        def _setRight(self, ckmer, side: SIDE, distance: int, number: int, problem: bool = None):
            assert ckmer in self._graph._ckmers.keys()
            assert distance>0
            assert number>=0
            #set or update connection
            if (ckmer in self._right.keys()) and (side in self._right[ckmer].keys()):
                if not self._right[ckmer][side]["problem"]:
                    assert self._right[ckmer][side]["distance"] == distance
                    assert self._right[ckmer][side]["number"] == number
                elif not problem is None:
                    self._right[ckmer][side]["problem"] = problem
            else:    
                if not ckmer in self._right.keys():
                    self._right[ckmer] = {}
                self._right[ckmer][side] = {
                    "distance": distance, "number": number, 
                    "problem": problem
                }
            #check existence orientated versions
            for keyFrom in self._orientated:
                if keyFrom[1]=="forward":
                    if side=="left":
                        keyTo = (ckmer, "forward")
                    else:
                        keyTo = (ckmer, "backward")
                else:
                    if side=="left":
                        keyTo = (ckmer, "backward")
                    else:
                        keyTo = (ckmer, "forward")
                if not keyTo in self._graph._ckmers[ckmer]._orientated:
                    self._graph._ckmers[ckmer]._orientate(keyTo[1])
            for keyTo in self._graph._ckmers[ckmer]._orientated:
                if keyTo[1]=="forward":
                    if side=="left":
                        keyFrom = (ckmer, "forward")
                    else:
                        keyFrom = (ckmer, "backward")
                else:
                    if side=="left":
                        keyFrom = (ckmer, "backward")
                    else:
                        keyFrom = (ckmer, "forward")
                if not keyFrom in self._orientated:
                    self._graph._ckmers[ckmer]._orientate(keyTo[1])
                
        def _orientate(self, orientation):
            ckmerKey = (self._ckmer, orientation)
            if not ckmerKey in self._graph._orientatedCkmers:
                orientatedCkmer = self._graph._createOrientatedCkmer(self._ckmer, orientation, self._number, self._split)
            if not ckmerKey in self._orientated:
                self._orientated.add(ckmerKey)
            
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
            self._graph._resetDistances()
            if not self._ckmer in self._graph._ckmers.keys():
                self._canonicalKmer = self._graph.CanonicalKmer(graph, ckmer, number, split)
            else:
                self._canonicalKmer = self._graph._ckmers[self._ckmer]
            self._canonicalKmer._orientated.add(self._key)
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
            assert self._incomingArm is None
            self._incomingArm = arm
        
        def _setOutgoingArm(self, arm):
            assert self._outgoingArm is None
            self._outgoingArm = arm
            
        def _setIncomingDeadEnd(self, orientatedCkmer, distance: int, path: str):
            assert len(self._incoming)==0
            self._incomingDeadEnd = (orientatedCkmer, distance, path)
            #consistency check
            if not path is None:
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
            if not path is None:
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
            assert path is None or len(path)==distance+self._graph._k
            assert self._incomingDeadEnd is None
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
                if not path is None:
                    assert path == computedPath
                else:
                    path = computedPath
            #set or update connection
            if orientatedCkmer in self._incoming.keys():
                if not self._incoming[orientatedCkmer]["problem"]:
                    assert self._incoming[orientatedCkmer]["distance"] == distance
                    assert self._incoming[orientatedCkmer]["number"] == number
            else:    
                self._incoming[orientatedCkmer] = {
                    "distance": distance, "number": number, 
                    "problem": problem, "path": None
                }
                self._graph._resetDistances()
            if not problem is None:
                self._incoming[orientatedCkmer]["problem"] = problem
            if not path is None:
                self._incoming[orientatedCkmer]["path"] = path    
            #update start and end based on connections
            if orientatedCkmer in self._graph._end:
                self._setPostEnd(1)
            elif not self._graph._orientatedCkmers[orientatedCkmer]._postEnd is None:
                self._setPostEnd(self._graph._orientatedCkmers[orientatedCkmer]._postEnd+1)
            if self._key in self._graph._start:
                self._graph._orientatedCkmers[orientatedCkmer]._setPreStart(1)
            elif not self._preStart is None:
                self._graph._orientatedCkmers[orientatedCkmer]._setPreStart(self._preStart+1)
            #update canonical k-mers
            if self._orientation=="forward":
                if orientatedCkmer[1]=="forward":
                    self._canonicalKmer._setLeft(orientatedCkmer[0], "right", distance, number, problem)
                else:
                    self._canonicalKmer._setLeft(orientatedCkmer[0], "left", distance, number, problem)
            else:
                if orientatedCkmer[1]=="forward":
                    self._canonicalKmer._setRight(orientatedCkmer[0], "right", distance, number, problem)
                else:
                    self._canonicalKmer._setRight(orientatedCkmer[0], "left", distance, number, problem)
            
        def _setOutgoing(self, orientatedCkmer, distance: int, number: int, problem: bool = None, path: str = None):
            assert orientatedCkmer in self._graph._orientatedCkmers.keys()
            assert distance>0
            assert number>=0
            assert path is None or len(path)==distance+self._graph._k
            assert self._outgoingDeadEnd is None
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
                if not path is None:
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
                self._graph._resetDistances()
            if not problem is None:
                self._outgoing[orientatedCkmer]["problem"] = problem
            if not path is None:
                self._outgoing[orientatedCkmer]["path"] = path    
            #update start and end based on connections
            if orientatedCkmer in self._graph._start:
                self._setPreStart(1)
            elif not self._graph._orientatedCkmers[orientatedCkmer]._preStart is None:
                self._setPreStart(self._graph._orientatedCkmers[orientatedCkmer]._preStart+1)
            if self._key in self._graph._end:
                self._graph._orientatedCkmers[orientatedCkmer]._setPostEnd(1)
            elif not self._postEnd is None:
                self._graph._orientatedCkmers[orientatedCkmer]._setPostEnd(self._postEnd+1)
            #update canonical k-mers
            if self._orientation=="forward":
                if orientatedCkmer[1]=="forward":
                    self._canonicalKmer._setRight(orientatedCkmer[0], "left", distance, number, problem)
                else:
                    self._canonicalKmer._setRight(orientatedCkmer[0], "right", distance, number, problem)
            else:
                if orientatedCkmer[1]=="forward":
                    self._canonicalKmer._setLeft(orientatedCkmer[0], "left", distance, number, problem)
                else:
                    self._canonicalKmer._setLeft(orientatedCkmer[0], "right", distance, number, problem)
            
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
                
        def armId(self):
            if self.arm():
                return [arm.id() for arm in self._arm]
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
            if (self._preStart is None) or (self._preStart>steps):
                self._preStart = steps
                for direction in self._orientatedBases.keys():
                    orientatedBaseKey = self._orientatedBases[direction]
                    self._graph._orientatedBases[orientatedBaseKey]._setPreStart(steps)
                for incomingOrientatedCkmer in self._incoming.keys():
                    self._graph._orientatedCkmers[incomingOrientatedCkmer]._setPreStart(steps+1)
                    
        def _setPostEnd(self, steps: int):
            if (self._postEnd is None) or (self._postEnd>steps):
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
            elif (self._preStart is None) or (self._preStart>steps):
                self._preStart = steps
                for orientatedCkmer in self._orientatedCkmers:
                    self._graph._orientatedCkmers[orientatedCkmer]._setPreStart(steps)

        def _setPostEnd(self,steps):
            if self.candidate():
                pass
            elif (self._postEnd is None) or (self._postEnd>steps):
                self._postEnd = steps   
                for orientatedCkmer in self._orientatedCkmers:
                    self._graph._orientatedCkmers[orientatedCkmer]._setPostEnd(steps)
            
        def _setOrder(self):
            self._order = None
            for orientatedCkmer in self._orientatedCkmers:
                if not self._graph._orientatedCkmers[orientatedCkmer]._order is None:
                    if self._order is None:
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
                if not self._graph._orientatedCkmers[orientatedCkmer]._order is None:
                    if self._order is None:
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
        
        TYPE = Literal["incoming","outgoing"]
        
        def __init__(self, graph, orientatedCkmer, type: TYPE):
            self._type = type
            self._graph = graph
            self._orientatedCkmers = set()
            self._connection = None
            self._maxFreq = 0
            self._size = 0
            self._id = None
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
            self._id = len(self._graph._arms)
            self._graph._arms.append(self)
            
        def __repr__(self):
            if self._type=="incoming":
                text = "Incoming Arm"
            elif self._type=="outgoing":
                text = "Outgoing Arm"
            else:
                text = "Arm"
            text = "{}, size {} with {} nodes, maximum frequency {}".format(
                text, self.size(), self.n(), self.maxFreq())
            return text
            
        def _add(self, orientatedCkmer):
            assert orientatedCkmer in self._graph._orientatedCkmers
            if self._type == "incoming":
                self._graph._orientatedCkmers[orientatedCkmer]._setArm(self)
            elif self._type == "outgoing":
                self._graph._orientatedCkmers[orientatedCkmer]._setArm(self)
            self._orientatedCkmers.add(orientatedCkmer)
            self._maxFreq = max(self._maxFreq,self._graph._orientatedCkmers[orientatedCkmer]._number)
            self._size = None
            
        def armType(self):
            return self._type
        
        def maxFreq(self):
            return self._maxFreq
        
        def id(self):
            return self._id
        
        def n(self):
            return len(self._orientatedCkmers)
        
        def size(self):
            if self._size is None:
                orientatedCkmerList = set(self._orientatedCkmers)
                orientatedCkmerList.add(self._connection)
                orientatedCkmerList = list(orientatedCkmerList)
                m = np.zeros((len(orientatedCkmerList),len(orientatedCkmerList)), dtype = "uint8") 
                for i in range(len(orientatedCkmerList)):
                    orientatedCkmer = self._graph._orientatedCkmers[orientatedCkmerList[i]]
                    for outgoing,properties in orientatedCkmer._outgoing.items():
                        if outgoing in orientatedCkmerList:
                            j = orientatedCkmerList.index(outgoing)
                            m[i,j] = properties["distance"]
                G = nx.from_numpy_array(m)
                self._size = nx.eccentricity(G,orientatedCkmerList.index(self._connection))
            return self._size
        
        def connection(self):
            return self._connection
            
            
            