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
        self._orientatedCkmers : dict = {}
        self._orientatedBases : dict = {}
        self._selected : set = set()
        self._connected = None
        self._directDistances = None
        self._directNumbers = None
        self._arms = []
        self._connections = []
            
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
                    ockmer1 = (ckmer,"forward")
                    ockmer2 = (ckmer,"backward")
                    if ockmer1 in self._orientatedCkmers:
                        selection.add(ockmer1)
                    if ockmer2 in self._orientatedCkmers:
                        selection.add(ockmer2)
            elif hasattr(config["selection"], "__iter__"):
                for kmer in config["selection"]:
                    if isinstance(kmer,str) and (len(kmer)==self._k):
                        ckmer = General.canonical(kmer)
                        ockmer1 = (ckmer,"forward")
                        ockmer2 = (ckmer,"backward")
                        if ockmer1 in self._orientatedCkmers:
                            selection.add(ockmer1)
                        if ockmer2 in self._orientatedCkmers:
                            selection.add(ockmer2)
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
            "edgeArrowSize": 0.5,
            "edgeFontSizeConnection": 10,
            "edgeStyleConnection":  "dashed",
            "edgeColorConnection": "orange",
            "edgePenWidthConnection": 5
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
                    #indirect connections
                    for connection in self._orientatedCkmers[orientatedCkmer1]._outgoingConnections:
                        orientatedCkmer2 = connection._orientatedCkmerTo
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

                            edge_label = self._visualize_connection_label(connection, {"size": config["edgeFontSizeConnection"]}, config)
                            edge_direction = "forward"
                            edge_constraint="true"
                            edge_color = config["edgeColorConnection"]
                            edge_style = config["edgeStyleConnection"]
                            edge_arrowsize = config["edgeArrowSize"]
                            edge_penwidth = config["edgePenWidthConnection"]
                            #create the edge
                            g.edge(node_key1, node_key2, label=edge_label,
                                   constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth)) 
                    #direct connections            
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
                    #indirect connections
                    for connection in self._orientatedCkmers[orientatedCkmer1]._outgoingConnections:
                        orientatedCkmer2 = connection._orientatedCkmerTo
                        if orientatedCkmer2 in connectedCandidates[i]["connected"]:
                            node_key2 = self._visualize_node_key(config["prefix"],orientatedCkmer2)
                            #only pass once
                            edgeKey = (node_key1, node_key2)
                            if edgeKey in edges:
                                continue
                            else:
                                edges.add(edgeKey)

                            edge_label = self._visualize_connection_label(connection, {"size": config["edgeFontSizeConnection"]}, config)
                            edge_direction = "forward"
                            edge_constraint="true"
                            edge_color = config["edgeColorConnection"]
                            edge_style = config["edgeStyleConnection"]
                            edge_arrowsize = config["edgeArrowSize"]
                            edge_penwidth = config["edgePenWidthConnection"]
                            #create the edge
                            g.edge(node_key1, node_key2, label=edge_label,
                                   constraint=edge_constraint, style=edge_style, arrowsize=str(edge_arrowsize),
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth)) 
                    #direct connections 
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
        else:
            processedArms = set()
            #show arms
            if config["showArms"]:
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
            "edgeNumberFontSize": "8",
            "edgeFontSizeConnection": 10,
            "edgeStyleConnection":  "dashed",
            "edgeColorConnection": "orange",
            "edgePenWidthConnection": 5
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
                #indirect connections
                for connection in self._orientatedCkmers[orientatedCkmer]._incomingConnections:
                    connectedCkmer = connection._orientatedCkmerFrom
                    if connectedCkmer in orientatedCkmerNodes.keys():     
                        edge_label = self._visualize_connection_label(connection, {"size": config["edgeFontSizeConnection"]}, config)
                        edge_direction = "forward"
                        edge_constraint="true"
                        edge_color = config["edgeColorConnection"]
                        edge_style = config["edgeStyleConnection"]
                        edge_penwidth = config["edgePenWidthConnection"]
                        #choose preferred nodes for the edges
                        if backward_key in orientatedCkmerNodes[orientatedCkmer]:
                            node = orientatedCkmerNodes[orientatedCkmer][backward_key]
                        else:
                            node = orientatedCkmerNodes[orientatedCkmer][forward_key]
                        if forward_key in orientatedCkmerNodes[connectedCkmer]:
                            connectedNode = orientatedCkmerNodes[connectedCkmer][forward_key]
                        else:
                            connectedNode = orientatedCkmerNodes[connectedCkmer][backward_key]
                        #create the edge only once
                        edgeKey = (connectedNode, node)
                        if not edgeKey in edges:
                            edges.add(edgeKey)
                            g.edge(connectedNode, node, label=edge_label,
                                   constraint=edge_constraint, style=edge_style,
                                   dir=edge_direction, color=edge_color, penwidth=str(edge_penwidth)) 
                #direct connections  
                for connectedCkmer in self._orientatedCkmers[orientatedCkmer]._incoming.keys():
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
        else:
            processedArms = set()
            #show arms
            if config["showArms"]:
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

    def _visualize_connection_label(self, connection, fontSize:dict = {}, config:dict = {}):
        edge_label = "<"
        if connection._sequence:
            s = len(connection._sequence) - connection._graph._k
            edge_label += "<font point-size=\"{}\" color=\"grey\">{}</font>".format(fontSize.get("size","0"), s)
        elif connection._sequenceFrom or connection._sequenceTo:
            s1 = (len(connection._sequenceFrom) - connection._graph._k) if connection._sequenceFrom else 0
            s2 = (len(connection._sequenceTo) - connection._graph._k) if connection._sequenceTo else 0
            if s1>0:
                edge_label += "<font point-size=\"{}\" color=\"grey\">{}</font><br/>".format(fontSize.get("size","0"), s1)
            edge_label += "<font point-size=\"{}\" color=\"grey\">&#183;&#183;&#183;</font>".format(fontSize.get("size","0"))
            if s2>0:
                edge_label += "<br/><font point-size=\"{}\" color=\"grey\">{}</font>".format(fontSize.get("size","0"), s2)
        else:
            edge_label += "<font point-size=\"{}\" color=\"grey\">??</font>".format(fontSize.get("size","0"))
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
            connected[i,i] = True
        self._connected = pd.DataFrame(connected, index=orientatedCkmerList, columns=orientatedCkmerList)
        
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

    def getArm(self, armId):
        if armId>=0 and armId<=len(self._arms):
            return self._arms[armId]
        else:
            return None

    def _detectConnectedArmsCandidates(self):
        return []

    def _resetConnections(self):
        for orientatedCkmer in self._orientatedCkmers:
            self._orientatedCkmers[orientatedCkmer]._incomingConnections = []
            self._orientatedCkmers[orientatedCkmer]._outgoingConnections = []
        self._connections = []

    def createConnection(self,fromOrientatedCkmer,toOrientatedCkmer):
        assert fromOrientatedCkmer in self._orientatedCkmers
        assert toOrientatedCkmer in self._orientatedCkmers
        for connection in self._connections:
            if (connection in self._orientatedCkmers[fromOrientatedCkmer]._outgoingConnections
                and
                connection in self._orientatedCkmers[toOrientatedCkmer]._incomingConnections):
                return connection
        return self.Connection(self,fromOrientatedCkmer,toOrientatedCkmer)
    
    def getConnections(self):
        return self._connections
        
    def _resetDistances(self):
        if not self._connected is None:
            self._logger.debug("reset computed connections for the De Bruijn Graph")  
            self._connected = None
        
    def _createOrientatedCkmer(self, ckmer: str, orientation: str, number: int, split: str):
        ckmerKey = (ckmer,orientation)
        if ckmerKey in self._orientatedCkmers:
            return self._orientatedCkmers[ckmerKey]
        else:
            self._resetDistances()
            return OrientatedCkmer(self,ckmer,orientation,number,split)

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
            self._incomingConnections = []
            self._outgoingConnections = []
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

        def _setIncomingConnection(self, connection):
            if not connection in self._incomingConnections:
                self._incomingConnections.append(connection)
        
        def _setOutgoingConnection(self, connection):
            if not connection in self._outgoingConnections:
                self._outgoingConnections.append(connection)
            
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

    class Connection:
        """
        internal object representing non-direct connection
        """
        
        def __init__(self, graph, orientatedCkmerFrom, orientatedCkmerTo):
            self._graph = graph
            self._sequenceFrom = None
            self._orientatedCkmersFrom = None
            self._sequenceTo = None
            self._orientatedCkmersTo = None
            self._sequence = None
            self._orientatedCkmers = None
            self._id = None
            #register from
            assert orientatedCkmerFrom in self._graph._orientatedCkmers
            assert self._graph._orientatedCkmers[orientatedCkmerFrom].candidate()
            self._graph._orientatedCkmers[orientatedCkmerFrom]._setOutgoingConnection(self)
            self._orientatedCkmerFrom = orientatedCkmerFrom
            #register to
            assert orientatedCkmerTo in self._graph._orientatedCkmers
            assert self._graph._orientatedCkmers[orientatedCkmerTo].candidate()
            self._graph._orientatedCkmers[orientatedCkmerTo]._setIncomingConnection(self)
            self._orientatedCkmerTo = orientatedCkmerTo
            #register connection
            self._id = len(self._graph._connections)
            self._graph._connections.append(self)
            
        def __repr__(self):
            text = "Connection"
            if not self._sequence is None:
                text = "{} - {}bp".format(text,len(self._sequence))
            elif not (self._sequenceFrom is None and self._sequenceTo is None):
                if self._sequenceFrom is None:
                    text = "{} - right {}bp known".format(text,len(self._sequenceTo))
                elif self._sequenceTo is None:
                    text = "{} - left {}bp known".format(text,len(self._sequenceFrom))
                else:
                    text = "{} - left {}bp and right {}bp known".format(text,len(self._sequenceFrom),len(self._sequenceTo))
            return text

        def setSequence(self,sequence:str,orientatedCkmers:list):
            self._sequence = sequence
            self._orientatedCkmers = orientatedCkmers
            self._sequenceFrom = None
            self._orientatedCkmersFrom = None
            self._sequenceTo = None
            self._orientatedCkmersTo = None

        def setSequenceFrom(self,sequence:str, orientatedCkmers:list):
            self._sequence = None
            self._orientatedCkmers = None
            self._sequenceFrom = sequence
            self._orientatedCkmersFrom = orientatedCkmers

        def setSequenceTo(self,sequence:str, orientatedCkmers:list):
            self._sequence = None
            self._orientatedCkmers = None
            self._sequenceTo = sequence
            self._orientatedCkmersTo = orientatedCkmers

        def getOrientatedCkmers(self):
            orientatedCkmers = []
            if self._orientatedCkmers:
                for ockmer in self._orientatedCkmers:
                    if ockmer in self._graph._orientatedCkmers:
                        if not self._graph._orientatedCkmers[ockmer].candidate():
                            orientatedCkmers.append(ockmer)
            else:
                if self._orientatedCkmersFrom:
                    for ockmer in self._orientatedCkmersFrom:
                        if ockmer in self._graph._orientatedCkmers:
                            if not self._graph._orientatedCkmers[ockmer].candidate():
                                orientatedCkmers.append(ockmer)
                if self._orientatedCkmersTo:
                    for ockmer in self._orientatedCkmersTo:
                        if ockmer in self._graph._orientatedCkmers:
                            if not self._graph._orientatedCkmers[ockmer].candidate():
                                orientatedCkmers.append(ockmer)
            return orientatedCkmers

        def getOrientatedCkmersFrom(self):
            orientatedCkmers = []
            if self._orientatedCkmersFrom:
                for ockmer in self._orientatedCkmersFrom:
                    if ockmer in self._graph._orientatedCkmers:
                        if not self._graph._orientatedCkmers[ockmer].candidate():
                            orientatedCkmers.append(ockmer)
            elif self._orientatedCkmers:
                for ockmer in self._orientatedCkmers:
                    if ockmer in self._graph._orientatedCkmers:
                        if not self._graph._orientatedCkmers[ockmer].candidate():
                            orientatedCkmers.append(ockmer)
                    else:
                        break
            return orientatedCkmers

        def getOrientatedCkmersTo(self):
            orientatedCkmers = []
            if self._orientatedCkmersTo:
                for ockmer in self._orientatedCkmersTo:
                    if ockmer in self._graph._orientatedCkmers:
                        if not self._graph._orientatedCkmers[ockmer].candidate():
                            orientatedCkmers.append(ockmer)
            elif self._orientatedCkmers:
                for ockmer in reversed(list(self._orientatedCkmers)):
                    if ockmer in self._graph._orientatedCkmers:
                        if not self._graph._orientatedCkmers[ockmer].candidate():
                            orientatedCkmers.append(ockmer)
                    else:
                        break
                orientatedCkmers = list(reversed(orientatedCkmers))
            return orientatedCkmers
            
        def id(self):
            return self._id
            
            
            