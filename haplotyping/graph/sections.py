import logging,requests,re, html
import haplotyping.graph.baseGraph as baseGraph
from haplotyping.general import General
from graphviz import Digraph

class Sections():
    """constructing linear ordered sections from candidates in the De Bruijn graph"""
    
    def __init__(self, graph: baseGraph.Graph):
        """
        construct sections
        """
        assert isinstance(graph, baseGraph.Graph)
        #logger
        self._logger = logging.getLogger(__name__)
        
        self._graph = graph
        self._parts = []
        
        #create sections
        for baseConnectedCandidates in self._graph.getConnectedCandidates(True):
            self._createSections(baseConnectedCandidates)
            
            
    def getOrientatedCkmers(self):
        orientatedCkmers = set()
        for i in range(self.partNumber()):
            orientatedCkmers.update(self.part(i+1).getOrientatedCkmers())
        return orientatedCkmers
        
    def getCandidates(self):
        candidates = []
        for orientatedCkmer in self.getOrientatedCkmers():
            if self._graph._orientatedCkmers[orientatedCkmer].candidate():
                candidates.append(orientatedCkmer)
        return set(candidates)
        
    def _createSections(self, connectedCandidates):
        #compute section start entries (not every connected start is necessarily a section start)
        shortestDistances = self._graph.getShortestDistances()
        startEntries = connectedCandidates["start"]
        problematicStartEntries = set()
        for orientatedCkmer in startEntries:
            baseLinkedOrientatedCkmers = set()
            for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                baseLinkedOrientatedCkmers.update(self._graph._orientatedBases[orientatedBase]._orientatedCkmers)
            for otherOrientatedCkmer in baseLinkedOrientatedCkmers:
                otherIncomingConnected = [k for k in connectedCandidates["connected"] 
                                          if shortestDistances[k][otherOrientatedCkmer]<0]
                if len(otherIncomingConnected)>0:
                    problematicStartEntries.add(orientatedCkmer)
                    break
        if len(problematicStartEntries)>0:
            if len(problematicStartEntries)==len(startEntries):
                self._logger.warning("all {} start entries sections problematic".format(len(startEntries)))
            else:
                startEntries = startEntries.difference(problematicStartEntries)
                self._logger.debug("ignore {} start entries, continue with {} for sections".format(
                    len(problematicStartEntries),len(startEntries)))
        #now try to create sections       
        part = self.Part(self._graph)   
        partOrientatedCkmers = set()
        sectionItem = self.SectionItem(startEntries, connectedCandidates["start"], connectedCandidates["end"],
                                       connectedCandidates["connected"],self._graph)
        part._append(sectionItem)
        partOrientatedCkmers.update(sectionItem._orientatedCkmers)
        while(len(sectionItem._end["ckmers"].difference(connectedCandidates["end"]))>0):
            newStartCkmers = (sectionItem._end["ckmers"]
                              .intersection(connectedCandidates["connected"])
                              .difference(sectionItem._start["ckmers"]))
            sectionItem = self.SectionItem(newStartCkmers, connectedCandidates["start"], connectedCandidates["end"],
                                       connectedCandidates["connected"],self._graph)
            part._append(sectionItem)
            partOrientatedCkmers.update(sectionItem._orientatedCkmers)
        #compute paths
        part._computePaths()
        self._logger.info("detected {} paths in part with {} sections".format(
            part._pathNumber,len(part._sectionItems)))
        #check consistency
        missingConnected = connectedCandidates["connected"].difference(partOrientatedCkmers)
        if len(missingConnected)>0:
            self._logger.warning("sections not covering complete set of connected candidates, {} missing".format(
                len(missingConnected)))
        #add finalized sections
        self._parts.append(part)
        #try to order logically
        self._parts = sorted(self._parts, key=lambda x: 0 if x._order==None else x._order)
                       
    def partNumber(self):
        return len(self._parts)
    
    def part(self, i=1):
        if i>0 and i<=len(self._parts):
            return self._parts[i-1]
        else:
            return None
        
    def visualize(self, *args, **kwargs):  
        
        g = Digraph(name="cluster_graph")
        graph_label = "Sectioned Graph"
        if self._graph._name:
            graph_label = "{} {}".format(graph_label,html.escape(self._graph._name))
        graph_label = "{} ({} part{})".format(graph_label,len(self._parts),"s" if len(self._parts)>1 else "")
        graph_label = "<<font point-size=\"12\" color=\"grey\">{}</font>>".format(graph_label)
        g.attr(label=graph_label, labelloc="t", nodesep="0", ranksep="0")        
                
        for i in range(len(self._parts)):
            kwargs["containerGraph"] = g
            kwargs["prefix"] = "graph_part_{}".format(i)
            kwargs["label"] = "Part {}".format(i+1)
            self._parts[i].visualize(*args, **kwargs)
        return g
                       
    
    class Part():
        
        def __init__(self, graph):
            """
            construct part
            """
            #logger
            self._logger = logging.getLogger(__name__)
            
            self._graph = graph
            self._sectionItems = []
            self._order : int = None
            self._pathNumber : int = None
                
        def __repr__(self):
            text = "Part graph"
            if self._graph._name:
                text = "{} {}".format(text,self._graph._name)
            text = "{}; containing {} section{} with {} path{}".format(text,len(self._sectionItems),
                   "s" if len(self._sectionItems)>1 else "", self._pathNumber, "s" if self._pathNumber>1 else "")
            return text
            
        def _append(self, item):
            self._sectionItems.append(item)
            if not item._order == None:
                self._order = (item._order if self._order==None else min(self._order,item._order))
                
        def _computePaths(self):
            self._pathNumber = 0
            numbers = {sectionItemPath.getOrientatedCkmers()[0]: 1 for sectionItemPath in self._sectionItems[0]._paths}
            for sectionItem in self._sectionItems:
                newNumbers = {sectionItemPath.getOrientatedCkmers()[-1]: 0 for sectionItemPath in sectionItem._paths}
                for sectionItemPath in sectionItem._paths:
                    orientatedCkmers = sectionItemPath.getOrientatedCkmers()
                    newNumbers[orientatedCkmers[-1]] = (newNumbers[orientatedCkmers[-1]] 
                                                               + numbers.get(orientatedCkmers[0],1))
                numbers = newNumbers
                for orientatedCkmer in numbers:
                    if not orientatedCkmer in sectionItem._end["ckmers"]:
                        self._pathNumber += numbers[orientatedCkmer]
            #add all remaining   
            for orientatedCkmer in numbers:
                if orientatedCkmer in sectionItem._end["ckmers"]:
                    self._pathNumber += numbers[orientatedCkmer]            
            
        def number(self):
            return len(self._sectionItems)
        
        def section(self, i):
            if i>0 and i<=self.number():
                return self._sectionItems[i-1]
        
        def getOrientatedCkmers(self):
            orientatedCkmers = set()
            for sectionItem in self._sectionItems:
                orientatedCkmers.update(sectionItem._orientatedCkmers)
            return orientatedCkmers
        
        def getCandidates(self):
            candidates = []
            for orientatedCkmer in self.getOrientatedCkmers():
                if self._graph._orientatedCkmers[orientatedCkmer].candidate():
                    candidates.append(orientatedCkmer)
            return set(candidates)
        
        def visualize(self, *args, **kwargs):  
        
            config = {
                "baseStyleIncomingArm": "filled",
                "baseFillColorIncomingArm": "lightblue",
                "basePenWidthIncomingArm": 1,
                "baseColorIncomingArm": "black",
                "baseStyleOutgoingArm": "filled",
                "baseFillColorOutgoingArm": "lightcoral",
                "basePenWidthOutgoingArm": 1,
                "baseColorOutgoingArm": "black",
                "baseStylePart": "filled",
                "baseFillColorPart": "beige",
                "basePenWidthPart": 1,
                "baseColorPart": "grey",
                "containerGraph": False,
                "showAllBases": False,
                "prefix": "graph",
                "label": None
            }
            for key, value in kwargs.items():
                if key in config.keys():
                    config[key] = value
                    
            #create graph
            if config["containerGraph"]:
                g = config["containerGraph"]
            else:
                g = Digraph()
            
            partGraph=g.subgraph(name="cluster_{}".format(config["prefix"]))      
            
            with partGraph as pg:
                if config["label"]:
                    graph_label = config["label"]
                    graph_label = "{} with {} section{}".format(graph_label, self.number(), 
                                                                "s" if self.number()>1 else "")
                else:
                    graph_label = "{} section{}".format(self.number(),
                                                                "s" if self.number()>1 else "")
                    if self._graph._name:
                        graph_label = "{} for Graph {}".format(graph_label,self._graph._name)
                graph_label = ("<" + 
                                  "<font point-size=\"12\" color=\"grey\">{}</font><br/>".format(graph_label) +
                                  "<font point-size=\"10\">{:,} paths</font><br/>".format(self._pathNumber) +
                                  ">")
                graph_style = config["baseStylePart"]
                graph_fillcolor = config["baseFillColorPart"]
                graph_penwidth = config["basePenWidthPart"]
                graph_color = config["baseColorPart"]
                pg.attr(label=graph_label, style=graph_style, fillcolor=graph_fillcolor, 
                        color=graph_color, penwidth=str(graph_penwidth), labelloc="t", nodesep="0", ranksep="0")

                #arms
                orientatedCkmersSections = self.getOrientatedCkmers()
                for j in range(len(self._graph._arms)):
                    if self._graph._arms[j].connection() in orientatedCkmersSections:
                        #only if visible
                        if not config["showAllBases"]:
                            orientatedCkmer = self._graph._arms[j].connection()
                            hideArm = False
                            for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                if not self._graph._orientatedBases[orientatedBase].candidate():
                                    hideArm = True
                            if hideArm:
                                continue
                        #create arm
                        armGraph=pg.subgraph(name="cluster_graph_arm_{}".format(self._graph._arms[j].key()))
                        with armGraph as ag:
                            arm_name = "Arm {}".format(j+1)
                            if self._graph._arms[j].armType()=="incoming":
                                arm_name = "Incoming {}".format(arm_name)
                                arm_style = config["baseStyleIncomingArm"]
                                arm_fillcolor = config["baseFillColorIncomingArm"]
                                arm_penwidth = config["basePenWidthIncomingArm"]
                                arm_color = config["baseColorIncomingArm"]
                            elif self._graph._arms[j].armType()=="outgoing":
                                arm_name = "Outgoing {}".format(arm_name)
                                arm_style = config["baseStyleOutgoingArm"]
                                arm_fillcolor = config["baseFillColorOutgoingArm"]
                                arm_penwidth = config["basePenWidthOutgoingArm"]
                                arm_color = config["baseColorOutgoingArm"]
                            else:
                                continue
                            arm_label = ("<" + 
                                         "<font point-size=\"12\" color=\"grey\">{}</font><br/>".format(arm_name) + 
                                         "<font point-size=\"10\">length: {} with {} nodes</font><br/>".format(
                                             self._graph._arms[j].size(),self._graph._arms[j].n()) + 
                                         "<font point-size=\"8\">maximum frequency: {}x</font><br/>".format(
                                             self._graph._arms[j].maxFreq()) + 
                                         ">")
                            ag.attr(label=arm_label, style=arm_style, fillcolor=arm_fillcolor, 
                                color=arm_color, penwidth=str(arm_penwidth), labelloc="t", nodesep="0", ranksep="0")

                            ag.node("arm_{}_{}".format(self._graph._arms[j].armType(),
                                                       self._graph._arms[j].key()), shape="point")
                #sections
                previousPrefix = ""
                for j in range(self.number()):
                    sectionItem = self.section(j+1)
                    section_prefix = "graph_section_{}".format(j)
                    section_label = "Section {} of {}".format(j+1,self.number())
                    sharedStart = (self._sectionItems[j-1]._orientatedCkmers.intersection(sectionItem._orientatedCkmers) 
                                   if j>0 else set())
                    sharedEnd = (self._sectionItems[j+1]._orientatedCkmers.intersection(sectionItem._orientatedCkmers) 
                                   if j<(self.number()-1) else set())
                    kwargs["containerGraph"] = pg
                    kwargs["prefix"] = section_prefix
                    kwargs["label"] = section_label
                    kwargs["hideDeadEndBefore"] = sharedStart
                    kwargs["hideDeadEndAfter"] = sharedEnd
                    #create section graph
                    self._sectionItems[j].visualize(*args, **kwargs)
                    #detect arm connections
                    for orientatedCkmer in sectionItem._orientatedCkmers:
                        if self._graph._orientatedCkmers[orientatedCkmer].arm():
                            for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                #only if visible
                                if not config["showAllBases"]:
                                    if not self._graph._orientatedBases[orientatedBase].candidate():
                                        continue 
                                #connect arm
                                ckmerKey = "{}_{}_{}_{}_{}".format(section_prefix,
                                            orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                if self._graph._orientatedCkmers[orientatedCkmer].incomingArmType():
                                    for armKey in self._graph._orientatedCkmers[orientatedCkmer].armKey():                     
                                        fullArmKey = "arm_incoming_{}".format(armKey)
                                        if not orientatedCkmer in sharedStart:
                                            pg.edge(fullArmKey,ckmerKey,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                                elif self._graph._orientatedCkmers[orientatedCkmer].outgoingArmType():
                                    for armKey in self._graph._orientatedCkmers[orientatedCkmer].armKey():                     
                                        fullArmKey = "arm_outgoing_{}".format(armKey)
                                        if not orientatedCkmer in sharedEnd:
                                            pg.edge(ckmerKey,fullArmKey,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                                #only one connection if multiple bases
                                break
                        else:
                            incomingArm = self._graph._orientatedCkmers[orientatedCkmer].incomingArm()
                            if incomingArm:
                                for orientatedBase in \
                                         self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                    ckmerKey = "{}_{}_{}_{}_{}".format(section_prefix,
                                            orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                    armKey = "arm_incoming_{}".format(incomingArm.key())
                                    if (not orientatedCkmer in sharedStart and 
                                        not len(incomingArm._orientatedCkmers.intersection(
                                            orientatedCkmersSections))>0):
                                        pg.edge(armKey,ckmerKey,style="dashed", 
                                           color="grey", rankdir="lr", constraint="true")
                                    break
                            outgoingArm = self._graph._orientatedCkmers[orientatedCkmer].outgoingArm()
                            if outgoingArm:
                                for orientatedBase in \
                                         self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                    ckmerKey = "{}_{}_{}_{}_{}".format(section_prefix,
                                            orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                    armKey = "arm_outgoing_{}".format(outgoingArm.key())
                                    if (not orientatedCkmer in sharedEnd and 
                                        not len(outgoingArm._orientatedCkmers.intersection(
                                            orientatedCkmersSections))>0):
                                        pg.edge(ckmerKey,armKey,style="dashed", 
                                           color="grey", rankdir="lr", constraint="true")
                                    break
                    #link start/end sections
                    for orientatedCkmer in sharedStart:
                        for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                            key1 = "{}_{}_{}_{}_{}".format(previousPrefix,
                                    orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                            key2 = "{}_{}_{}_{}_{}".format(section_prefix,
                                    orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                            pg.edge(key1,key2,style="dashed", color="grey", rankdir="tb", constraint="true")
                    previousPrefix = section_prefix                        
                        
            return g
            
        
    class SectionItem():
    
        def __init__(self, initialCkmers, startCkmers, endCkmers, connectedCkmers, graph):
            """
            construct section item
            """
            #logger
            self._logger = logging.getLogger(__name__)
            
            self._start = {"ckmers": set(), "bases": set()}
            self._end = {"ckmers": set(), "bases": set()}
            self._orientatedCkmers = set()
            self._orientatedBases = set()
            
            self._connectedCkmers = connectedCkmers
            self._connectedBases = set()
            self._startCkmers = startCkmers
            self._endCkmers = endCkmers
            self._graph = graph
            self._order : int = None
            self._paths : None
            self._pathNumber : int = None
            self._optionalPathNumber : int = None
            self._requiredPathNumber : int = None            
            
            self._checkBackwardList = set()
            self._checkForwardList = {}           
            
            #compute connected bases
            for connectedCkmer in connectedCkmers:
                self._connectedBases.update(self._graph._orientatedCkmers[connectedCkmer]._orientatedBases.values())

            self._logger.debug("create section with {} starting k-mers".format(len(startCkmers)))
            for initialCkmer in initialCkmers:
                assert initialCkmer in graph._orientatedCkmers
                self._addOrientatedCkmer(initialCkmer,start=True)
                assert initialCkmer in self._orientatedCkmers
                assert initialCkmer in self._start["ckmers"]

            while len(self._checkBackwardList)>0 or len(self._checkForwardList)>0:
            
                self._checkBackwardList = self._checkBackwardList.difference(self._orientatedCkmers)
                while len(self._checkBackwardList)>0:
                    checkBackwardList = self._checkBackwardList.intersection(self._connectedCkmers)
                    self._checkBackwardList = set()
                    self._logger.debug("check backward {}".format(len(checkBackwardList)))
                    for checkBackwardItem in checkBackwardList:
                        self._addOrientatedCkmer(checkBackwardItem)
                        
                #check if section finished
                allEndingCkmers = set()
                for endingCkmers in self._checkForwardList.values():
                    allEndingCkmers.update(endingCkmers)
                allEndingCkmers = allEndingCkmers.difference(self._start["ckmers"])
                allEndingBases = set()
                for endingCkmer in allEndingCkmers:
                    allEndingBases.update(self._graph._orientatedCkmers[endingCkmer]._orientatedBases.values())
                if len(allEndingBases)==1:
                    for endingCkmer in allEndingCkmers:
                        self._end["ckmers"].add(endingCkmer)
                    self._end["bases"].update(allEndingBases)
                    self._logger.debug("close section with {} k-mers".format(len(self._end["ckmers"])))
                    self._checkForwardList = {}                    
            
                if len(self._checkForwardList)>0:                                
                    checkForwardList = list(self._checkForwardList.keys())
                    #first check for section end
                    forwardBases = set()
                    for checkForwardItem in checkForwardList:
                        assert checkForwardItem in self._graph._orientatedCkmers
                        forwardBases.update(self._graph._orientatedCkmers[checkForwardItem]._orientatedBases.values()) 
                    self._logger.debug("check forward {}".format(len(checkForwardList)))
                    for checkForwardItem in checkForwardList:                        
                        self._addOrientatedCkmer(checkForwardItem)
                        
            #compute paths
            self._computePaths()
            
            #set order
            for orientatedCkmer in self._orientatedCkmers:
                if not self._graph._orientatedCkmers[orientatedCkmer]._order == None:
                    self._order = (self._graph._orientatedCkmers[orientatedCkmer]._order if self._order==None else
                                   min(self._graph._orientatedCkmers[orientatedCkmer]._order,self._order))
                    
                
        def __repr__(self):
            text = "Section object graph"
            if self._graph._name:
                text = "{} {}".format(text,self._graph._name)
            text = "{}; containing {} paths".format(text,self._pathNumber)
            return text
                            
        def _addOrientatedCkmer(self, orientatedCkmer, start=False):
            for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                assert orientatedBase in self._graph._orientatedBases
                if start:
                    self._start["bases"].add(orientatedBase)
                self._orientatedBases.add(orientatedBase)
                assert orientatedCkmer in self._graph._orientatedBases[orientatedBase]._orientatedCkmers
                for orientatedBaseCkmer in self._graph._orientatedBases[orientatedBase]._orientatedCkmers:
                    assert orientatedBaseCkmer in self._graph._orientatedCkmers
                    self._orientatedCkmers.add(orientatedBaseCkmer)
                    if orientatedBaseCkmer in self._checkForwardList.keys():
                        del self._checkForwardList[orientatedBaseCkmer]
                    if start:
                        self._start["ckmers"].add(orientatedBaseCkmer)
                    if orientatedCkmer in self._endCkmers:
                        self._end["ckmers"].add(orientatedCkmer)
                        self._end["bases"].add(orientatedBase)
                    else:
                        for outgoingCkmer in self._graph._orientatedCkmers[orientatedBaseCkmer]._outgoing:
                            outgoingBases = self._graph._orientatedCkmers[outgoingCkmer]._orientatedBases.values()
                            if len(self._connectedBases.intersection(outgoingBases))>0:
                                if not outgoingCkmer in self._orientatedCkmers:
                                    if outgoingCkmer in self._checkForwardList:
                                        self._checkForwardList[outgoingCkmer].add(orientatedBaseCkmer)
                                    else:
                                        self._checkForwardList[outgoingCkmer] = set([orientatedBaseCkmer])
                    if not orientatedBase in self._start["bases"]:
                        for incomingCkmer in self._graph._orientatedCkmers[orientatedBaseCkmer]._incoming:
                            if incomingCkmer in self._connectedCkmers:
                                if not incomingCkmer in self._orientatedCkmers:
                                    self._checkBackwardList.add(incomingCkmer)
                                
        def _computePaths(self):
            self._paths = []
            #starters are connected, within section and sectionStarter or real starter
            sectionStarters = self._connectedCkmers.intersection(self._orientatedCkmers).intersection(self._start["ckmers"])
            realStarters = self._connectedCkmers.intersection(self._orientatedCkmers).intersection(self._startCkmers)
            partialPaths = [{"distance": 0,
                             "sequence": (orientatedCkmer[0] if orientatedCkmer[1]=="forward" 
                                      else General.reverse_complement(orientatedCkmer[0])),
                             "list": [orientatedCkmer]} for orientatedCkmer in 
                                        sectionStarters.union(realStarters)]
            while len(partialPaths)>0:
                newPartialPaths = []
                for partialPath in partialPaths:
                    if partialPath["list"][-1] in self._end["ckmers"]:
                            self._paths.append(self.SectionItemPath(partialPath["list"],
                                                               partialPath["sequence"],partialPath["distance"],
                                                               self._graph))
                    else:
                        orientatedCkmerInfo = self._graph._orientatedCkmers[partialPath["list"][-1]]
                        newPaths = 0
                        for orientatedCkmer in orientatedCkmerInfo._outgoing:
                            if orientatedCkmer in self._orientatedCkmers and orientatedCkmer in self._connectedCkmers:
                                outgoingInfo = orientatedCkmerInfo._outgoing[orientatedCkmer]
                                assert partialPath["sequence"][-self._graph._k:]==outgoingInfo["path"][:self._graph._k]
                                newPartialPath = {"distance": partialPath["distance"] + outgoingInfo["distance"],
                                                  "sequence": partialPath["sequence"] + outgoingInfo["path"][self._graph._k:],
                                                  "list": partialPath["list"] + [orientatedCkmer]}
                                newPartialPaths.append(newPartialPath)
                                newPaths+=1
                        if newPaths==0:
                            self._paths.append(self.SectionItemPath(partialPath["list"],
                                                               partialPath["sequence"],partialPath["distance"],
                                                               self._graph))
                partialPaths = newPartialPaths 
            #update statistics
            self._pathNumber = len(self._paths)
            self._requiredPathNumber = 0
            self._optionalPathNumber = 0
            frequencies = {orientatedCkmer:[] for orientatedCkmer in 
                           self._connectedCkmers.intersection(self._orientatedCkmers)}
            for i in range(len(self._paths)):
                for orientatedCkmer in self._paths[i].getOrientatedCkmers():
                    frequencies[orientatedCkmer].append(i)
            coveredOrientatedCkmers = set()
            for orientatedCkmer in frequencies:
                if len(frequencies[orientatedCkmer])==1:
                    if not self._paths[frequencies[orientatedCkmer][0]].required():
                        self._paths[frequencies[orientatedCkmer][0]]._setRequired()
                        self._requiredPathNumber+=1
                        coveredOrientatedCkmers.update(self._paths[frequencies[orientatedCkmer][0]].getOrientatedCkmers())
            for i in range(len(self._paths)):
                if not self._paths[i].required():
                    if len(set(self._paths[i].getOrientatedCkmers()).difference(coveredOrientatedCkmers))==0:
                        self._paths[i]._setOptional()
                        self._optionalPathNumber+=1
        
        def pathNumber(self):
            return self._pathNumber
        
        def path(self, i: int):
            if i>0 and i<=self._pathNumber:
                return self._paths[i-1]
            else:
                return None
            
        def getOrientatedCkmers(self):
            return self._orientatedCkmers                        
        
        def getCandidates(self):
            candidates = []
            for orientatedCkmer in self.getOrientatedCkmers():
                if self._graph._orientatedCkmers[orientatedCkmer].candidate():
                    candidates.append(orientatedCkmer)
            return set(candidates)
            
        def visualize(self, *args, **kwargs):  
            
            config = {
                "baseStyleSection": "filled",
                "baseFillColorSection": "grey98",
                "basePenWidthSection": 1,
                "baseColorSection": "grey",
                "containerGraph": False,
                "prefix": "graph_section",
                "label": None
            }
            for key, value in kwargs.items():
                if key in config.keys():
                    config[key] = value
                    
                    #create graph
            if config["containerGraph"]:
                g = config["containerGraph"]
            else:
                g = Digraph()
                
            sectionGraph=g.subgraph(name="cluster_{}".format(config["prefix"]))
            with sectionGraph as sg:
                if config["label"]:
                    section_label = config["label"]
                else:
                    section_label = "Section Graph"
                    if self._graph._name:
                        section_label = "{} {}".format(section_label,self._graph._name)
                path_stats = "{} path{}".format(self._pathNumber, "s" if self._pathNumber>1 else "")
                if self._pathNumber==self._requiredPathNumber and self._optionalPathNumber==0:
                    path_stats = "{}, all required".format(path_stats)
                elif self._requiredPathNumber>0 or self._optionalPathNumber>0:
                    path_stats = "{}, {} required and {} optional".format(path_stats, self._requiredPathNumber, 
                                                                          self._optionalPathNumber)
                section_label = ("<" + 
                                 "<font point-size=\"10\">{}</font><br/>".format(section_label) +
                                 "<font point-size=\"8\">{}</font><br/>".format(path_stats) +
                                 ">")                    
                section_style = config["baseStyleSection"]
                section_fillcolor = config["baseFillColorSection"]
                section_penwidth = config["basePenWidthSection"]
                section_color = config["baseColorSection"]
                sg.attr(label=section_label, style=section_style, fillcolor=section_fillcolor, 
                        color=section_color, penwidth=str(section_penwidth),
                        labelloc="t", nodesep="0", ranksep="0")    

                kwargs["containerGraph"] = sg
                kwargs["prefix"] = config["prefix"]
                kwargs["restrictedListOfOrientatedCkmers"] = self._orientatedCkmers
                self._graph.visualize(*args, **kwargs)
            return g
        
        
        class SectionItemPath:
            
            def __init__(self, orientatedCkmerList, pathSequence, distance, graph):                
                self._orientatedCkmerList = orientatedCkmerList
                self._pathSequence = pathSequence
                self._distance = distance
                self._required = False
                self._optional = False
                self._graph = graph
                
            def distance(self):
                return self._distance
            
            def sequence(self):
                return self._pathSequence
            
            def getOrientatedCkmers(self):
                return self._orientatedCkmerList
            
            def getCandidates(self):
                candidates = []
                for orientatedCkmer in self.getOrientatedCkmers():
                    if self._graph._orientatedCkmers[orientatedCkmer].candidate():
                        candidates.append(orientatedCkmer)
                return set(candidates)
            
            def n(self):
                return len(self._orientatedCkmerList)
            
            def _setRequired(self):
                self._required = True
                self._optional = False
                
                
            def required(self):
                return self._required
            
            def _setOptional(self):
                if not self._required:
                    self._optional = True                
                
            def optional(self):
                return self._optional
                
            
                
                
                
                      
        
            
        
        