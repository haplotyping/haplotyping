import logging,requests,re
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
        self._sections = []
        
        #create sections
        for connectedCandidates in self._graph.getConnectedCandidates():
            self._createSections(connectedCandidates)
            
            
    def _createSections(self, connectedCandidates):
        sectionsItem = self.SectionsItem(self._graph)   
        sectionCkmers = set()
        sectionItem = self.SectionItem(connectedCandidates["start"], connectedCandidates["end"],
                                       connectedCandidates["connected"],self._graph)
        sectionsItem._append(sectionItem)
        sectionCkmers.update(sectionItem._orientatedCkmers)
        while(len(sectionItem._end["ckmers"].intersection(connectedCandidates["end"]))==0):
            newStartCkmers = sectionItem._end["ckmers"].intersection(connectedCandidates["connected"])
            sectionItem = self.SectionItem(newStartCkmers, connectedCandidates["end"],
                                       connectedCandidates["connected"],self._graph)
            sectionsItem._append(sectionItem)
            sectionCkmers.update(sectionItem._orientatedCkmers)
        sectionsItem._computePaths()
        assert len(connectedCandidates["connected"].difference(sectionCkmers))==0
        self._sections.append(sectionsItem)
        self._sections = sorted(self._sections, key=lambda x: 0 if x._order==None else x._order)
                       
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
            "baseStyleSections": "filled",
            "baseFillColorSections": "beige",
            "basePenWidthSections": 1,
            "baseColorSections": "grey",
            "baseStyleSection": "filled",
            "baseFillColorSection": "grey98",
            "basePenWidthSection": 1,
            "baseColorSection": "grey"
        }
        for key, value in kwargs.items():
            if key in config.keys():
                config[key] = value
        
        g = Digraph(name="cluster_graph")
        graph_label = "Sectioned Graph"
        if self._graph._name:
            graph_label = "{} '{}'".format(graph_label,self._graph._name)
        if len(self._sections)>1:
            graph_label = "{} ({} parts)".format(graph_label,len(self._sections))
        graph_label = ("<" + "<font point-size=\"12\" color=\"grey\">" 
                       + graph_label.replace("<", "").replace(">", "") 
                       + "</font>" + ">")
        g.attr(label=graph_label, labelloc="t", nodesep="0", ranksep="0")        
        
        for i in range(len(self._sections)):
            sectionsGraph=g.subgraph(name="cluster_graph_sections_{}".format(i))
            with sectionsGraph as sg:
                sections_label = "{} sections".format(self._sections[i].number())
                if len(self._sections)>1:
                    sections_label = "Part {}, {}".format(i+1,sections_label)
                sections_label = ("<" + 
                                  "<font point-size=\"12\" color=\"grey\">{}</font><br/>".format(sections_label) +
                                  "<font point-size=\"10\">{:,} paths</font><br/>".format(self._sections[i]._pathNumber) +
                                  ">")
                sections_style = config["baseStyleSections"]
                sections_fillcolor = config["baseFillColorSections"]
                sections_penwidth = config["basePenWidthSections"]
                sections_color = config["baseColorSections"]
                sg.attr(label=sections_label, style=sections_style, fillcolor=sections_fillcolor, 
                        color=sections_color, penwidth=str(sections_penwidth), labelloc="t", nodesep="0", ranksep="0")
                #arms
                orientatedCkmersSections = self._sections[i].orientatedCkmers()
                for j in range(len(self._graph._arms)):
                    if self._graph._arms[j].connection() in orientatedCkmersSections:
                        armGraph=sg.subgraph(name="cluster_graph_sections_{}_arm_{}".format(i,self._graph._arms[j].key()))
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
                for j in range(self._sections[i].number()):
                    sectionItem = self._sections[i].get(j)
                    sectionGraph=sg.subgraph(name="cluster_graph_{}_section_{}".format(i,j))
                    with sectionGraph as ssg:
                        section_prefix = "graph_{}_section_{}".format(i,j)
                        section_label = "Section {} of {}".format(j+1,self._sections[i].number())
                        section_label = ("<" + 
                                         "<font point-size=\"10\">{}</font><br/>".format(section_label) +
                                         "<font point-size=\"8\">{} paths</font><br/>".format(len(sectionItem._paths)) +
                                         ">")
                        section_style = config["baseStyleSection"]
                        section_fillcolor = config["baseFillColorSection"]
                        section_penwidth = config["basePenWidthSection"]
                        section_color = config["baseColorSection"]
                        ssg.attr(label=section_label, style=section_style, fillcolor=section_fillcolor, 
                                color=section_color, penwidth=str(section_penwidth),
                                labelloc="t", nodesep="0", ranksep="0")
                        sharedStart = (self._sections[i].get(j-1)._orientatedCkmers.intersection(sectionItem._orientatedCkmers) 
                                       if j>0 else set())
                        sharedEnd = (self._sections[i].get(j+1)._orientatedCkmers.intersection(sectionItem._orientatedCkmers) 
                                       if j<(self._sections[i].number()-1) else set())
                        kwargs["containerGraph"] = ssg
                        kwargs["prefix"] = section_prefix
                        kwargs["hideDeadEndBefore"] = sharedStart
                        kwargs["hideDeadEndAfter"] = sharedEnd
                        #create section graph
                        self._sections[i].get(j).visualize(*args, **kwargs)
                        #detect arm connections
                        for orientatedCkmer in sectionItem._orientatedCkmers:
                            if self._graph._orientatedCkmers[orientatedCkmer].arm():
                                for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                    ckmer_key = "{}_{}_{}_{}_{}".format(section_prefix,
                                                orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                    if self._graph._orientatedCkmers[orientatedCkmer].incomingArmType():
                                        arm_key = "arm_incoming_{}".format(
                                            self._graph._orientatedCkmers[orientatedCkmer].armKey())
                                        if not orientatedCkmer in sharedStart:
                                            g.edge(arm_key,ckmer_key,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                                    elif self._graph._orientatedCkmers[orientatedCkmer].outgoingArmType():
                                        arm_key = "arm_outgoing_{}".format(
                                            self._graph._orientatedCkmers[orientatedCkmer].armKey())
                                        if not orientatedCkmer in sharedEnd:
                                            g.edge(ckmer_key,arm_key,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                                    break
                            else:
                                incomingArm = self._graph._orientatedCkmers[orientatedCkmer].incomingArm()
                                if incomingArm:
                                    for orientatedBase in \
                                             self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                        ckmer_key = "{}_{}_{}_{}_{}".format(section_prefix,
                                                orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                        arm_key = "arm_incoming_{}".format(incomingArm.key())
                                        if (not orientatedCkmer in sharedStart and 
                                            not len(incomingArm._orientatedCkmers.intersection(
                                                orientatedCkmersSections))>0):
                                            g.edge(arm_key,ckmer_key,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                                            break
                                outgoingArm = self._graph._orientatedCkmers[orientatedCkmer].outgoingArm()
                                if outgoingArm:
                                    for orientatedBase in \
                                             self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                        ckmer_key = "{}_{}_{}_{}_{}".format(section_prefix,
                                                orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                        arm_key = "arm_outgoing_{}".format(outgoingArm.key())
                                        if (not orientatedCkmer in sharedEnd and 
                                            not len(outgoingArm._orientatedCkmers.intersection(
                                                orientatedCkmersSections))>0):
                                            g.edge(ckmer_key,arm_key,style="dashed", 
                                                   color="grey", rankdir="lr", constraint="true")
                                            break
                        #link start/end sections
                        for orientatedCkmer in sharedStart:
                            for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                                key1 = "{}_{}_{}_{}_{}".format(previousPrefix,
                                        orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                key2 = "{}_{}_{}_{}_{}".format(section_prefix,
                                        orientatedBase[0],orientatedBase[1],orientatedCkmer[0],orientatedCkmer[1])
                                g.edge(key1,key2,style="dashed", color="grey", rankdir="tb", constraint="true")
                        previousPrefix = section_prefix
        return g
                       
    
    class SectionsItem():
        
        def __init__(self, graph):
            """
            construct sections item
            """
            #logger
            self._logger = logging.getLogger(__name__)
            
            self._graph = graph
            self._sectionItems = []
            self._order : int = None
            self._pathNumber : int = None
            
        def _append(self, item):
            self._sectionItems.append(item)
            if not item._order == None:
                self._order = (item._order if self._order==None else min(self._order,item._order))
                
        def _computePaths(self):
            numbers = {orientatedCkmer: 1 for orientatedCkmer in self._sectionItems[0]._start["ckmers"]}
            for sectionItem in self._sectionItems:
                newNumbers = {orientatedCkmer:0 for orientatedCkmer in sectionItem._end["ckmers"]}
                for sectionItemPath in sectionItem._paths:
                    orientatedCkmers = sectionItemPath.orientatedCkmers()
                    newNumbers[orientatedCkmers[-1]] = (newNumbers[orientatedCkmers[-1]] 
                                                               + numbers[orientatedCkmers[0]])
                numbers = newNumbers
            self._pathNumber = sum(numbers.values())
            
        def number(self):
            return len(self._sectionItems)
        
        def get(self, i):
            assert i>=0 and i<self.number()
            return self._sectionItems[i]
        
        def orientatedCkmers(self):
            orientatedCkmers = set()
            for sectionItem in self._sectionItems:
                orientatedCkmers.update(sectionItem._orientatedCkmers)
            return orientatedCkmers
            
        
    class SectionItem():
    
        def __init__(self, startCkmers, endCkmers, connectedCkmers, graph):
            """
            construct section item
            """
            #logger
            self._logger = logging.getLogger(__name__)
            
            self._start = {"ckmers":set(), "bases": set()}
            self._end = {"ckmers":set(), "bases": set()}
            self._orientatedCkmers = set()
            self._orientatedBases = set()
            
            self._connectedCkmers = connectedCkmers
            self._endCkmers = endCkmers
            self._graph = graph
            self._order : int = None
            self._paths : None
            self._pathNumber : int = None
            self._optionalPathNumber : int = None
            self._requiredPathNumber : int = None            
            
            self._checkBackwardList = set()
            self._checkForwardList = {}           
            
            self._logger.debug("create section with {} starting k-mers".format(len(startCkmers)))
            for startCkmer in startCkmers:
                assert startCkmer in graph._orientatedCkmers
                self._addOrientatedCkmer(startCkmer,start=True)
                assert startCkmer in self._orientatedCkmers
                assert startCkmer in self._start["ckmers"]

            while len(self._checkBackwardList)>0 or len(self._checkForwardList)>0:
            
                self._checkBackwardList = self._checkBackwardList.difference(self._orientatedCkmers)
                while len(self._checkBackwardList)>0:
                    checkBackwardList = self._checkBackwardList.intersection(self._connectedCkmers)
                    self._checkBackwardList = set()
                    self._logger.debug("check backward {}".format(len(checkBackwardList)))
                    for checkBackwardItem in checkBackwardList:
                        self._addOrientatedCkmer(checkBackwardItem)
                        
                #check if section finished
                if len(self._endCkmers.intersection(self._orientatedCkmers))==0:
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
                    
                
        def _addOrientatedCkmer(self, orientatedCkmer, start=False):
            assert orientatedCkmer in self._connectedCkmers
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
                            if outgoingCkmer in self._connectedCkmers:
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
            partialPaths = [{"distance": 0,
                             "sequence": (orientatedCkmer[0] if orientatedCkmer[1]=="forward" 
                                      else General.reverse_complement(orientatedCkmer[0])),
                             "list": [orientatedCkmer]} for orientatedCkmer in 
                                        self._start["ckmers"].intersection(self._connectedCkmers)]
            while len(partialPaths)>0:
                newPartialPaths = []
                for partialPath in partialPaths:
                    if partialPath["list"][-1] in self._end["ckmers"]:
                            self._paths.append(self.SectionItemPath(partialPath["list"],
                                                               partialPath["sequence"],partialPath["distance"]))
                    else:
                        orientatedCkmerInfo = self._graph._orientatedCkmers[partialPath["list"][-1]]
                        for orientatedCkmer in orientatedCkmerInfo._outgoing:
                            if orientatedCkmer in self._orientatedCkmers and orientatedCkmer in self._connectedCkmers:
                                outgoingInfo = orientatedCkmerInfo._outgoing[orientatedCkmer]
                                assert partialPath["sequence"][-self._graph._k:]==outgoingInfo["path"][:self._graph._k]
                                newPartialPath = {"distance": partialPath["distance"] + outgoingInfo["distance"],
                                                  "sequence": partialPath["sequence"] + outgoingInfo["path"][self._graph._k:],
                                                  "list": partialPath["list"] + [orientatedCkmer]}
                                newPartialPaths.append(newPartialPath)
                partialPaths = newPartialPaths 
            #update statistics
            self._pathNumber = len(self._paths)
            self._requiredPathNumber = 0
            self._optionalPathNumber = 0
            frequencies = {orientatedCkmer:[] for orientatedCkmer in 
                           self._connectedCkmers.intersection(self._orientatedCkmers)}
            for i in range(len(self._paths)):
                for orientatedCkmer in self._paths[i].orientatedCkmers():
                    frequencies[orientatedCkmer].append(i)
            coveredOrientatedCkmers = set()
            for orientatedCkmer in frequencies:
                if len(frequencies[orientatedCkmer])==1:
                    self._paths[frequencies[orientatedCkmer][0]]._setRequired()
                    self._requiredPathNumber+=1
                    coveredOrientatedCkmers.update(self._paths[i].orientatedCkmers())
            uncoveredOrientatedCkmers = self._connectedCkmers.difference(coveredOrientatedCkmers)
            if len(uncoveredOrientatedCkmers)==0:
                for i in range(len(self._paths)):
                    if not self._paths[i].required():
                        self._paths[i]._setOptional()
                        self._optionalPathNumber+=1
        
        def visualize(self, *args, **kwargs):            
            kwargs["restrictedListOfCkmers"] = self._orientatedCkmers
            return self._graph.visualize(*args, **kwargs)
        
        
        class SectionItemPath:
            
            def __init__(self, orientatedCkmerList, pathSequence, distance):                
                self._orientatedCkmerList = orientatedCkmerList
                self._pathSequence = pathSequence
                self._distance = distance
                
            def distance(self):
                return self._distance
            
            def sequence(self):
                return self._sequence
            
            def orientatedCkmers(self):
                return self._orientatedCkmerList
            
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
                
            
                
                
                
                      
        
            
        
        