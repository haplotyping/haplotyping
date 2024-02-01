import logging,html
import haplotyping.graph.baseGraph as baseGraph
from haplotyping.graph.api import APIGraph
from haplotyping.general import General
import haplotyping
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
        
        self._datasetFrequencies = None
        self._datasetVarieties = None
        
        #create sections
        for baseConnectedCandidates in self._graph.getConnectedCandidates(True):
            self._createSections(baseConnectedCandidates)
            
            
    def getOrientatedCkmers(self):
        orientatedCkmers = set()
        for i in range(self.partNumber()):
            orientatedCkmers.update(self.part(i).getOrientatedCkmers())
        return orientatedCkmers
        
    def getCandidates(self):
        candidates = []
        for orientatedCkmer in self.getOrientatedCkmers():
            if self._graph._orientatedCkmers[orientatedCkmer].candidate():
                candidates.append(orientatedCkmer)
        return set(candidates)
    
    def processDatasetFrequencies(self):
        if isinstance(self._graph, APIGraph):
            self._datasetFrequencies = self._graph.getDatasetFrequencies()
            self._datasetVarieties = self._graph.getDatasetVarieties()      
            for part in self._parts:
                part._processDatasetFrequencies()                
        else:
            self._logger.warning("unsupported for this graph type")
        
    def _createSections(self, connectedCandidates):
        #compute section start entries (not every connected start is necessarily a section start)
        connected = self._graph.getConnected()
        startEntries = connectedCandidates["start"]
        problematicStartEntries = set()
        for orientatedCkmer in startEntries:
            baseLinkedOrientatedCkmers = set()
            for orientatedBase in self._graph._orientatedCkmers[orientatedCkmer]._orientatedBases.values():
                baseLinkedOrientatedCkmers.update(self._graph._orientatedBases[orientatedBase]._orientatedCkmers)
            for otherOrientatedCkmer in baseLinkedOrientatedCkmers:
                otherIncomingConnected = [k for k in connectedCandidates["connected"] 
                                          if connected[otherOrientatedCkmer][k]]
                if len(otherIncomingConnected)>1:
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
        part = self.SectionsPart(self)   
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
        self._parts = sorted(self._parts, key=lambda x: 0 if x._order is None else x._order)
                       
    def partNumber(self):
        return len(self._parts)
    
    def part(self, i=0):
        if i>=0 and i<len(self._parts):
            return self._parts[i]
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
                
        edges = set()
        nodes = {}

        partConfig = kwargs.copy()
        for i in range(len(self._parts)):
            partConfig["containerGraph"] = g
            partConfig["prefix"] = "graph_part_{}".format(i)
            partConfig["label"] = "Part {}".format(i)
            newEdges, newNodes = self._parts[i].visualize(*args, **partConfig)
            edges.update(newEdges)
            for key,value in newNodes.items():
                if key in nodes:
                    for name,nodeKey in newNodes[key].items():
                        nodes[key][name] = nodeKey
                else:
                    nodes[key] = value.copy()
        return g
                       
    
    class SectionsPart():
        
        def __init__(self, sections):
            """
            construct part
            """
            #logger
            self._logger = logging.getLogger(__name__)
            
            self._sections = sections
            self._graph = sections._graph
            self._sectionItems = []
            self._order : int = None
            self._pathNumber : int = None
                
            self._datasets = set()
            self._varieties = set()
                
        def __repr__(self):
            text = "Part graph"
            if self._graph._name:
                text = "{} {}".format(text,self._graph._name)
            text = "{}; containing {} section{} with {} path{}".format(text,len(self._sectionItems),
                   "s" if len(self._sectionItems)>1 else "", self._pathNumber, "s" if self._pathNumber>1 else "")
            return text
            
        def _append(self, item):
            self._sectionItems.append(item)
            if not item._order is None:
                self._order = (item._order if self._order is None else min(self._order,item._order))
                
        def _computePaths(self):
            self._pathNumber = 0
            numbers = {sectionItemPath.getOrientatedCkmers()[0]: 1 for sectionItemPath in self._sectionItems[0]._paths}
            for sectionItem in self._sectionItems:
                extended = set()
                newNumbers = {sectionItemPath.getOrientatedCkmers()[-1]: 0 for sectionItemPath in sectionItem._paths}
                for sectionItemPath in sectionItem._paths:            
                    orientatedCkmers = sectionItemPath.getOrientatedCkmers()
                    if orientatedCkmers[0] in numbers.keys():
                        extended.add(orientatedCkmers[0])
                    newNumbers[orientatedCkmers[-1]] = (newNumbers[orientatedCkmers[-1]] 
                                                               + numbers.get(orientatedCkmers[0],1))
                for orientatedCkmer in numbers:
                    if not orientatedCkmer in extended:
                        self._pathNumber += numbers[orientatedCkmer]
                numbers = newNumbers
            #add all remaining   
            for orientatedCkmer in numbers:
                if orientatedCkmer in sectionItem._end["ckmers"]:
                    self._pathNumber += numbers[orientatedCkmer]         
            
        def _processDatasetFrequencies(self, ignoreDatasets=[]):
            for sectionItem in self._sectionItems:
                for sectionItemPath in sectionItem._paths:
                    orientatedCkmers = sectionItemPath.getOrientatedCkmers()
                    sectionItemPath._datasets = set(self._sections._datasetFrequencies.index)
                    for orientatedCkmer in orientatedCkmers:
                        ckmerDatasets=self._sections._datasetFrequencies.index[
                            self._sections._datasetFrequencies[orientatedCkmer[0]]>0]
                        sectionItemPath._datasets = sectionItemPath._datasets.intersection(ckmerDatasets)
                    sectionItemPath._varieties = set([self._sections._datasetVarieties[ds]["uid"] 
                                                      for ds in sectionItemPath._datasets 
                                                      if ds in self._sections._datasetVarieties])
                    
        def sectionNumber(self):
            return len(self._sectionItems)
        
        def section(self, i):
            if i>=0 and i<self.sectionNumber():
                return self._sectionItems[i]

        def pathNumber(self):
            return self._pathNumber
            
        def paths(self):
            return self.SectionsPartPathIterator(self)
        
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

            def connectArm(arm_key, arm, **kwargs):
                #create and connect arm            
                orientatedCkmer = arm.connection()
                orientatedCkmersArm = arm._orientatedCkmers
                if orientatedCkmer in nodes or len(orientatedCkmersArm.intersection(nodes))>0:
                    if len(orientatedCkmersArm.intersection(nodes))>0:
                        for orientatedCkmerArm in orientatedCkmersArm.intersection(nodes):
                            if arm.armType()=="incoming":
                                if(len(set(self._graph._orientatedCkmers[orientatedCkmerArm]._incoming.keys())
                                               .intersection(nodes))>0):
                                    continue
                                if orientatedCkmerArm in endNodes:
                                    pg.edge(arm_key,endNodes[orientatedCkmerArm],style="dashed", 
                                                       color="grey", rankdir="lr", constraint="true")
                                else:
                                    for ckmerKey in nodes[orientatedCkmerArm].values():
                                        pg.edge(arm_key,ckmerKey,style="dashed", 
                                                       color="grey", rankdir="lr", constraint="true")
                            elif arm.armType()=="outgoing":
                                if(len(set(self._graph._orientatedCkmers[orientatedCkmerArm]._outgoing.keys())
                                               .intersection(nodes))>0):
                                    continue
                                if orientatedCkmerArm in startNodes:
                                    pg.edge(startNodes[orientatedCkmerArm],arm_key,style="dashed", 
                                                       color="grey", rankdir="lr", constraint="true")
                                else:
                                    for ckmerKey in nodes[orientatedCkmerArm].values():
                                        pg.edge(ckmerKey,arm_key,style="dashed", 
                                                       color="grey", rankdir="lr", constraint="true")
                    else:
                        for ckmerKey in nodes[orientatedCkmer].values():
                            if arm.armType()=="incoming":
                                pg.edge(arm_key,ckmerKey,style="dashed", 
                                               color="grey", rankdir="lr", constraint="true")
                            elif arm.armType()=="outgoing":
                                pg.edge(ckmerKey,arm_key,style="dashed", 
                                               color="grey", rankdir="lr", constraint="true")
 
            def getNodeKeys(arm,nodes,startNodes,endNodes):
                keyList = []
                orientatedCkmer = arm.connection()
                orientatedCkmersArm = arm._orientatedCkmers
                if orientatedCkmer in nodes or len(orientatedCkmersArm.intersection(nodes))>0:
                    if len(orientatedCkmersArm.intersection(nodes))>0:
                        for orientatedCkmerArm in orientatedCkmersArm.intersection(nodes):
                            if arm.armType()=="incoming":
                                if(len(set(self._graph._orientatedCkmers[orientatedCkmerArm]._incoming.keys())
                                               .intersection(nodes))>0):
                                    continue
                                if orientatedCkmerArm in endNodes:
                                    keyList.append(endNodes[orientatedCkmerArm])
                                else:
                                    for ckmerKey in nodes[orientatedCkmerArm].values():
                                        keyList.append(ckmerKey)
                            elif arm.armType()=="outgoing":
                                if(len(set(self._graph._orientatedCkmers[orientatedCkmerArm]._outgoing.keys())
                                               .intersection(nodes))>0):
                                    continue
                                if orientatedCkmerArm in startNodes:
                                    keyList.append(startNodes[orientatedCkmerArm])
                                else:
                                    for ckmerKey in nodes[orientatedCkmerArm].values():
                                        keyList.append(ckmerKey)
                    else:
                        for ckmerKey in nodes[orientatedCkmer].values():
                            if arm.armType()=="incoming":
                                keyList.append(ckmerKey)
                            elif arm.armType()=="outgoing":
                                keyList.append(ckmerKey)
                return keyList
            
            initConfig = {
                "showArms": False,
                "showPotentialConnectedArms": False,
                "baseStylePart": "filled",
                "baseFillColorPart": "beige",
                "basePenWidthPart": 1,
                "baseColorPart": "grey",
                "containerGraph": False,
                "prefix": "graph",
                "label": None
            }
            config = kwargs.copy()
            for key, value in initConfig.items():
                if not key in config.keys():
                    config[key] = value
                    
            #create graph
            if config["containerGraph"]:
                g = config["containerGraph"]
            else:
                g = Digraph()

            edges = set()
            nodes = {}
            startNodes = {}
            endNodes = {}
            
            partGraph=g.subgraph(name="cluster_{}".format(config["prefix"]))      
            
            with partGraph as pg:
                if config["label"]:
                    graph_label = config["label"]
                    graph_label = "{} with {} section{}".format(graph_label, self.sectionNumber(), 
                                                                "s" if self.sectionNumber()>1 else "")
                else:
                    graph_label = "{} section{}".format(self.sectionNumber(),
                                                                "s" if self.sectionNumber()>1 else "")
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
                
                #sections
                previousPrefix = ""
                newEdges = set()
                newNodes = {}
                sectionConfig = config.copy()
                sectionConfig["containerGraph"] = pg
                sectionConfig["showPotentialConnectedArms"] = False
                sectionConfig["showArms"] = False
                sectionConfig["showAllBases"] = False
                sectionConfig["showAllNodes"] = False
                for j in range(self.sectionNumber()):
                    sectionItem = self.section(j)
                    section_prefix = "graph_section_{}".format(j)
                    section_label = "Section {}".format(j)
                    sharedStart = (self._sectionItems[j-1]._orientatedCkmers.intersection(sectionItem._orientatedCkmers) 
                                   if j>0 else set())
                    sharedEnd = (self._sectionItems[j+1]._orientatedCkmers.intersection(sectionItem._orientatedCkmers) 
                                   if j<(self.sectionNumber()-2) else set())
                    sectionConfig["prefix"] = section_prefix
                    sectionConfig["label"] = section_label
                    sectionConfig["hideDeadEndBefore"] = sharedStart
                    sectionConfig["hideDeadEndAfter"] = sharedEnd
                    #create section graph
                    previousEdges = newEdges.copy()
                    previousNodes = newNodes.copy()
                    newEdges, newNodes = self._sectionItems[j].visualize(*args, **sectionConfig)
                    edges.update(newEdges)
                    for key,value in newNodes.items():
                        if key in nodes:
                            for name,nodeKey in newNodes[key].items():
                                nodes[key][name] = nodeKey
                        else:
                            nodes[key] = value.copy()
                    #update start and end
                    for key,value in newNodes.items():
                        if key in sharedStart:
                            for name,nodeKey in newNodes[key].items():
                                startNodes[key] = nodeKey
                        if key in sharedEnd:
                            for name,nodeKey in newNodes[key].items():
                                endNodes[key] = nodeKey
                    #link start/end sections
                    for orientatedCkmer in sharedStart:
                        if orientatedCkmer in previousNodes and orientatedCkmer in newNodes:
                            for ckmerKey1 in previousNodes[orientatedCkmer].values():
                                for ckmerKey2 in newNodes[orientatedCkmer].values():
                                    edge_key = (ckmerKey1, ckmerKey2,)
                                    #single edges
                                    if not edge_key in edges:
                                        edges.add(edge_key)
                                        pg.edge(ckmerKey1, ckmerKey2, style="dashed", color="grey", dir="none", rankdir="tb", constraint="true")
                    previousPrefix = section_prefix                          
                #arms
                if config["showArms"]:
                    arms = self._graph.getArms()
                    processedArms = set()
                    orientatedCkmersSections = self.getOrientatedCkmers()
                    #show potential connected arms
                    if config["showPotentialConnectedArms"]:
                        connectedCandidateArms = self._graph._detectConnectedArmsCandidates()
                        for i in range(len(connectedCandidateArms)):
                            arm1 = self._graph.getArm(connectedCandidateArms[i][0])
                            arm2 = self._graph.getArm(connectedCandidateArms[i][1])
                            #only if visible
                            orientatedCkmer1 = arm1.connection()
                            orientatedCkmer2 = arm2.connection()
                            orientatedCkmersArm1 = arm1._orientatedCkmers
                            orientatedCkmersArm2 = arm2._orientatedCkmers
                            if ((orientatedCkmer1 in nodes or len(orientatedCkmersArm1.intersection(nodes))>0) and
                                (orientatedCkmer2 in nodes or len(orientatedCkmersArm2.intersection(nodes))>0)):
                                processedArms.add(arm1.id())
                                processedArms.add(arm2.id())
                                arm_key1, arm_key2 = self._graph._visualizeConnectedArms(pg,i,arm1,arm2,**config)
                                connectArm(arm_key1,arm1,**config)
                                connectArm(arm_key2,arm2,**config)
                    for j in range(len(arms)):
                        arm = arms[j]
                        #only if no potential connected arms
                        if arm.id() in processedArms:
                            continue
                        else:
                            processedArms.add(arm.id())
                        #only if visible
                        armNodeKeys = getNodeKeys(arm,nodes,startNodes,endNodes)
                        if len(armNodeKeys)>0:
                            arm_key = self._graph._visualizeArm(pg,arm,**config)
                            for armNodeKey in armNodeKeys:
                                if arm.armType()=="incoming":
                                    pg.edge(arm_key,armNodeKey,style="dashed", 
                                                               color="grey", rankdir="lr", constraint="true")
                                elif arm.armType()=="outgoing":
                                    pg.edge(armNodeKey,arm_key,style="dashed", 
                                                               color="grey", rankdir="lr", constraint="true")
                #show potential connected arms
                elif config["showPotentialConnectedArms"]:
                    connectedCandidateArms = self._graph._detectConnectedArmsCandidates()
                    for i in range(len(connectedCandidateArms)):
                        arm1 = self._graph.getArm(connectedCandidateArms[i][0])
                        arm2 = self._graph.getArm(connectedCandidateArms[i][1])
                        orientatedCkmer1 = arm1.connection()
                        orientatedCkmer2 = arm2.connection()
                        orientatedCkmersArm1 = arm1._orientatedCkmers
                        orientatedCkmersArm2 = arm2._orientatedCkmers
                        if ((orientatedCkmer1 in nodes or len(orientatedCkmersArm1.intersection(nodes))>0) and
                            (orientatedCkmer2 in nodes or len(orientatedCkmersArm2.intersection(nodes))>0)):
                            node_keys1 = getNodeKeys(arm1,nodes,startNodes,endNodes)
                            node_keys2 = getNodeKeys(arm2,nodes,startNodes,endNodes)
                            for armNodeKey1 in node_keys1:
                                for armNodeKey2 in node_keys2:
                                    self._graph._visualizeConnectedArmsConnection(pg, armNodeKey1,armNodeKey2, **config)
            if config["containerGraph"]:
                return (list(edges), nodes)
            else:
                return g
            
        
        class SectionsPartPathIterator():
            
            def __init__(self, sectionsPart):
                self._sectionsPart = sectionsPart
                self._n = 0
                self._path = [-1,[]]
                self._starts = []
                #recompute total
                self._sectionsPart._computePaths()
            
            def __iter__(self):
                self._n = 0
                self._path = [-1,[]]
                self._starts = []
                previousEnds = set()
                for j in range(len(self._sectionsPart._sectionItems)):
                    currentEnds = set()
                    for i in range(len(self._sectionsPart._sectionItems[j]._paths)):
                        orientatedCkmers = self._sectionsPart._sectionItems[j]._paths[i].getOrientatedCkmers()
                        if not orientatedCkmers[0] in previousEnds:
                            self._starts.append((j,i))
                        currentEnds.add(orientatedCkmers[-1])
                    previousEnds = currentEnds
                return self
        
            def __next__(self):
                self._computeNextPath()
                orientatedCkmerList = []
                pathSequence = ""
                distance = 0
                for i in range(len(self._path[1])):
                    if self._path[1][i]>=0:
                        path = self._sectionsPart._sectionItems[i]._paths[self._path[1][i]]
                        if len(orientatedCkmerList)>0:
                            newOrientatedCkmerList = path.getOrientatedCkmers()
                            newPathSequence = path.sequence()
                            assert orientatedCkmerList[-1] == newOrientatedCkmerList[0]
                            assert pathSequence[-self._sectionsPart._graph._k:] == newPathSequence[:self._sectionsPart._graph._k]
                            orientatedCkmerList.extend(newOrientatedCkmerList[1:])
                            pathSequence+=newPathSequence[self._sectionsPart._graph._k:]
                        else:
                            orientatedCkmerList.extend(path.getOrientatedCkmers())
                            pathSequence+=path.sequence()
                        distance+=path.distance()
                return self._sectionsPart.SectionPartPath(orientatedCkmerList, pathSequence, distance, self._sectionsPart._graph)
        
            def _computeNextPath(self):
                findNextPath = True
                nItems = len(self._sectionsPart._sectionItems)
                while findNextPath:
                    findNextPath = False
                    if (len(self._path[1])==0) or (self._path[0]==nItems-1):
                        self._path[0]+=1
                        if self._path[0]>=len(self._starts):
                            raise StopIteration
                        else:
                            #initialise first path
                            startItem = self._starts[self._path[0]]
                            self._path[1]=[-1] * nItems
                            self._path[1][startItem[0]] = startItem[1]
                    else:
                        startItem = self._starts[self._path[0]]
                    #only if not already finished
                    if not (startItem[0]==nItems-1):
                        #find last
                        lastItem = None
                        for i in range(startItem[0],nItems):
                            if self._path[1][i]>=0:
                                lastItem = i
                        #go to next path
                        for i in range(lastItem,startItem[0],-1):
                            self._path[1][i]+=1
                            if self._path[1][i]>=len(self._sectionsPart._sectionItems[i]._paths):
                                self._path[1][i]=-1
                            else:
                                break
                        #all paths checked for startitem
                        if (lastItem>startItem[0]) and (self._path[1][startItem[0]+1]<0):
                            #no other paths for this startitem
                            self._path[1] = []
                            findNextPath = True
                        else:
                            #check
                            previousEnd = self._sectionsPart._sectionItems[startItem[0]]._paths[startItem[1]].getOrientatedCkmers()[-1]
                            for i in range(startItem[0]+1,nItems):
                                if self._path[1][i]<0:
                                    self._path[1][i]=0
                                orientatedCkmers = self._sectionsPart._sectionItems[i]._paths[self._path[1][i]].getOrientatedCkmers()
                                if not orientatedCkmers[0]==previousEnd:
                                    for j in range(i+1,nItems):
                                        self._path[1][j] = -1
                                    findNextPath = True
                                    #check if this is end of path, then ok
                                    currentPathStarts = set([path.getOrientatedCkmers()[0] for 
                                                                path in self._sectionsPart._sectionItems[i]._paths])
                                    if not previousEnd in currentPathStarts:
                                        self._path[1][i] = -1
                                        findNextPath = False
                                    break
                                else:
                                    previousEnd = orientatedCkmers[-1]
                            #reset to next start
                            if max(self._path[1][startItem[0]+1:])<0:
                                self._path[1] = []
                                findNextPath = True
        
                self._n += 1
    
        class SectionPartPath():
                
            def __init__(self, orientatedCkmerList, pathSequence, distance, graph):                
                self._orientatedCkmerList = orientatedCkmerList
                self._pathSequence = pathSequence
                self._distance = distance
                self._graph = graph
                
            def distance(self):
                return self._distance
            
            def sequence(self):
                return self._pathSequence
            
            def kmers(self):
                kmers = []
                for i in range(len(self._pathSequence)-self._graph._k+1):
                    kmers.append(self._pathSequence[i:self._graph._k+i])
                return kmers
            
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
            
            self._datasets = set()
            self._varieties = set()
            
            #compute connected bases
            for connectedCkmer in connectedCkmers:
                self._connectedBases.update(self._graph._orientatedCkmers[connectedCkmer]._orientatedBases.values())

            self._logger.debug("create section with {} starting k-mers".format(len(initialCkmers)))
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
                if not self._graph._orientatedCkmers[orientatedCkmer]._order is None:
                    self._order = (self._graph._orientatedCkmers[orientatedCkmer]._order if self._order is None else
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
                    if orientatedBaseCkmer in self._endCkmers:
                        self._end["ckmers"].add(orientatedBaseCkmer)
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
                             "list": [orientatedCkmer]} 
                            for orientatedCkmer in 
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
                            #ignore cycles
                            if orientatedCkmer in partialPath["list"]:                                
                                pass
                            elif orientatedCkmer in self._orientatedCkmers and orientatedCkmer in self._connectedCkmers:
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
            if i>=0 and i<self._pathNumber:
                return self._paths[i]
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
            
            initConfig = {
                "baseStyleSection": "filled",
                "baseFillColorSection": "grey98",
                "basePenWidthSection": 1,
                "baseColorSection": "grey",
                "containerGraph": False,
                "prefix": "graph_section",
                "label": None
            }
            config = kwargs.copy()
            for key, value in initConfig.items():
                if not key in config.keys():
                    config[key] = value

            #create graph
            if config["containerGraph"]:
                g = config["containerGraph"]
            else:
                g = Digraph()
                
            sectionGraph=g.subgraph(name="cluster_{}".format(config["prefix"]))
            sectionConfig = config.copy()
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

                sectionConfig["containerGraph"] = sg
                sectionConfig["restrictedListOfOrientatedCkmers"] = self._orientatedCkmers
                edges,nodes = self._graph.visualize(*args, **sectionConfig)
            if config["containerGraph"]:
                return (list(edges), nodes)
            else:
                return g

        class SectionItemPath():
            
            def __init__(self, orientatedCkmerList, pathSequence, distance, graph):         
                #initialise
                self._orientatedCkmerList = orientatedCkmerList
                self._pathSequence = pathSequence
                self._distance = distance
                self._graph = graph
                #other variables
                self._required = False
                self._optional = False
                self._datasets = set()
                self._varieties = set()

            def distance(self):
                return self._distance
            
            def sequence(self):
                return self._pathSequence

            def kmers(self):
                kmers = []
                for i in range(len(self._pathSequence)-self._graph._k+1):
                    kmers.append(self._pathSequence[i:self._graph._k+i])
                return kmers
            
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
                
            
                
                
                
                      
        
            
        
        