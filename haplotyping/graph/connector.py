import logging
from haplotyping.general import General
from haplotyping.graph.api import APIGraph

class ArmConnector:
        
    def __init__(self, graph, **args):
        self._graph = graph
        assert isinstance(self._graph, APIGraph)
        #logger
        self._logger = logging.getLogger(__name__)
        #parameters
        self.numberOfCandidates = args.get("numberOfCandidates",3)
        self.numberOfKmers = args.get("numberOfKmers",10)
        self.maximumLength = args.get("maximumLength",1000)
        self.maximumSteps = args.get("maximumSteps",5)
        self.maximumPaths = args.get("maximumPaths",0)
        self.minimumReadDepth = args.get("minimumReadDepth",2)
        self.minimumGlueSize = args.get("minimumGlueSize",10)

        #main function
        self._graph._detectArms()
        armPairs = self._graph._detectConnectedArmsCandidates()
        for pairCounter in range(len(armPairs)):
            self._logger.debug("trying to connect {} of {} pair(s) of arms".format(pairCounter+1,len(armPairs)))
            arm1 = self._graph.getArm(armPairs[pairCounter][0])
            arm2 = self._graph.getArm(armPairs[pairCounter][1])
            connection = self._graph.createConnection(arm1.connection(),arm2.connection())
            solution,trustedPath1,trustedPath2 = self._connectPair(arm1,arm2)
            if solution:
                self._logger.debug("solution for {} of {} pair(s) of arms".format(pairCounter+1,len(armPairs)))
                sequence,orientatedCkmers = self._createPathSequence(solution, True)
                if sequence:
                    connection.setSequence(sequence,orientatedCkmers)
                else:
                    sequencePath1, orientatedCkmersPath1 = self._createPathSequence(trustedPath1)
                    sequencePath2, orientatedCkmersPath2 = self._createPathSequence(trustedPath2)
                    connection.setSequenceFrom(sequencePath1,orientatedCkmersPath1)
                    connection.setSequenceTo(sequencePath2,orientatedCkmersPath2)
            else:
                sequencePath1, orientatedCkmersPath1 = self._createPathSequence(trustedPath1)
                sequencePath2, orientatedCkmersPath2 = self._createPathSequence(trustedPath2)
                connection.setSequenceFrom(sequencePath1,orientatedCkmersPath1)
                connection.setSequenceTo(sequencePath2,orientatedCkmersPath2)
                self._logger.debug("no solution for {} of {} pair(s) of arms".format(
                    pairCounter+1,len(armPairs)))

    def _createPathSequence(self,path,abortOnProblem=False):
        sequence = path[0][0][0] if path[0][0][1]=="forward" else General.reverse_complement(path[0][0][0])
        orientatedCkmers = [(path[0][0][0],path[0][0][1])]
        k = self._graph._k
        previousPosition = path[2][0]
        previousKmer = sequence
        for i in range(1,len(path[0])):
            position = path[2][i]
            newKmer = path[0][i][0] if path[0][i][1]=="forward" else General.reverse_complement(path[0][i][0])
            d = (position-previousPosition)
            if d<(k-2) and (sequence[position:] == newKmer[0:k-d]):
                sequence = sequence + newKmer[k-d:]
            else:
                data = self._graph._api.getKmerPath(self._graph._uid,previousKmer,newKmer,distance=d+k)
                if data:
                    sequence = sequence + data[k:]
                else:
                    if abortOnProblem:
                        return None,None
                    else:
                        return sequence,orientatedCkmers
            previousPosition = position
            previousKmer = newKmer
            orientatedCkmers.append((path[0][i][0],path[0][i][1]))
        return sequence,orientatedCkmers

    def _processArmReads(self,reads,arm,newArmKmers={}):
        #first raw filtering
        filteredReads = []
        candidateKmers = set([entry[0] for entry in self._graph.getCandidates()])
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
                    if ockmer in self._graph.getCandidates() or ockmer in arm._orientatedCkmers:
                        if orientation is None:
                            orientation="forward"
                        elif not orientation=="forward":
                            orientation = "conflict"
                            continue
                    elif rckmer in self._graph.getCandidates() or rckmer in arm._orientatedCkmers:
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

    def _expandGlueResponses(self,arm,glueResponse):
        for item in glueResponse:
            if not "u" in item[1]:
                if arm.armType()=="outgoing":
                    expanded = True
                    while expanded:
                        expanded = False
                        kmer = item[0][-1]
                        direct = self._graph._api.getSplitDirect(self._graph._uid,kmer[0])
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
                        direct = self._graph._api.getSplitDirect(self._graph._uid,kmer[0])
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

    def _glue(self,reads,arm,initialList=None):
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
                            if glueList[i][3][k]<self.minimumReadDepth:
                                if glueList[i][1][k]=="u":
                                    reportEntry["unconfirmed"] = True
                                elif glueList[i][1][k]=="a":
                                    if arm.armType()=="incoming":
                                        if (not "u" in glueList[i][1]) or (k>glueList[i][1].index("a")+self.numberOfKmers-1):
                                            reportEntry["ignore"] = True
                                    elif arm.armType()=="outgoing":
                                        if (not "u" in glueList[i][1]) or (k<glueList[i][1].index("u")-self.numberOfKmers):
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
    
        self._expandGlueResponses(arm,glueList)                
        
        return glueList

    def _getPathLengths(self,glueResponse):
        if len(glueResponse)>0:
            minLength = min([item[2][-1] for item in glueResponse])
            maxLength = max([item[2][-1] for item in glueResponse])
        else:
            minLength = 0
            maxLength = 0
        length = "{}".format(maxLength) if minLength==maxLength else "{}-{}".format(minLength,maxLength)
        return minLength,maxLength,length

    def _getTrustedPath(self,arm,glueResponse):
        trustedPath = None
        entryIndex = []
        if arm.armType()=="incoming":
            connection = arm.connection()
            for item in glueResponse:
                hasConnection = False
                entry = [[],[],[],[]]
                for i in range(len(item[0])):
                    if item[1][i]=="u":
                        if len(entry[0])>0:
                            break
                        else:
                            continue
                    entry[0].append(item[0][i])
                    entry[1].append(item[1][i])
                    entry[2].append(item[2][i])
                    entry[3].append(item[3][i])
                    if (entry[0][-1][0],entry[0][-1][1]) == connection:
                        hasConnection = True
                        break
                entry[2] = [x - entry[2][0] for x in entry[2]]
                if hasConnection and not entry[0] in entryIndex:
                    entryIndex.append(entry[0])
                    if trustedPath is None or sum(entry[3])>sum(trustedPath[3]):
                        trustedPath = entry
        elif arm.armType()=="outgoing":
            connection = arm.connection()
            for item in glueResponse:
                hasConnection = False
                entry = [[],[],[],[]]
                for i in range(len(item[0])):
                    if item[1][i]=="u":
                        break
                    elif (item[0][i][0],item[0][i][1]) == connection:
                        hasConnection = True
                    elif not hasConnection:
                        continue
                    entry[0].append(item[0][i])
                    entry[1].append(item[1][i])
                    entry[2].append(item[2][i])
                    entry[3].append(item[3][i])
                if hasConnection and not entry[0] in entryIndex:
                    entryIndex.append(entry[0])
                    if trustedPath is None or sum(entry[3])>sum(trustedPath[3]):
                        trustedPath = entry
        return trustedPath

    def _getSolution(self,arm1,arm2,glueResponse1,glueResponse2):
        solution = None
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
                            if p1==p2 and (self.minimumGlueSize==0 or len(p1)>=self.minimumGlueSize):
                                entry = [[],[],[],[]]
                                connection1 = arm1.connection()
                                connection2 = arm2.connection()
                                hasConnection1 = False
                                hasConnection2 = False
                                for m in range(len(path1[0])):
                                    if (path1[0][m][0],path1[0][m][1])==connection1:
                                        hasConnection1=True
                                        entry[0] = path1[0][m:]
                                        entry[1] = path1[1][m:]
                                        entry[2] = path1[2][m:]
                                        entry[3] = path1[3][m:i]
                                        for k in range(i,i+len(p1)):
                                            entry[3].append(max(path1[3][k],path2[3][j-i+k]))
                                        break
                                for m in range(j+len(p1),len(path2[0])):
                                    entry[0].append(path2[0][m])
                                    entry[1].append(path2[1][m])
                                    entry[2].append(path1[2][-1] + path2[2][m] - path2[2][j+len(p1)-1])
                                entry[2] = [x-entry[2][0] for x in entry[2]]
                                entry[3]= entry[3] + path2[3][j+len(p1):]
                                if min(entry[3])>=self.minimumReadDepth:
                                    if solution is None:
                                        solution = entry
                                    elif sum(entry[3])>sum(solution[3]):
                                        solution = entry
                            break
        return solution

    def _expandArm(self, arm, glueResponse, allNewKmersArm, n, processedReads):
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
        newReadsArm = self._graph._api.getSplitReads(self._graph._uid,newKmersArm.keys())
        allNewKmersArm.update(newKmersArm)
        responsenew = self._processArmReads(newReadsArm,arm,allNewKmersArm)
        responsefiltered = []
        for processedRead in responsenew:
            entry = tuple([item[0] for item in processedRead[0]])
            if not entry in processedReads:
                processedReads.add(entry)
                responsefiltered.append(processedRead)
        glueResponse = self._glue(responsefiltered,arm,glueResponse)
        return glueResponse,allNewKmersArm,len(newKmersArm),len(responsefiltered)

    def _connectPair(self, arm1, arm2):
        arm1Kmers = set([item[0] for item in arm1._orientatedCkmers])
        arm2Kmers = set([item[0] for item in arm2._orientatedCkmers])
        assert arm1.armType()=="outgoing"
        assert arm2.armType()=="incoming"
        #get orientated k-mers
        candidatesArm1 = set([entry for entry in [arm1.connection()] 
                                  if entry in self._graph.getCandidates()])
        candidatesArm2 = set([entry for entry in [arm2.connection()] 
                                  if entry in self._graph.getCandidates()])
        for i in range(self.numberOfCandidates):
            newCandidatesArm1 = set()
            for entry in candidatesArm1:
                newCandidatesArm1.update([newEntry for newEntry in self._graph._orientatedCkmers[entry]._incoming
                                          if newEntry in self._graph.getCandidates()])
            candidatesArm1.update(newCandidatesArm1)
            newCandidatesArm2 = set()
            for entry in candidatesArm2:
                newCandidatesArm2.update([newEntry for newEntry in self._graph._orientatedCkmers[entry]._outgoing
                                          if newEntry in self._graph.getCandidates()])
            candidatesArm2.update(newCandidatesArm2)
        #get plain k-mers
        candidateKmersArm1 = set([entry[0] for entry in candidatesArm1])
        candidateKmersArm2 = set([entry[0] for entry in candidatesArm2])

        readsArm1 = self._graph._api.getSplitReads(self._graph._uid,arm1Kmers.union(candidateKmersArm1))
        response1 = self._processArmReads(readsArm1,arm1)
        processedReads1 = set([tuple([item[0] for item in processedRead[0]]) for processedRead in response1])
        glueResponse1 = self._glue(response1,arm1)
        glueResponse1 = [entry for entry in glueResponse1 if "c" in entry[1]]
        for entry in glueResponse1:
            for i in range(len(entry[0])):
                if entry[0][i][0] in arm1Kmers or entry[0][i][0] in candidateKmersArm1:
                    entry[3][i]=max(entry[3][i],self.minimumReadDepth)
        minLength1,maxLength1,length1 = self._getPathLengths(glueResponse1)
        self._logger.debug("starting with {} outgoing path(s) of length {}".format(len(glueResponse1),length1))
        
        readsArm2 = self._graph._api.getSplitReads(self._graph._uid,arm2Kmers.union(candidateKmersArm2))
        response2 = self._processArmReads(readsArm2,arm2)
        processedReads2 = set([tuple([item[0] for item in processedRead[0]]) for processedRead in response2])
        glueResponse2 = self._glue(response2,arm2)
        glueResponse2 = [entry for entry in glueResponse2 if "c" in entry[1]]
        for entry in glueResponse2:
            for i in range(len(entry[0])):
                if entry[0][i][0] in arm2Kmers or entry[0][i][0] in candidateKmersArm2:
                    entry[3][i]=max(entry[3][i],self.minimumReadDepth)
        minLength2,maxLength2,length2 = self._getPathLengths(glueResponse2)
        self._logger.debug("starting with {} incoming path(s) of length {}".format(len(glueResponse2),length2))
        
        allNewKmersArm1 = {}
        allNewKmersArm2 = {}
        counter = 0
        solutions = []
        while True:
            counter+=1
            glueResponse1,allNewKmersArm1,newKmersArm1,newReads1 = self._expandArm(
                arm1,glueResponse1,allNewKmersArm1,self.numberOfKmers,processedReads1)
            minLength1,maxLength1,length1 = self._getPathLengths(glueResponse1)
            self._logger.debug("step {}: expand {} k-mers with {} reads to {} outgoing path(s) of length {}".format(
                counter,newKmersArm1,newReads1,len(glueResponse1),length1))
            glueResponse2,allNewKmersArm2,newKmersArm2,newReads2 = self._expandArm(
                arm2,glueResponse2,allNewKmersArm2,self.numberOfKmers,processedReads2)
            minLength2,maxLength2,length2 = self._getPathLengths(glueResponse2)
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
            solution = self._getSolution(arm1,arm2,glueResponse1,glueResponse2)
            #check for solutions
            if solution:
                break
            #check limits
            if self.maximumSteps>0 and counter>=self.maximumSteps:
                break
            elif self.maximumLength>0 and maxLength1+maxLength2>self.maximumLength:
                break
            elif self.maximumPaths>0 and len(glueResponse1)+len(glueResponse2)>self.maximumPaths:
                break
        trustedPath1 = self._getTrustedPath(arm1,glueResponse1)
        trustedPath2 = self._getTrustedPath(arm2,glueResponse2)
        return solution,trustedPath1,trustedPath2
