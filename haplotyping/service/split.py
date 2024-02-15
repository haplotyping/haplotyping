import h5py, haplotyping, numpy as np

class Split:
    
    #---
    
    def _findItem(item,table,start=0,number=None, cache={}):
        if number==None:
            number=table.shape[0]
        itemBinary = bytearray(item,"utf8")
        newCache={}
        minRowId = start
        maxRowId = number-1
        for key in cache.keys():
            if cache[key]>=minRowId:
                if key<itemBinary:
                    minRowId=max(minRowId,cache[key])
                elif key>itemBinary:
                    maxRowId=min(maxRowId,cache[key])
                    newCache[key]=cache[key]
                else:
                    minRowId = cache[key]
                    maxRowId = cache[key]
                    newCache[key]=cache[key]
        currentRowId = minRowId + int((maxRowId-minRowId)/2) 
        while True:
            if start>=number:
                return (None,currentRowId,newCache,)
            currentRow = table[currentRowId]
            currentItem = currentRow[0]        
            if currentItem==itemBinary:
                newCache[currentItem] = currentRowId
                return (currentRow,currentRowId,newCache,)
            elif currentItem>itemBinary:
                newCache[currentItem] = currentRowId
                maxRowId = currentRowId-1            
            else:
                minRowId = currentRowId+1
            if maxRowId<minRowId:
                return (None,currentRowId,newCache,)
            else:
                currentRowId = minRowId + int((maxRowId-minRowId)/2)
                
    #old code to get correct entries without index
    def _findId(id,table,start=0,number=None, cache={}):
        if number==None:
            number=table.shape[0]
        newCache={}
        minRowId = start
        maxRowId = number-1
        for cacheId in cache.keys():
            if cache[cacheId]>=minRowId:
                if cacheId<id:
                    minRowId=max(minRowId,cache[cacheId])
                elif cacheId>id:
                    maxRowId=min(maxRowId,cache[cacheId])
                    newCache[cacheId]=cache[cacheId]
                else:
                    minRowId = cache[cacheId]
                    maxRowId = cache[cacheId]
                    newCache[cacheId]=cache[cacheId]
        currentRowId = minRowId + int((maxRowId-minRowId)/2) 
        while True:
            if start>=number:
                return ([],currentRowId,newCache,)
            currentRow = table[currentRowId]
            currentId = currentRow[0][0]    
            if currentId==id:
                currentRows = [currentRow]
                tmpId = currentRowId
                while tmpId>0:
                    tmpId-=1
                    tmpRow = table[tmpId]
                    if tmpRow[0][0]==id:
                        currentRows.append(tmpRow)
                    else:
                        break
                tmpId = currentRowId
                while tmpId<table.shape[0]-1:
                    tmpId+=1
                    tmpRow = table[tmpId]
                    if tmpRow[0][0]==id:
                        currentRows.append(tmpRow)
                    else:
                        break
                newCache[currentId] = currentRowId
                return (currentRows,currentRowId,newCache,)
            elif currentId>id:
                newCache[currentId] = currentRowId
                maxRowId = currentRowId-1            
            else:
                minRowId = currentRowId+1
            if maxRowId<minRowId:
                return ([],currentRowId,newCache,)
            else:
                currentRowId = minRowId + int((maxRowId-minRowId)/2)
                
    #---
    
    def _translate_type(value):
        if value=="l":
            return "left"
        elif value=="r":
            return "right"
        elif value=="b":
            return "both"
        else:
            return value
    
    def _base_result(row, h5file: h5py.File, expandCkmer=True):
        if row:
            response = {
                "base": row[0].decode("ascii"),
                "number": int(row[1])
            }
            branches = row.dtype[2].names
            ckmerTable = h5file.get("/split/ckmer")
            if expandCkmer:
                response["ckmers"] = []
                for i in range(len(branches)):
                    if row[2][i][0]>0:
                        ckmerRow = ckmerTable[row[2][i][1]]
                        response["ckmers"].append(Split._kmer_result(ckmerRow,h5file,False))
            else:                
                response["ckmers"] = {}
                for i in range(len(branches)):
                    if row[2][i][0]>0:
                        kmer = response["base"]+branches[i]
                        ckmer = haplotyping.General.canonical(kmer)
                        ckmerRow = ckmerTable[row[2][i][1]]
                        response["ckmers"][ckmer] = {
                            "number": int(ckmerRow[2]), 
                            "split": Split._translate_type(ckmerRow[1].decode("ascii"))
                        }
        else:
            response = None
        return response
                
    def _kmer_result(row, h5file: h5py.File, expandBase=True):
        if row:
            response = {
                "ckmer": row[0].decode("ascii"),
                "split": Split._translate_type(row[1].decode("ascii")),
                "number": int(row[2]),
                "rightSplitBase": {
                },
                "direct": {
                    "left": {
                        "distinct": int(row[4][1][0]),
                        "number": int(row[4][1][1])
                    },"right": {
                        "distinct": int(row[4][2][0]),
                        "number": int(row[4][2][1])
                    }
                },
                "cycle": int(row[6][0]),
                "reverse": int(row[7][0]),
                "paired": int(row[8][1]),                
            }
            baseTable = h5file.get("/split/base")
            if response["split"] in ["left","both"]:
                baseRow=baseTable[row[3][0]]
                if expandBase:
                    response["rightSplitBase"]["left"] = Split._base_result(baseRow,h5file,False)
                else:
                    response["rightSplitBase"]["left"] = baseRow[0].decode("ascii")
            if response["split"] in ["right","both"]:
                baseRow=baseTable[row[3][1]]
                if expandBase:
                    response["rightSplitBase"]["right"] = Split._base_result(baseRow,h5file,False)
                else:
                    response["rightSplitBase"]["right"] = baseRow[0].decode("ascii")            
        else:
            response = None
        return response
    
    def _kmer_direct_result(ckmerRow, directRows, h5file: h5py.File, expandBase=True):
        if (not ckmerRow==None) or (len(directRows)==0):
            response = {
                "ckmer": ckmerRow[0].decode("ascii"),
                "split": Split._translate_type(ckmerRow[1].decode("ascii")),
                "number": int(ckmerRow[2]),
                "direct": {
                    "left": [],
                    "right": []
                }
            } 
            ckmerTable = h5file.get("/split/ckmer")
            baseTable = h5file.get("/split/base")
            for directRow in directRows:
                fromDirection = directRow[0][1].decode("ascii")
                ckmerRow = ckmerTable[directRow[1][0]]
                item = {
                    "connection": {                    
                        "direction": Split._translate_type(directRow[1][1].decode("ascii")),
                        "distance": int(directRow[3]),
                        "number": int(directRow[2]),
                        "problem": int(directRow[4])
                    },
                    "ckmer": ckmerRow[0].decode("ascii"),
                    "split": Split._translate_type(ckmerRow[1].decode("ascii")),
                    "number": int(ckmerRow[2]),
                    "rightSplitBase": {
                    },                
                }
                if item["split"] in ["left","both"]:
                    baseRow=baseTable[ckmerRow[3][0]]
                    item["rightSplitBase"]["left"] = Split._base_result(baseRow,h5file,False)
                if item["split"] in ["right","both"]:
                    baseRow=baseTable[ckmerRow[3][1]]
                    item["rightSplitBase"]["right"] = Split._base_result(baseRow,h5file,False)
                if fromDirection=="l":
                    response["direct"]["left"].append(item)
                elif fromDirection=="r":
                    response["direct"]["right"].append(item)
        else:
            response = None
        return response

    def _data_kmer(kmerId,h5file,kmerDict={}):
        if not kmerId in kmerDict:
            ckmerTable = h5file.get("/split/ckmer")
            entry = ckmerTable[kmerId]
            return (entry[0].decode("ascii"),Split._translate_type(entry[1].decode("ascii")),
                                int(entry[2]),int(entry[4][0]),int(entry[4][1][0]+entry[4][2][0]))
        return kmerDict[kmerId]
    
    def _data_kmer_direct(kmerId,h5file,kmerDict={},directDict={}):
        if not kmerId in directDict:
            entry = Split._data_kmer(kmerId,h5file,kmerDict)
            directTable = h5file.get("/relations/direct")
            directList = directTable[entry[3]:entry[3]+entry[4]]
            kmerDirectDict = {}
            for directEntry in directList:
                assert directEntry[0][0]==kmerId
                fromSide = directEntry[0][1].decode()
                toKmer = directEntry[1][0]
                toSide = directEntry[1][1].decode()
                kmerDirectDict[fromSide] = kmerDirectDict.get(fromSide,{})
                kmerDirectDict[fromSide][toKmer] = (toSide,directEntry[3],)
            return kmerDirectDict
        return directDict[kmerId]
    
    def _kmer_direction_reverse(direction):
        if direction=="l":
            return "r"
        elif direction=="r":
            return "l"
        else:
            return "u"
    
    def _data_kmer_connect(start,end,h5file,kmerDict={},directDict={}):
        startLength = 0
        endLength = 0
        optionsStart = []
        optionsEnd = []
        startFinished = False
        endFinished = False
        nstart = start.copy()
        nend = end.copy()
        while not (startFinished and endFinished):
            if not startFinished:
                optionsStart = []
                kmerDict[nstart[-2]] = Split._data_kmer(nstart[-2],h5file,kmerDict)
                directDict[nstart[-2]] = Split._data_kmer_direct(nstart[-2],h5file,kmerDict,directDict)   
                if nstart[-1] in directDict[nstart[-2]]:
                    for key,value in directDict[nstart[-2]][nstart[-1]].items():
                        optionsStart.append([value[1],value[0],key,Split._kmer_direction_reverse(value[0])])
                    if not len(optionsStart)==1:
                        startFinished = True
                else:
                    startFinished = True
            if not endFinished:
                optionsEnd = []
                kmerDict[nend[-2]] = Split._data_kmer(nend[-2],h5file,kmerDict)
                directDict[nend[-2]] = Split._data_kmer_direct(nend[-2],h5file,kmerDict,directDict) 
                if nend[-1] in directDict[nend[-2]]:
                    for key,value in directDict[nend[-2]][nend[-1]].items():
                        optionsEnd.append([value[1],value[0],key,Split._kmer_direction_reverse(value[0])])
                    if not len(optionsEnd)==1:
                        endFinished = True
                else:
                    endFinished = True
            #try to find a connection from start options, return on first (should be only) match
            if len(optionsStart)>0:
                for entry in optionsStart:
                    if entry[2]==nend[-2] and entry[1]==nend[-1]:
                        solution = nstart + entry[0:-1] + nend[::-1][2:]
                        return solution, False, False
            #glue start
            if len(optionsStart)==1:
                nstart.extend(optionsStart[0])
                startLength+=optionsStart[0][0]
                optionsStart = []
            #try to find a connection from end options, return on first (should be only) match
            if len(optionsEnd)>0:
                for entry in optionsEnd:
                    if entry[2]==nstart[-2] and entry[1]==nstart[-1]:
                        solution = nstart + [entry[0]] + nend[::-1]
                        return solution,False, False
                if len(optionsStart)>0:
                    for entryStart in optionsStart:
                        for entryEnd in optionsEnd:
                            if entryStart[2]==entryEnd[2] and entryStart[3]==entryEnd[1]:                            
                                solution = nstart + entryStart + [entryEnd[0]] + nend[::-1]
                                return solution,False, False
            #glue end
            if len(optionsEnd)==1:
                nend.extend(optionsEnd[0])
                endLength+=optionsEnd[0][0]
                optionsEnd = []
        #fix returned options
        if len(start)<len(nstart):
            for i in range(len(optionsStart)):
                optionsStart[i] = nstart[len(start):]+optionsStart[i]
        if len(end)<len(nend):
            for i in range(len(optionsEnd)):
                optionsEnd[i] = nend[len(end):]+optionsEnd[i]
        return False, optionsStart, optionsEnd
    
    def _data_kmer_connection(start,end,h5file,kmerDict,directDict,connectionDict):
        if tuple(start) in connectionDict:
            if tuple(end) in connectionDict[tuple(start)]:
                return connectionDict[tuple(start)][tuple(end)]
        elif tuple(end) in connectionDict:
            if tuple(start) in connectionDict[tuple(end)]:
                connection = connectionDict[tuple(end)][tuple(start)]
                if connection:
                    return connection[::-1]
                else:
                    return False
            else:
                connectionDict[tuple(start)] = {}
        else:
            connectionDict[tuple(start)] = {}
        (connection, optionsStart, optionsEnd) = Split._data_kmer_connect(start,end,h5file,kmerDict,directDict)
        if connection:
            connectionDict[tuple(start)][tuple(end)] = connection
            return connection
        elif len(optionsStart)>0 and len(optionsEnd)>0:
            for i in range(len(optionsStart)):
                for j in range(len(optionsEnd)):
                    newStart = optionsStart[i][-2:]
                    newEnd = optionsEnd[j][-2:]
                    (newConnection, newOptionsStart, newOptionsEnd) = Split._data_kmer_connect(
                        newStart,newEnd,h5file,kmerDict,directDict)  
                    if newConnection:
                        connection = start + optionsStart[i][:-2] + newConnection + optionsEnd[j][::-1][2:] + end[::-1]
                        connectionDict[tuple(start)][tuple(end)] = connection
                        return connection
        connectionDict[tuple(start)][tuple(end)] = False
        return False
    
    
    def _kmer_read_result(kmerIds,readInfoList,readData,h5file,kmerDict={},directDict={}):
        problems = 0
        response = []
        connectionDict = {}
        k = int(h5file.get("/config").attrs["k"])
        #get expanded checkset
        expandedKmerIds = set()
        def _expand(id,direction):
            expandedKmerIds.add(id)
            while True:
                kmerDict[id] = Split._data_kmer(id,h5file,kmerDict)
                directDict[id] = Split._data_kmer_direct(id,h5file,kmerDict,directDict)
                if direction in directDict[id]:
                    expandedKmerIds.update(directDict[id][direction].keys())
                    if len(directDict[id][direction])==1:
                        for id,value in directDict[id][direction].items():
                            direction="l" if value[0]=="r" else "r"
                    else:
                        options = []
                        for id,value in directDict[id][direction].items():
                            options.append((id,"l" if value[0]=="r" else "r"))
                        return options
                else:
                    return []
                    
        for kmerId in kmerIds:
            for direction in ["r","l"]:
                id = kmerId
                options = _expand(id,direction)
                for option in options:
                    _expand(option[0],option[1])
        #get reads
        if len(kmerIds)>0:
            n = 0
            checkIds = set(kmerIds)
            for item in readInfoList:
                read = readData[n:n+item[0]]
                number = item[1]
                n+=item[0]
                if not any(x in expandedKmerIds for x in read):
                    continue
                elif len(read)>1:
                    initialConnections = []
                    readConnection = None
                    #more options for initial connection
                    for sd in ["l","r"]:
                        for ed in ["l","r"]:
                            start = [read[0],sd]
                            end = [read[1],ed]
                            connection = Split._data_kmer_connection(start,end,h5file,kmerDict,directDict,connectionDict)
                            if connection and not connection in initialConnections:
                                initialConnections.append(connection)
                    for readConnection in initialConnections:
                        for i in range(2,len(read)):
                            assert readConnection[-1]==read[i-1]
                            for ed in ["l","r"]:
                                start = [read[i-1],Split._kmer_direction_reverse(readConnection[-2])]
                                end = [read[i],ed]
                                connection = Split._data_kmer_connection(start,end,h5file,kmerDict,directDict,connectionDict)
                                if connection:
                                    assert readConnection[-1]==connection[0]
                                    readConnection = readConnection + connection[1:]
                                    break
                            if not connection:
                                readConnection = None
                                break
                        if readConnection:
                            break
                    if not readConnection:
                        problems+=1
                    else:
                        ids = [readConnection[i] for i in range(0,len(readConnection),4)]
                        #only if relevant
                        if any(x in kmerIds for x in ids):
                            kmerList = []
                            position = 0
                            for i in range(0,len(readConnection),4):
                                entry = Split._data_kmer(readConnection[i],h5file,kmerDict)
                                if i>0:
                                    orientation = ("forward" if readConnection[i-1]=="l" else 
                                                   ("backward" if readConnection[i-1]=="r" else "unknown"))
                                    position+=readConnection[i-2]
                                else:
                                    orientation = ("forward" if readConnection[i+1]=="r" else 
                                                   ("backward" if readConnection[i+1]=="l" else "unknown"))
                                kmerInfo = {"ckmer": entry[0], 
                                            "split": entry[1], 
                                            "number": int(entry[2]),
                                            "position": int(position),
                                            "orientation": orientation
                                           }                        
                                kmerList.append(kmerInfo)
                            length = sum([readConnection[i+2] for i in range(0,len(readConnection)-1,4)])+k
                            response.append({"kmers": kmerList, "length": int(length), "number": int(number)})
        return response,problems
    
    def _kmer_paired_result(kmerId,pairedList,h5file,kmerDict={}):
        if len(pairedList)>0:
            ckmerTable = h5file.get("/split/ckmer")
            response = []
            for item in pairedList:
                if item[0]==kmerId:
                    if not item[1] in kmerDict:
                        kmerDict[item[1]] = ckmerTable[item[1]][0].decode("ascii")
                    response.append(kmerDict[item[1]])
        else:
            response = []
        return response,kmerDict
    
    #---
    
    def _info(h5file: h5py.File):
        response = {}
        configTable = h5file.get("/config")
        for k in configTable.attrs.keys():
            value = configTable.attrs[k]
            if np.issubdtype(type(value), np.integer):
                response[k] = int(value)
            else:
                response[k] = str(value)
        return response

    def _kmer_distribution(h5file: h5py.File):
        histogramTable = h5file.get("/histogram/kmer")
        configTable = h5file.get("/config")
        response = {}
        response["k"] = int(configTable.attrs["k"])
        response["minimumKmerFrequencies"] = int(configTable.attrs["minimumKmerFrequencies"])
        response["maximumKmerFrequencies"] = int(configTable.attrs["maximumKmerFrequencies"])
        response["numberKmers"] = int(configTable.attrs["numberKmers"])
        response["totalKmerFrequencies"] = int(configTable.attrs["totalKmerFrequencies"])
        response["frequencies"] = {int(item[0]): int(item[1]) for item in histogramTable}
        return response

    def _kmer_split_distribution(h5file: h5py.File):
        histogramTable = h5file.get("/histogram/ckmer")
        configTable = h5file.get("/config")
        response = {}
        response["k"] = int(configTable.attrs["k"])
        response["minimumCanonicalSplitFrequency"] = int(configTable.attrs["minimumCanonicalSplitFrequency"])
        response["maximumCanonicalSplitFrequency"] = int(configTable.attrs["maximumCanonicalSplitFrequency"])
        response["numberCanonicalSplit"] = int(configTable.attrs["numberCanonicalSplit"])
        response["numberCanonicalSplitBoth"] = int(configTable.attrs["numberCanonicalSplitBoth"])
        response["numberCanonicalSplitLeft"] = int(configTable.attrs["numberCanonicalSplitLeft"])
        response["numberCanonicalSplitRight"] = int(configTable.attrs["numberCanonicalSplitRight"])
        response["totalCanonicalSplitFrequencies"] = int(configTable.attrs["totalCanonicalSplitFrequencies"])
        response["frequencies"] = {int(item[0]): int(item[1]) for item in histogramTable}
        return response

    def _kmer_base_distribution(h5file: h5py.File):
        histogramTable = h5file.get("/histogram/base")
        configTable = h5file.get("/config")
        response = {}
        response["k"] = int(configTable.attrs["k"])
        response["numberRightSplitBases"] = int(configTable.attrs["numberRightSplitBases"])
        response["numberRightSplitKmers"] = int(configTable.attrs["numberRightSplitKmers"])
        response["frequencies"] = {int(item[0]): int(item[1]) for item in histogramTable}
        return response

    def _kmer_split_direct_distribution(h5file: h5py.File):
        histogramTable = h5file.get("/histogram/distance")
        configTable = h5file.get("/config")
        response = {}
        response["k"] = int(configTable.attrs["k"])
        response["maximumCycleLength"] = int(configTable.attrs["maximumCycleLength"])
        response["maximumReversalLength"] = int(configTable.attrs["maximumReversalLength"])
        response["numberCycles"] = int(configTable.attrs["numberCycles"])
        response["numberReversals"] = int(configTable.attrs["numberReversals"])
        response["frequencies"] = {int(item[0]): int(item[1]) for item in histogramTable}
        return response
    
    def _kmer_info(h5file: h5py.File, kmer: str):
        ckmer = haplotyping.General.canonical(kmer)
        ckmerTable = h5file.get("/split/ckmer")
        (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable)
        return Split._kmer_result(ckmerRow,h5file)
    
    def _kmers_info(h5file: h5py.File, kmers: list):
        ckmerList = set()
        for kmer in kmers:
            ckmerList.add(haplotyping.General.canonical(kmer))
        ckmerList = list(ckmerList)
        ckmerList.sort()
        response = []
        ckmerTable = h5file.get("/split/ckmer")
        number = ckmerTable.shape[0]
        start = 0
        cache = {}
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable,start,number,cache)
            if ckmerRow:
                response.append(Split._kmer_result(ckmerRow,h5file))
                start = id+1
            else:
                start = id 
        return response
    
    def _kmer_direct(h5file: h5py.File, kmer: str):
        ckmer = haplotyping.General.canonical(kmer)
        ckmerTable = h5file.get("/split/ckmer")
        (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable)
        if ckmerRow:
            directTable = h5file.get("/relations/direct")
            directRows = directTable[ckmerRow[4][0]:ckmerRow[4][0]+(ckmerRow[4][1][0]+ckmerRow[4][2][0])]
            return Split._kmer_direct_result(ckmerRow,directRows,h5file)
        else:
            return None
    
    def _kmers_direct(h5file: h5py.File, kmers: list):
        ckmerList = set()
        for kmer in kmers:
            ckmerList.add(haplotyping.General.canonical(kmer))
        ckmerList = list(ckmerList)
        ckmerList.sort()
        response = []
        ckmerTable = h5file.get("/split/ckmer")
        directTable = h5file.get("/relations/direct")
        number = ckmerTable.shape[0]
        start = 0
        cache = {}
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable,start,number,cache)
            if ckmerRow:
                directRows = directTable[ckmerRow[4][0]:ckmerRow[4][0]+(ckmerRow[4][1][0]+ckmerRow[4][2][0])]
                response.append(Split._kmer_direct_result(ckmerRow,directRows,h5file))
                start = id+1
            else:
                start = id 
        return list(filter(None, response))
   
    def _kmer_read(h5file: h5py.File, kmer: str):
        ckmer = haplotyping.General.canonical(kmer)
        ckmerTable = h5file.get("/split/ckmer")
        readPartitionTable = h5file.get("/relations/readPartition")
        readDataTable = h5file.get("/relations/readData")
        readInfoTable = h5file.get("/relations/readInfo")
        (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable)
        if ckmerRow:
            kmerDict = {}
            directDict = {}
            partitionRow = readPartitionTable[ckmerRow[5]]
            readDataList = readDataTable[partitionRow[0][0]:partitionRow[0][0]+partitionRow[0][1]]
            readInfoList = readInfoTable[partitionRow[1][0]:partitionRow[1][0]+partitionRow[1][1]]
            reads,problems = Split._kmer_read_result([id],readInfoList,readDataList,h5file)
        else:
            reads = []
        return reads
        
    def _kmers_read(h5file: h5py.File, kmers: list, additional: list):
        ckmerList = set()
        ckmerSet = set()
        for kmer in kmers:
            ckmerSet.add(haplotyping.General.canonical(kmer))
        ckmerList.update(ckmerSet)
        for kmer in additional:
            ckmerList.add(haplotyping.General.canonical(kmer))
        ckmerList = list(ckmerList)
        ckmerList.sort()
        ckmerTable = h5file.get("/split/ckmer")
        readPartitionTable = h5file.get("/relations/readPartition")
        readDataTable = h5file.get("/relations/readData")
        readInfoTable = h5file.get("/relations/readInfo")
        number = ckmerTable.shape[0]
        start = 0
        cache = {}
        partitions = set()
        kmerIds = []
        kmerDict = {}
        directDict = {}
        #get partition data from k-mers
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmerList[i],ckmerTable,start,number,cache)
            if ckmerRow:
                if ckmerList[i] in ckmerSet:
                    kmerIds.append(id)
                partitions.add(ckmerRow[5])
                start = id+1
            else:
                start = id 
        #collect reads from partitions
        reads = []
        problems = 0
        for p in partitions:
            partitionRow = readPartitionTable[p]
            readDataList = readDataTable[partitionRow[0][0]:partitionRow[0][0]+partitionRow[0][1]]
            readInfoList = readInfoTable[partitionRow[1][0]:partitionRow[1][0]+partitionRow[1][1]]
            newReads,newProblems = Split._kmer_read_result(kmerIds,readInfoList,readDataList,h5file)
            problems+=newProblems
            reads.extend(newReads)
        return reads
    
    def _kmer_paired(h5file: h5py.File, kmer: str):
        ckmer = haplotyping.General.canonical(kmer)
        ckmerTable = h5file.get("/split/ckmer")
        pairedTable = h5file.get("/relations/paired")
        (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable)
        if ckmerRow:
            kmerDict = {id: ckmerRow[0].decode("ascii")}
            pairedList = pairedTable[ckmerRow[8][0]:ckmerRow[8][0]+ckmerRow[8][1]]
            paired,kmerDict = Split._kmer_paired_result(id,pairedList,h5file,kmerDict)
        else:
            paired = []
        return paired
        
    def _kmers_paired(h5file: h5py.File, kmers: list):
        ckmerList = set()
        for kmer in kmers:
            ckmerList.add(haplotyping.General.canonical(kmer))
        ckmerList = list(ckmerList)
        ckmerList.sort()
        ckmerTable = h5file.get("/split/ckmer")
        pairedTable = h5file.get("/relations/paired")
        number = ckmerTable.shape[0]
        start = 0
        cache = {}
        kmerDict = {}
        response = {}
        #get paired data for k-mers
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmerList[i],ckmerTable,start,number,cache)
            if ckmerRow:
                kmerDict[id] = ckmerRow[0].decode("ascii")
                if ckmerRow[8][1]>0:
                    pairedList = pairedTable[ckmerRow[8][0]:ckmerRow[8][0]+ckmerRow[8][1]]
                    response[kmerDict[id]],kmerDict = Split._kmer_paired_result(id,pairedList,h5file,kmerDict)
                start = id+1
            else:
                start = id 
        return response
    
    #---
    
    def info(location_split: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._info(h5file)
    
   #---
    
    def kmer_distribution(location_split: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_distribution(h5file)

    def kmer_split_distribution(location_split: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_split_distribution(h5file)

    def kmer_base_distribution(location_split: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_base_distribution(h5file)

    def kmer_split_direct_distribution(location_split: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_split_direct_distribution(h5file)
    
   #---
    
    def kmer_info(location_split: str, kmer: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_info(h5file,kmer)
    
    def kmer_list_info(location_split: str, kmers: list):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmers_info(h5file,kmers)
    
    def kmer_sequence_info(location_split: str, sequence: str):
        with h5py.File(location_split, mode="r") as h5file:            
            configTable = h5file.get("/config")
            k = int(configTable.attrs["k"])
            kmers = [sequence[i:i+k] for i in range(len(sequence)-(k-1))] 
            return Split._kmers_info(h5file,kmers)
        
    def kmer_direct(location_split: str, kmer: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_direct(h5file,kmer)
    
    def kmer_list_direct(location_split: str, kmers: list):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmers_direct(h5file,kmers)
        
    def kmer_read(location_split: str, kmer: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_read(h5file,kmer)
    
    def kmer_list_read(location_split: str, kmers: list, additional: list):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmers_read(h5file,kmers,additional)
        
    def kmer_paired(location_split: str, kmer: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_paired(h5file,kmer)
    
    def kmer_list_paired(location_split: str, kmers: list):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmers_paired(h5file,kmers)
        
    #---
    
    def base_info(location_split: str, base: str):
        with h5py.File(location_split, mode="r") as h5file:            
            baseTable = h5file.get("/split/base")
            (baseRow,id,cache) = Split._findItem(base,baseTable)
            return Split._base_result(baseRow,h5file)
        
    def base_list_info(location_split: str, bases: list):
        baseList = []
        for base in set(bases):
            baseList.append(base)
        baseList.sort()
        response = []
        with h5py.File(location_split, mode="r") as h5file:            
            baseTable = h5file.get("/split/base")
            number = baseTable.shape[0]
            start = 0
            cache = {}
            for i in range(len(baseList)):
                base = baseList[i]
                (baseRow,id,cache) = Split._findItem(base,baseTable,start,number,cache)
                if not baseRow==None:
                    response.append(Split._base_result(baseRow,h5file))
                    start = id+1
                else:
                    start = id 
        return response
        
        