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
    
    def _kmer_read_result(kmerIds,readInfoList,readData,h5file,kmerDict={}):
        if len(kmerIds)>0:
            ckmerTable = h5file.get("/split/ckmer")
            response = []
            n = 0
            checkIds = set(kmerIds)
            for item in readInfoList:
                read = readData[n:n+item[0]]
                if len(checkIds.intersection(read))>0:
                    kmerList = []
                    for kmerId in read:
                        if not kmerId in kmerDict:
                            entry = ckmerTable[kmerId]
                            kmerDict[kmerId] = (entry[0].decode("ascii"),Split._translate_type(entry[1].decode("ascii")),int(entry[2]))
                        kmerList.append(kmerDict[kmerId])
                    response.append({"kmers": kmerList, "number": int(item[1])})
                n+=item[0]
        else:
            response = []
        return response,kmerDict
    
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
    
    def _kmers_connected(h5file: h5py.File, kmers: list):
        ckmerList = set()
        for kmer in kmers:
            ckmerList.add(haplotyping.General.canonical(kmer))
        ckmerList = list(ckmerList)
        ckmerList.sort()
        response = []
        partitions = set()
        kmerIds = set()
        ckmerTable = h5file.get("/split/ckmer")
        readPartitionTable = h5file.get("/relations/readPartition")
        readDataTable = h5file.get("/relations/readData")
        readInfoTable = h5file.get("/relations/readInfo")
        #get partitions
        number = ckmerTable.shape[0]
        start = 0
        cache = {}
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable,start,number,cache)
            if ckmerRow:
                partitions.add(int(ckmerRow[5]))
                kmerIds.add(id)
                start = id+1
            else:
                start = id 
        #check reads from partitions
        partitions = list(partitions)
        partitions.sort()
        for i in range(len(partitions)):
            partition = readPartitionTable[partitions[i]]
            readData = readDataTable[partition[0][0]:partition[0][0]+partition[0][1]]
            readInfo = readInfoTable[partition[1][0]:partition[1][0]+partition[1][1]]
            
            print(len(readData),"readData")
            for item in readData:
                print(item)
            print(len(readInfo),"readInfo")
            for item in readInfo:
                print(item)
            print(partition)
            print(readPartitionTable[partitions[i]+1])
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
            kmerDict = {id: (ckmerRow[0].decode("ascii"),Split._translate_type(ckmerRow[1].decode("ascii")),int(ckmerRow[2]))}
            partitionRow = readPartitionTable[ckmerRow[5]]
            readDataList = readDataTable[partitionRow[0][0]:partitionRow[0][0]+partitionRow[0][1]]
            readInfoList = readInfoTable[partitionRow[1][0]:partitionRow[1][0]+partitionRow[1][1]]
            reads,kmerDict = Split._kmer_read_result([id],readInfoList,readDataList,h5file,kmerDict)
        else:
            reads = []
        return reads
        
    def _kmers_read(h5file: h5py.File, kmers: list):
        ckmerList = set()
        for kmer in kmers:
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
        #get partition data from k-mers
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmerList[i],ckmerTable,start,number,cache)
            if ckmerRow:
                kmerDict[id] = (ckmerRow[0].decode("ascii"),Split._translate_type(ckmerRow[1].decode("ascii")),int(ckmerRow[2]))
                kmerIds.append(id)
                partitions.add(ckmerRow[5])
                start = id+1
            else:
                start = id 
        #collect reads from partitions
        reads = []
        for p in partitions:
            partitionRow = readPartitionTable[p]
            readDataList = readDataTable[partitionRow[0][0]:partitionRow[0][0]+partitionRow[0][1]]
            readInfoList = readInfoTable[partitionRow[1][0]:partitionRow[1][0]+partitionRow[1][1]]
            newReads,kmerDict = Split._kmer_read_result(kmerIds,readInfoList,readDataList,h5file,kmerDict)
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
    
    def kmer_list_read(location_split: str, kmers: list):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmers_read(h5file,kmers)
        
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
        
        