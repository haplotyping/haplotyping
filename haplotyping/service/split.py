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
                "connected": {
                    "number": int(row[5][1]),
                    "connected": int(row[5][2])
                },
                "paired": int(row[6][1]),
                "cycle": int(row[7][0]),
                "reverse": int(row[8][0]),
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
        if ckmerRow and directRows:
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
    
    def _kmer_connected_result(ckmerId, ckmerRow, connectionsIndexRows: list, connectionsPairedRows: list, h5file: h5py.File):
        response = {"ckmer": ckmerRow[0].decode("ascii"), "connections":[], "paired": {}}
        ckmerTable = h5file.get("/split/ckmer")
        if len(connectionsIndexRows)>0:
            connectionsDataTable = h5file.get("/connections/data")
            for connectionsIndexRow in connectionsIndexRows:
                assert connectionsIndexRow[0]==ckmerId
                direct = connectionsIndexRow[3]
                number = connectionsIndexRow[4]
                connections = connectionsDataTable[connectionsIndexRow[1]:connectionsIndexRow[1]+connectionsIndexRow[2]]
                assert connectionsIndexRow[0] in connections
                connectionsCkmers = []
                for ckmerLink in connections:
                    if not ckmerLink==ckmerId:
                        connectionsCkmers.append(ckmerTable[ckmerLink][0].decode("ascii"))
                response["connections"].append({"connections": connectionsCkmers, 
                                                "direct": int(connectionsIndexRow[3]),
                                                "number": int(connectionsIndexRow[4])})
        for connectionsPairedRow in connectionsPairedRows:
            assert connectionsPairedRow[0]==ckmerId
            linkedCkmer = ckmerTable[connectionsPairedRow[1]][0].decode("ascii")
            response["paired"][linkedCkmer] = int(connectionsPairedRow[2])                
        return response
    
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
            (directRows,directId,directCache) = Split._findId(id,directTable)
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
                (directRows,directId,directCache) = Split._findId(id,directTable)
                response.append(Split._kmer_direct_result(ckmerRow,directRows,h5file))
                start = id+1
            else:
                start = id 
        return list(filter(None, response))
    
    def _kmer_connected(h5file: h5py.File, kmer: str):
        ckmer = haplotyping.General.canonical(kmer)
        ckmerTable = h5file.get("/split/ckmer")
        (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable)
        if ckmerRow:
            if ckmerRow[5][1]>0:   
                connectionsIndexTable = h5file.get("/connections/index")                
                connectionsIndexRows = connectionsIndexTable[ckmerRow[5][0]:ckmerRow[5][0]+ckmerRow[5][1]]               
            else:
                connectionsIndexRows = []
            if ckmerRow[6][1]>0:   
                connectionsPairedTable = h5file.get("/connections/paired")                
                connectionsPairedRows = connectionsPairedTable[ckmerRow[6][0]:ckmerRow[6][0]+ckmerRow[6][1]]
            else:
                connectionsPairedRows = []
            return Split._kmer_connected_result(id, ckmerRow, connectionsIndexRows, connectionsPairedRows, h5file)
        else:
            return None
    
    def _kmers_connected(h5file: h5py.File, kmers: list):
        #get canonical k-mers
        ckmerList = set()
        for kmer in kmers:
            ckmerList.add(haplotyping.General.canonical(kmer))
        ckmerList = list(ckmerList)
        ckmerList.sort()
        #get ids
        ckmerTable = h5file.get("/split/ckmer")
        number = ckmerTable.shape[0]
        start = 0
        cache = {}
        #construct response
        response = []
        for i in range(len(ckmerList)):
            ckmer = ckmerList[i]
            (ckmerRow,id,cache) = Split._findItem(ckmer,ckmerTable,start,number,cache)
            if ckmerRow:
                start = id+1
                if ckmerRow[5][1]>0:   
                    connectionsIndexTable = h5file.get("/connections/index")                
                    connectionsIndexRows = connectionsIndexTable[ckmerRow[5][0]:ckmerRow[5][0]+ckmerRow[5][1]]               
                else:
                    connectionsIndexRows = []
                if ckmerRow[6][1]>0:   
                    connectionsPairedTable = h5file.get("/connections/paired")                
                    connectionsPairedRows = connectionsPairedTable[ckmerRow[6][0]:ckmerRow[6][0]+ckmerRow[6][1]]
                else:
                    connectionsPairedRows = []
                response.append(Split._kmer_connected_result(id, ckmerRow, connectionsIndexRows, 
                                                             connectionsPairedRows, h5file))
            else:
                start = id 
        return response
            
   
    #---
    
    def info(location_split: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._info(h5file)
    
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
    
    def kmer_connected(location_split: str, kmer: str):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmer_connected(h5file,kmer)
    
    def kmer_list_connected(location_split: str, kmers: list):
        with h5py.File(location_split, mode="r") as h5file:            
            return Split._kmers_connected(h5file,kmers)
    
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
        
        