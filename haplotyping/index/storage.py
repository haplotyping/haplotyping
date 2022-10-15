import haplotyping.index.database
from multiprocessing import shared_memory, current_process
from queue import Empty
import os, re, pickle, tables, statistics, logging
import numpy as np, math
from contextlib import ExitStack

class Storage:
    
    """
    Internal use, storage and processing
    """
    
    def create_merge_storage(pytablesStorage, numberOfKmers, maximumFrequency, maximumConnectionLength,
                             nCycle=None, nReversal=None, nDirect=None, 
                             nConnections=None, nPaired=None):
        pytablesStorage.create_table(pytablesStorage.root, 
            "cycle",{
            "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "minimumLength": tables.UInt16Col(pos=1),
            "number": haplotyping.index.Database.getTablesUint(maximumFrequency,2),
        }, "Cycles", expectedrows=nCycle)
        pytablesStorage.create_table(pytablesStorage.root, 
            "reversal",{
            "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "minimumLength": tables.UInt16Col(pos=1),
            "number": haplotyping.index.Database.getTablesUint(maximumFrequency,2),
        }, "Reversals", expectedrows=nReversal)
        pytablesStorage.create_table(pytablesStorage.root, 
            "direct",{
            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "fromDirection": tables.StringCol(itemsize=1,pos=1),
            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,2),
            "toDirection": tables.StringCol(itemsize=1,pos=3),
            "number": haplotyping.index.Database.getTablesUint(maximumFrequency,4),
            "distance": tables.UInt16Col(pos=5),
            "splitDirection": tables.StringCol(itemsize=1,pos=6),
            "reverseBase": haplotyping.index.Database.getTablesUint(numberOfKmers,7),
            "forwardBase": haplotyping.index.Database.getTablesUint(numberOfKmers,8),
            "problematic": tables.UInt8Col(pos=9),
        }, "Direct relations", expectedrows=nDirect)
        pytablesStorage.create_table(pytablesStorage.root, 
            "deleteDirect",{
            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "fromDirection": tables.StringCol(itemsize=1,pos=1),
            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,2),
            "toDirection": tables.StringCol(itemsize=1,pos=3),
        }, "Incorrect direct relations", expectedrows=nDirect)
        pytablesStorage.create_table(pytablesStorage.root, 
            "connections",{
            "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "dataLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
            "length": haplotyping.index.Database.getTablesUint(maximumConnectionLength,2),
            "number": haplotyping.index.Database.getTablesUint(maximumFrequency,3),
            "direct": tables.UInt8Col(pos=4),
        }, "Indirect relations", expectedrows=nConnections)
        pytablesStorage.create_earray(pytablesStorage.root, 
            "data",haplotyping.index.Database.getTablesUintAtom(numberOfKmers),(0,), 
            "Data", expectedrows=nConnections)
        pytablesStorage.create_table(pytablesStorage.root, 
            "paired",{
            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
            "number": haplotyping.index.Database.getTablesUint(maximumFrequency,2),
        }, "Paired relations", expectedrows=nPaired)
        pytablesStorage.create_table(pytablesStorage.root, 
            "deletePaired",{
            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
        }, "Incorrect paired relations", expectedrows=nPaired)
        
    def worker_automaton(shutdown_event,queue_automaton,queue_index,queue_finished,k,automatonKmerSize,automatonFile):
        
        logger = logging.getLogger(__name__)
        
        def compute_matches(sequence,automatonSplits):
            rsequence = haplotyping.General.reverse_complement(sequence)
            boundary = (len(sequence)-k+automatonKmerSize-1)
            correction_reverse = len(sequence) + automatonKmerSize - 1 - k
            correction_forward = automatonKmerSize - 1
            #use automaton to check reverse complement sequence
            rdict = {correction_reverse - end_index: (number, startLinks,)
                     for (end_index, (number,startLinks)) in automatonSplits.iter(rsequence) 
                     if end_index <= boundary}
            if len(rdict.keys())>0:
                #use automaton to check sequence
                flist = [(end_index - correction_forward, number, startLinks)
                         for (end_index, (number,startLinks)) in automatonSplits.iter(sequence)
                         if end_index <= boundary]
                #combine results
                clist = [(pos, (number, startLinks,) ,rdict[pos]) 
                         for (pos, number, startLinks) in flist 
                         if pos in rdict.keys()]
            else:
                clist = []
            return clist
        
        #load automaton into memory
        with open(automatonFile, "rb") as f:
            automatonSplits = pickle.load(f)
        logger.debug("automaton: fsm loaded from {} MB".format(round(os.stat(automatonFile).st_size/1048576)))
            
        while not shutdown_event.is_set():
            try:
                item = queue_automaton.get(block=True, timeout=1)
                if item==None:
                    break
                elif isinstance(item,tuple) and len(item)==2:
                    queue_index.put((
                        (item[0],compute_matches(item[0],automatonSplits),),
                        (item[1],compute_matches(item[1],automatonSplits),),
                    ))
                    pass
                elif isinstance(item,str):
                    queue_index.put((
                        (item,compute_matches(item,automatonSplits),),
                    ))
            except Empty:
                continue
        del automatonSplits
        logger.debug("automaton: fsm released")
        queue_finished.put("automaton")
            
                
    
    def worker_index(shutdown_event,queue_index,queue_matches,queue_finished,k,shm_name):

        logger = logging.getLogger(__name__)
        problemPattern = re.compile(r"["+"".join(haplotyping.index.Database.letters)+
                                         "][^"+"".join(haplotyping.index.Database.letters)+"]+")
        
        def problemStartPositions(sequence):
            problemStartPositions = []        
            for m in re.finditer(problemPattern, sequence):
                problemStartPositions.append(m.span()[0]+1)
            return problemStartPositions

        def compute_matches(sequence,clist):
            history = {}
            matchesList = []
            problems = problemStartPositions(sequence)
            relevantProblem = None if len(problems)==0 else problems[0]
            matches = []
            
            previousLink = -1
            previousPos = -1
            
            positions = []
            for pos, (forward_number,forward_startLinks), (reverse_number,reverse_startLinks) in clist:
                positions.append(pos)
                #check if really match
                kmer = sequence[pos:pos+k]            
                #only possible if is reversed
                if forward_number==0:  
                    if reverse_number==0:
                        continue
                    else:
                        rkmer = haplotyping.General.reverse_complement(kmer)
                        if rkmer>kmer:
                            continue
                        else:
                            startLinks = max(forward_startLinks,reverse_startLinks)
                            endLinks = reverse_startLinks+reverse_number
                            foundMatch = False
                            for i in range(startLinks,endLinks):
                                indexLink = i*k
                                if rkmer==shm.buf[indexLink:indexLink+k].tobytes().decode():
                                    link = i
                                    orientation = "r"
                                    foundMatch = True
                                    break
                            if not foundMatch:
                                continue
                else:
                    rkmer = haplotyping.General.reverse_complement(kmer)
                    #only possible if forward
                    if reverse_number==0:
                        if rkmer<kmer:
                            continue
                        else:
                            startLinks = max(forward_startLinks,reverse_startLinks)
                            endLinks = forward_startLinks+forward_number
                            foundMatch = False
                            for i in range(startLinks,endLinks):
                                indexLink = i*k
                                if kmer==shm.buf[indexLink:indexLink+k].tobytes().decode():
                                    link = i
                                    orientation = "c"
                                    foundMatch = True
                                    break
                            if not foundMatch:
                                continue
                    else:
                        foundMatch = False
                        #if forward
                        if kmer<rkmer:
                            startLinks = forward_startLinks
                            endLinks = forward_startLinks+forward_number
                            for i in range(startLinks,endLinks):
                                indexLink = i*k
                                if kmer==shm.buf[indexLink:indexLink+k].tobytes().decode():
                                    link = i
                                    orientation = "c"
                                    foundMatch = True
                                    break
                            if not foundMatch:
                                #no match
                                continue
                        #if reverse
                        else:
                            startLinks = reverse_startLinks
                            endLinks = reverse_startLinks+reverse_number
                            foundMatch = False
                            for i in range(startLinks,endLinks):
                                indexLink = i*k
                                if rkmer==shm.buf[indexLink:indexLink+k].tobytes().decode():
                                    link = i
                                    orientation = "r"
                                    foundMatch = True
                                    break
                            if not foundMatch:
                                continue

                history[link]=[orientation,pos]
                #check if a problem did occur between last match and current
                if relevantProblem and len(matches)>0 and pos>relevantProblem:
                    matchesList.append(matches)
                    matches=[]  
                matches.append([pos,link,orientation])
                previousPos = pos
                previousLink = link
                #compute where the next problem will occur
                if relevantProblem and pos>relevantProblem:
                    relevantProblem = None
                    for problem in problems:
                        if problem>pos:
                            relevantProblem=problem
                            break                   
            if len(matches)>0:
                matchesList.append(matches)
            return (matchesList, len(matchesList)==1)
        
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("index: shared memory of {} MB used".format(round(shm.size/1048576)))
        try:
            while not shutdown_event.is_set():
                try:
                    item = queue_index.get(block=True, timeout=1)
                    if item==None:
                        break
                    elif isinstance(item,tuple):
                        if len(item)==1:
                            (matches, direct, ) = compute_matches(item[0][0],item[0][1])
                            if len(matches)==0 or (len(matches)==1 and len(matches[0])==1):
                                pass
                            else:
                                queue_matches.put(((matches, direct,),))
                        elif len(item)==2:
                            sequence0 = item[0][0]
                            sequence1 = item[1][0]
                            (matches0, direct0, ) = compute_matches(item[0][0],item[0][1])
                            (matches1, direct1, ) = compute_matches(item[1][0],item[1][1])
                            if len(matches0)==0 and len(matches1)==0:
                                pass
                            elif len(matches0)==0:
                                if len(matches1)==1 and len(matches1[0])==1:
                                    pass
                                else:
                                    queue_matches.put(((matches1, direct1,),))
                            elif len(matches1)==0:
                                if len(matches0)==1 and len(matches0[0])==1:
                                    pass
                                else:
                                    queue_matches.put(((matches0, direct0, ), ))
                            else:
                                queue_matches.put(((matches0, direct0, ),
                                                   (matches1, direct1, )))
                except Empty:
                    continue
        finally:
            shm.close()
        logger.debug("index: shared memory released")
        queue_finished.put("index")
            
    
    """
    Process matches into direct connections and queue full list of matches for indirect connections
    """
    def worker_matches(shutdown_event,queue_matches,queue_connections,queue_storage,queue_finished,
                       filenameBase,numberOfKmers,maximumFrequency,estimatedMaximumReadLength,arrayNumberDirect):
        
        logger = logging.getLogger(__name__)
        try:
            curr_proc = current_process()
            pytablesFileWorker = filenameBase+"_tmp_direct_{}.process.h5".format(curr_proc.name)
            if os.path.exists(pytablesFileWorker):
                os.remove(pytablesFileWorker)

            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                
                logger.debug("matches: store direct connections (dimension: {})".format(arrayNumberDirect))
                
                #define correct maxValues based on previous results             
                dtype = [("cycle", [("number", haplotyping.index.Database.getUint(maximumFrequency)),
                                    ("minimum", haplotyping.index.Database.getUint(
                                                      2*estimatedMaximumReadLength))]),
                         ("reversal", [("number", haplotyping.index.Database.getUint(maximumFrequency)),
                                       ("minimum", haplotyping.index.Database.getUint(
                                                      2*estimatedMaximumReadLength))]),
                         ("direct", [("d"+str(i),[("type","uint8"),
                                                  ("link",haplotyping.index.Database.getUint(numberOfKmers)),
                                                  ("distance",haplotyping.index.Database.getUint(
                                                      2*estimatedMaximumReadLength)),
                                                  ("number",haplotyping.index.Database.getUint(maximumFrequency))]) 
                                     for i in range(arrayNumberDirect)]),]
                connections = np.ndarray((numberOfKmers,), dtype=dtype, order="C")
                connections.fill(((0,0,),(0,0,),tuple((0,0,0,0,) for i in range(arrayNumberDirect))))
                
                logger.debug("created memory storage for direct connections: {} MB".format(round(connections.nbytes/1048576)))
                
                #define correct maxValues based on previous results
                tableDirectOther = pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                    "directOther",{
                    "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                    "fromDirection": tables.UInt8Col(pos=1),
                    "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,2),
                    "toDirection": tables.UInt8Col(pos=3),
                    "number": haplotyping.index.Database.getTablesUint(maximumFrequency,4),
                    "distance": haplotyping.index.Database.getTablesUint(2*estimatedMaximumReadLength,5),
                }, "Temporary to dump other direct relations", track_times=False)

                def store_cycle(fromLink, distance):
                    (previousNumber,previousDistance,) = connections[fromLink][0]
                    if (previousDistance==0) or (distance<previousDistance):
                        connections[fromLink][0] = (previousNumber+1,distance,)
                    else:
                        connections[fromLink][0] = (previousNumber+1,previousDistance,)

                def store_reversal(fromLink, distance):
                    (previousNumber,previousDistance,) = connections[fromLink][1]
                    if (previousDistance==0) or (distance<previousDistance):
                        connections[fromLink][1] = (previousNumber+1,distance,)
                    else:
                        connections[fromLink][1] = (previousNumber+1,previousDistance,)

                def store_direct(fromLink, fromDirection, toLink, toDirection, distance):
                    numberFilled = arrayNumberDirect
                    connectionType = 1 + fromDirection + (2*toDirection)
                    directRow = connections[fromLink][2]
                    for i in range(arrayNumberDirect):
                        if directRow[i][0]==0:
                            numberFilled = i
                            break
                        elif (directRow[i][1]==toLink and directRow[i][0]==connectionType 
                              and directRow[i][2]==distance):
                            directRow[i][3]+=1 
                            return
                    if numberFilled<arrayNumberDirect:
                        directRow[numberFilled] = (connectionType,toLink,distance,1,)
                    else:
                        #try to smartly swap with other direct link
                        #doing this will increase the probability that only valid connections are stored
                        #in these arrays, and mainly invalid connections end up in the additional table
                        sameFromDirectionTypes = [1 + fromDirection + (2*x) for x in range(2)]
                        sameFromDirectionDistances = []
                        sameFromDirectionEntries = []
                        otherFromDirectionDistances = []
                        otherFromDirectionEntries = []
                        swapEntry = None
                        #get distances
                        for i in range(numberFilled):
                            if directRow[i][0] in sameFromDirectionTypes:
                                sameFromDirectionDistances.append(directRow[i][2])
                                sameFromDirectionEntries.append(i)
                            elif directRow[i][0]>0:
                                otherFromDirectionDistances.append(directRow[i][2])
                                otherFromDirectionEntries.append(i)
                            else:
                                break
                        #too much in same direction
                        if 2*len(sameFromDirectionDistances)>=arrayNumberDirect:
                            sameFromDirectionDistance = sorted(statistics.multimode(sameFromDirectionDistances))[0]
                            if distance==sameFromDirectionDistance:
                                candidates=[]
                                for i,d in zip(sameFromDirectionEntries,sameFromDirectionDistances):
                                    if not d==sameFromDirectionDistance:
                                        candidates.append(i)
                                if len(candidates)>0:
                                    minimum = min([directRow[i][3] for i in candidates])
                                    for i in candidates:
                                        if directRow[i][3]==minimum:
                                            swapEntry = i
                                            break
                        #too much in other direction
                        elif 2*len(otherFromDirectionDistances)>arrayNumberDirect:
                            otherFromDirectionDistance = sorted(statistics.multimode(otherFromDirectionDistances))[0]
                            candidates=[]
                            for i,d in zip(otherFromDirectionEntries,otherFromDirectionDistances):
                                if not d==otherFromDirectionDistance:
                                    candidates.append(i)
                            if len(candidates)>0:
                                minimum = min([directRow[i][3] for i in candidates])
                                for i in candidates:
                                    if directRow[i][3]==minimum:
                                        swapEntry = i
                                        break
                        #only swap if candidate has been found
                        dumpDirectRow = tableDirectOther.row
                        dumpDirectRow["fromLink"] = fromLink
                        if not swapEntry==None:
                            dumpDirectRow["fromDirection"] = (directRow[swapEntry][0] - 1) & 1
                            dumpDirectRow["toLink"] = directRow[swapEntry][1]
                            dumpDirectRow["toDirection"] = int(((directRow[swapEntry][0] - 1) & 2)/2)
                            dumpDirectRow["number"] = directRow[swapEntry][3]
                            dumpDirectRow["distance"] = directRow[swapEntry][2]
                            directRow[swapEntry] = (connectionType,toLink,distance,1,)
                        else:
                            dumpDirectRow["fromDirection"] = fromDirection
                            dumpDirectRow["toLink"] = toLink
                            dumpDirectRow["toDirection"] = toDirection
                            dumpDirectRow["number"] = 1
                            dumpDirectRow["distance"] = distance
                        dumpDirectRow.append()

                def process_matches(matchesList):
                    history = {}
                    links = []
                    startPos = 0
                    endPos = 0
                    for matches in matchesList:
                        if len(matches)>0:
                            previousPos = matches[0][0]
                            previousLink = matches[0][1]
                            startPos = previousPos
                            #0: left, 1:right
                            previousDirection = (1 if matches[0][2]=="c" else 0)
                            links.append((previousLink,previousDirection))
                            #cycles and reversals
                            if previousLink in history.keys():
                                if history[previousLink][0]==previousDirection:
                                    store_cycle(previousLink,previousPos-history[previousLink][1])
                                else:
                                    store_reversal(previousLink,previousPos-history[previousLink][1])
                            history[previousLink]=(previousDirection,previousPos,)
                            #loop over matches
                            for i in range(1,len(matches)):
                                currentPos = matches[i][0]
                                currentLink = matches[i][1]
                                endPos = currentPos
                                #0: left, 1:right
                                currentDirection = (0 if matches[i][2]=="c" else 1)
                                links.append((currentLink,currentDirection))
                                store_direct(previousLink, previousDirection, 
                                             currentLink, currentDirection, currentPos-previousPos)
                                store_direct(currentLink, currentDirection, 
                                             previousLink, previousDirection, currentPos-previousPos)
                                previousPos = currentPos
                                previousLink = currentLink
                                previousDirection = (currentDirection+1)%2
                                #cycles and reversals
                                if currentLink in history.keys():
                                    if history[currentLink][0]==previousDirection:
                                        store_cycle(currentLink,currentPos-history[currentLink][1])
                                    else:
                                        store_reversal(currentLink,currentPos-history[currentLink][1])
                                history[previousLink]=(previousDirection,previousPos,)
                    return (links,max(0,endPos-startPos),)

                while not shutdown_event.is_set():
                    try:
                        item = queue_matches.get(block=True, timeout=1)
                        if item==None:
                            break
                        elif isinstance(item,tuple):
                            if len(item)==1:
                                matchesList = item[0][0]
                                directConnected = item[0][1]
                                (links,length,) = process_matches(matchesList)
                                if len(links)>2:
                                    queue_connections.put(((links,length,directConnected,)))
                            elif len(item)==2:
                                matchesList0 = item[0][0]
                                directConnected0 = item[0][1]
                                matchesList1 = item[1][0]
                                directConnected1 = item[1][1]
                                (links0,length0,) = process_matches(matchesList0)
                                (links1,length1,) = process_matches(matchesList1)
                                if len(links0)==0:
                                    if len(links1)>2:
                                        queue_connections.put(((links1,length1,directConnected1,)))
                                elif len(links1)==0:
                                    if len(links0)>2:
                                        queue_connections.put(((links0,length0,directConnected0,)))
                                else:
                                    queue_connections.put((
                                        (links0,length0,directConnected0,),
                                        (links1,length1,directConnected1,),
                                    ))
                    except Empty:
                        continue   
                pytablesStorageWorker.flush()
                tableDirectOther.cols.fromLink.create_csindex()
                pytablesStorageWorker.flush()
                pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                                              name="direct", obj=connections, expectedrows=numberOfKmers)
                logger.debug("matches: memory storage saved, other connections indexed")
                pytablesStorageWorker.flush()
                
            #now the file can be released for merge
            queue_storage.put(pytablesFileWorker) 
        finally:
            pass
        queue_finished.put("matches")
            
    """
    Process indirect connections: filtering and efficient storage
    """
    def worker_connections(shutdown_event,queue_connections,queue_storage,queue_finished,
                       filenameBase,numberOfKmers,arrayNumberConnection,maximumFrequency,shm_name):
        
        logger = logging.getLogger(__name__)
        
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("connections: shared memory of {} MB used".format(round(shm.size/1048576)))
            
        try:
            curr_proc = current_process()
            pytablesFileWorker = filenameBase+"_tmp_connections_"+str(curr_proc.name)+".process.h5"
            if os.path.exists(pytablesFileWorker):
                os.remove(pytablesFileWorker)

            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                
                logger.debug("connections: store additional connections (dimension: {})".format(arrayNumberConnection))

                #connections - define correct maxValues based on previous results             
                dtype = [("number", haplotyping.index.Database.getUint(arrayNumberConnection+1)),
                         ("connections", [("c"+str(i),[
                                                        ("link",haplotyping.index.Database.getUint(numberOfKmers)),
                                                        ("length",haplotyping.index.Database.getUint(arrayNumberConnection)),
                                                        ("direct","uint8"),
                                                        ("hash","int64"),
                                                        ("number",haplotyping.index.Database.getUint(maximumFrequency))
                                                      ]) 
                                                      for i in range(arrayNumberConnection)])]
                connections = np.ndarray((numberOfKmers,), dtype=dtype, order="C")
                connections.fill((0,tuple((0,0,0,0,0,) for i in range(arrayNumberConnection))))
                
                #paired entries - define correct maxValues based on previous results   
                #can't however really estimate the required size because this depends 
                #on the (unknown) distance between paired reads
                dtype = [("number", haplotyping.index.Database.getUint(arrayNumberConnection+1)),
                         ("paired", [("p"+str(i),[("link",haplotyping.index.Database.getUint(numberOfKmers)),
                                                  ("number",haplotyping.index.Database.getUint(maximumFrequency))]) 
                                     for i in range(arrayNumberConnection)]),]
                paired = np.ndarray((numberOfKmers,), dtype=dtype, order="C")
                paired.fill((0,tuple((0,0,) for i in range(arrayNumberConnection))))
                
                logger.debug("created memory storage for additional connections: {} MB".format(
                    round(connections.nbytes/1048576)))
                logger.debug("created memory storage for paired entries: {} MB".format(
                    round(paired.nbytes/1048576)))
                
                #define correct maxValues based on previous results
                tablePairedOther = pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                    "pairedOther",{
                    "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                    "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1)
                }, "Temporary to dump other paired relations", track_times=False)
                
                mapSplitDirection = [{b"l":"l", b"r":"r", b"b":"b"},{b"l":"r", b"r":"l", b"b":"b"}]
                
                def store_connections(connectionsList,direct,tableConnectionsData,tableConnectionsDataNumber):
                    #compact pairs with equal dosage
                    compactedConnectionsList = []
                    for id, (link,direction) in enumerate(connectionsList):
                        props = kmer_properties[link] 
                        newEntry = (link,mapSplitDirection[direction][props[0]],props[1],id,)
                        #contract pairs with equal dosage (rightsplitter -> leftsplitter)                        
                        if (newEntry[1]=="l" and (len(compactedConnectionsList)>0) and 
                            compactedConnectionsList[-1][1]=="r"):
                            compactedConnectionsList[-1] = (compactedConnectionsList[-1],newEntry,)
                        else:
                            compactedConnectionsList.append(newEntry)
                    #select entries with (potential) minimal dosage
                    potentialDirection = None
                    potentialId = None
                    compactedSelection = [] 
                    # entry: [ link, splitType, frequency, #id connectionsList]
                    for id,entry in enumerate(compactedConnectionsList):
                        newDirection = entry[0][1] if len(entry)==2 else entry[1]
                        if potentialId==None:
                            potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                            potentialId = id
                        else:
                            if potentialDirection=="l":
                                if newDirection=="l":
                                    pass
                                else:
                                    compactedSelection.append(potentialId)
                                    potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                                    potentialId = id
                            elif potentialDirection=="r":
                                if newDirection=="l":
                                    raise Exception("can't happen")
                                else:
                                    potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                                    potentialId = id
                            else:
                                if newDirection=="l":
                                    pass
                                else:
                                    compactedSelection.append(potentialId)
                                    potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                                    potentialId = id
                    if not potentialId==None:
                        compactedSelection.append(potentialId)
                    #compute final selection and linking entry
                    finalSelection = []
                    linkingId = None
                    linkingNumber = None
                    def recompute_linking_entry(linkingId,linkingNumber,newId,newNumber):
                        if linkingId==None or linkingNumber>newNumber:
                            linkingId = newId
                            linkingNumber = newNumber
                        elif linkingNumber==newNumber:
                            if connectionsList[linkingId][0]>connectionsList[newId][0]:
                                linkingId = newId
                                linkingNumber = newNumber
                        return (linkingId,linkingNumber,)
                    
                    for id in compactedSelection:
                        if len(compactedConnectionsList[id])==2:
                            if id==0:
                                finalSelection.append(compactedConnectionsList[id][1][-1])
                                (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][1][-1],
                                                                compactedConnectionsList[id][1][2])
                            elif id==len(compactedConnectionsList)-1:
                                finalSelection.append(compactedConnectionsList[id][0][-1])
                                (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][0][-1],
                                                                compactedConnectionsList[id][0][2])
                            else:
                                if compactedConnectionsList[id][0][2]>compactedConnectionsList[id][1][2]:
                                    finalSelection.append(compactedConnectionsList[id][0][-1])
                                    (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][0][-1],
                                                                compactedConnectionsList[id][0][2])
                                elif compactedConnectionsList[id][1][2]>compactedConnectionsList[id][0][2]:
                                    finalSelection.append(compactedConnectionsList[id][1][-1])
                                    (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][1][-1],
                                                                compactedConnectionsList[id][1][2])
                                elif compactedConnectionsList[id][0][0]>compactedConnectionsList[id][1][0]:
                                    finalSelection.append(compactedConnectionsList[id][0][-1])
                                    (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][0][-1],
                                                                compactedConnectionsList[id][0][2])
                                else:
                                    finalSelection.append(compactedConnectionsList[id][1][-1])
                                    (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][1][-1],
                                                                compactedConnectionsList[id][1][2])
                        else:
                            finalSelection.append(compactedConnectionsList[id][-1])
                            (linkingId,linkingNumber) = recompute_linking_entry(linkingId,linkingNumber,
                                                                compactedConnectionsList[id][-1],
                                                                compactedConnectionsList[id][2])

                    linkingConnectedKmer = None
                    
                    #single entry not informative
                    if len(compactedSelection)<=1:
                        pass
                    #directly neighhbouring entries not informative
                    elif len(compactedSelection)==2 and (compactedSelection[1]-compactedSelection[0])==1:
                        pass
                    elif len(finalSelection)>0:
                        #get final list of connected k-mers
                        connectedKmers = [connectionsList[id][0] for id in finalSelection]
                        linkingConnectedKmer = connectionsList[linkingId][0]
                        assert linkingConnectedKmer in connectedKmers
                        #remove neighbouring duplicates
                        connectedKmers = sorted(set(connectedKmers), key=connectedKmers.index)
                        #always sorted in same direction
                        if connectedKmers[0]>connectedKmers[-1]:
                            connectedKmers = tuple(connectedKmers[::-1])
                        else:
                            connectedKmers = tuple(connectedKmers)
                        connectionsRow = connections[linkingConnectedKmer]
                        #don't store if too much k-mers are involved or if non-informative
                        if connectionsRow[0]>arrayNumberConnection or len(connectedKmers)<=1:
                            pass
                        else:
                            lengthValue = len(connectedKmers)
                            directValue = 1 if direct else 0
                            hashValue = hash(connectedKmers)
                            stored = False
                            #check hash-based if already stored
                            for i in range(connectionsRow[0]):
                                if connectionsRow[1][i][3] == hashValue:                                    
                                    if connectionsRow[1][i][1] == lengthValue:                                        
                                        connectionsRow[1][i][2] == max(directValue,connectionsRow[1][i][2])
                                        connectionsRow[1][i][4] += 1    
                                        stored = True
                                    else:
                                        #conflict
                                        connectionsRow[0] = arrayNumberConnection+1
                                        for j in range(connectionsRow[0]):
                                            connectionsRow[1][j] = (0,0,0,0,0,)
                                    break
                            if not stored and connectionsRow[0]<=arrayNumberConnection:
                                if connectionsRow[0]<arrayNumberConnection:
                                    linkValue = tableConnectionsDataNumber
                                    tableConnectionsData.append(connectedKmers)
                                    tableConnectionsDataNumber+=len(connectedKmers)
                                    connectionsRow[1][connectionsRow[0]] = (linkValue,lengthValue,
                                                                            directValue,hashValue,1,)
                                    connectionsRow[0]+=1
                                    stored = True
                                else:
                                    #no space to store, reset all
                                    for j in range(connectionsRow[0]):
                                        connectionsRow[1][j] = (0,0,0,0,0,)
                                    connectionsRow[0] = arrayNumberConnection+1
                                    
                            #store (in memory)
                            connections[linkingConnectedKmer] = connectionsRow
                            
                    return (tableConnectionsDataNumber,linkingConnectedKmer)
                
                def store_paired(linkedCkmer0,linkedCkmer1):
                    #only store if potentially necessary
                    if linkedCkmer0==linkedCkmer1:
                        pass
                    elif connections[linkedCkmer0][0]>arrayNumberConnection:
                        pass
                    elif connections[linkedCkmer1][0]>arrayNumberConnection:
                        pass
                    else:
                        pairedRow = paired[linkedCkmer0]
                        stored = False
                        for i in range(pairedRow[0]):
                            if pairedRow[1][i][0] == linkedCkmer1:                                    
                                pairedRow[1][i][1] += 1    
                                stored = True                                
                        if not stored:
                            if pairedRow[0]<arrayNumberConnection:
                                pairedRow[1][pairedRow[0]] = (linkedCkmer1,1,)
                                pairedRow[0]+=1                                
                            else:
                                dumpPairedRow = tablePairedOther.row
                                dumpPairedRow["fromLink"] = linkedCkmer0
                                dumpPairedRow["toLink"] = linkedCkmer1
                                dumpPairedRow.append()
                        #store (in memory)
                        paired[linkedCkmer0] = pairedRow
                    
                #get shared memory
                shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
                shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maximumFrequency)).type
                kmer_properties = np.ndarray((numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                                   ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm.buf)
            
                tableConnectionsData = pytablesStorageWorker.create_earray(pytablesStorageWorker.root, "data", 
                               haplotyping.index.Database.getTablesUintAtom(numberOfKmers), (0,), "Data")
                tableConnectionsDataNumber = 0

                while not shutdown_event.is_set():
                    try:
                        item = queue_connections.get(block=True, timeout=1)
                        if item==None:
                            break
                        elif isinstance(item,tuple):
                            if len(item)==1:
                                if len(item[0][0])>2:
                                    tableConnectionsDataNumber = store_connections(
                                                      item[0][0],item[0][2],
                                                      tableConnectionsData,tableConnectionsDataNumber)[0]
                            elif len(item)==2:
                                if len(item[0][0])==0:
                                    if len(item[1][0])>0:
                                        tableConnectionsDataNumber = store_connections(
                                                      item[1][0],item[1][2],
                                                      tableConnectionsData,tableConnectionsDataNumber)[0]
                                elif len(item[1][0])==0:
                                    if len(item[0][0])>0:
                                        tableConnectionsDataNumber = store_connections(
                                                      item[0][0],item[0][2],
                                                      tableConnectionsData,tableConnectionsDataNumber)[0]
                                else:
                                    (tableConnectionsDataNumber,linkedCkmer0) = store_connections(
                                                      item[0][0],item[0][2],
                                                      tableConnectionsData,tableConnectionsDataNumber)
                                    (tableConnectionsDataNumber,linkedCkmer1) = store_connections(
                                                      item[1][0],item[1][2],
                                                      tableConnectionsData,tableConnectionsDataNumber)
                                    #store paired connections
                                    if not (linkedCkmer0==None or linkedCkmer1==None):
                                        store_paired(linkedCkmer0,linkedCkmer1)
                                        store_paired(linkedCkmer1,linkedCkmer0)
                    except Empty:
                        continue   
                
                pytablesStorageWorker.flush()                
                tablePairedOther.cols.fromLink.create_csindex()
                pytablesStorageWorker.flush()                
                pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                    name="connections", obj=connections, expectedrows=numberOfKmers)
                logger.debug("matches: memory storage connections saved")
                pytablesStorageWorker.flush()
                pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                    name="paired", obj=paired, expectedrows=numberOfKmers)
                logger.debug("matches: memory storage paired saved")
                pytablesStorageWorker.flush()
                
            #now the file can be released for merge
            queue_storage.put(pytablesFileWorker) 
        finally:
            shm.close()
            logger.debug("index: shared memory released")

        queue_finished.put("connections")
                
            
    """
    Merge stored direct and indirect connections
    """                
    def worker_merges(shutdown_event,queue_ranges,queue_merges,storageDirectFiles,storageConnectionFiles,
                      filenameBase,numberOfKmers,maximumFrequency,minimumFrequency,maximumConnectionLength,shm_name):
        
        logger = logging.getLogger(__name__)

        def merge_direct_storage(pytablesFileWorker, storageDirectFiles, mergeStart, mergeNumber, 
                                 numberOfKmers, kmer_properties):

            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                tableCycle = pytablesStorageWorker.root.cycle
                tableReversal = pytablesStorageWorker.root.reversal
                tableDirect = pytablesStorageWorker.root.direct
                tableDeleteDirect = pytablesStorageWorker.root.deleteDirect                

                with ExitStack() as stack:
                    #get handlers
                    storageHandlers = [{"handler": stack.enter_context(tables.open_file(fname, mode="r"))} 
                                               for fname in storageDirectFiles]
                    
                    #get and position sorted iterators for other entries
                    for i in range(len(storageHandlers)):
                        storageHandlers[i]["otherRow"] = None
                        if storageHandlers[i]["handler"].root.directOther.shape[0]>0:
                            storageHandlers[i]["otherIterator"] = storageHandlers[i]["handler"].root.directOther.itersorted(
                                "fromLink",checkCSI=True)
                            for row in storageHandlers[i]["otherIterator"]:
                                if row["fromLink"]>=mergeStart:
                                    if row["fromLink"]<=mergeEnd:
                                        storageHandlers[i]["otherRow"] = row
                                    break
                        else:
                            storageHandlers[i]["otherIterator"] = None
                    
                    #direct connections
                    def add_directData(directDataEntry,fromDirection,toLink,toDirection,distance,number):
                        if not fromDirection in directDataEntry.keys():
                            directDataEntry[fromDirection] = {}
                        if not distance in directDataEntry[fromDirection].keys():
                            directDataEntry[fromDirection][distance] = {}
                        if not toLink in directDataEntry[fromDirection][distance].keys():
                            directDataEntry[fromDirection][distance][toLink] = {}
                        if not toDirection in directDataEntry[fromDirection][distance][toLink].keys():
                            props = kmer_properties[toLink]                                            
                            splitDirection = props[0].tobytes().decode()
                            reverseBase = None
                            forwardBase = None
                            if toDirection==0:
                                #from the left
                                if splitDirection in ["l","b"]:
                                    reverseBase = props[2]
                                    direction = "r" if splitDirection=="l" else "b"
                                if splitDirection in ["r","b"]:
                                    forwardBase = props[3]
                                    direction = "f" if splitDirection=="r" else "b"
                            elif toDirection==1:
                                #from the right
                                if splitDirection in ["r","b"]:
                                    reverseBase = props[3]
                                    direction = "r" if splitDirection=="r" else "b"
                                if splitDirection in ["l","b"]:
                                    forwardBase = props[2]
                                    direction = "f" if splitDirection=="l" else "b"
                            directDataEntry[fromDirection][distance][toLink][toDirection] = {"number": np.uint(number), 
                                         "direction": direction, "forwardBase": forwardBase, 
                                         "reverseBase": reverseBase, "problematic": 0}
                        else:
                            directDataEntry[fromDirection][distance][toLink][toDirection]["number"] = (
                                min(maximumFrequency,
                                    directDataEntry[fromDirection][distance][toLink][toDirection]["number"]+number))
                        return directDataEntry
                    
                    #delete probably wrong direct connections
                    def reduceDistances(fromLink,fromDirection,distanceData):
                        
                        def delete_connection(fromLink,fromDirection,toLink,toDirection):
                            tableDeleteDirectRow = tableDeleteDirect.row
                            tableDeleteDirectRow["fromLink"] = fromLink
                            tableDeleteDirectRow["fromDirection"] = fromDirection
                            tableDeleteDirectRow["toLink"] = toLink
                            tableDeleteDirectRow["toDirection"] = toDirection
                            tableDeleteDirectRow.append()
                        
                        #compute incorrect distances
                        incorrectDistances = []
                        for distance in distanceData.keys():
                            distanceNumber = 0
                            dataEntry = distanceData[distance] 
                            for toLink in dataEntry.keys():
                                for toDirection in dataEntry[toLink].keys():
                                    number = dataEntry[toLink][toDirection]["number"]
                                    if number>=minimumFrequency:
                                        distanceNumber+=number
                                    
                            if distanceNumber==0:
                                incorrectDistances.append(distance)
                        #remove incorrect distances
                        if len(incorrectDistances)<len(distanceData.keys()):
                            newValues = {}
                            for distance in distanceData.keys():
                                if distance in incorrectDistances:
                                    dataEntry = distanceData[distance] 
                                    for toLink in dataEntry.keys():
                                        for toDirection in dataEntry[toLink].keys():
                                            delete_connection(toLink,"l" if toDirection==0 else "r",
                                                  fromLink,"l" if fromDirection==0 else "r")                                   
                                else:
                                    newValues[distance] = distanceData[distance]
                            distanceData = newValues
                            
                        #compute possible splitting bases
                        reverseBases = {}
                        forwardBases = {}
                        relevantConnectionNumber = 0
                        distances = {}
                        for distance in distanceData.keys():
                            dataEntry = distanceData[distance] 
                            for toLink in dataEntry.keys():
                                for toDirection in dataEntry[toLink].keys():
                                    number = dataEntry[toLink][toDirection]["number"]
                                    if number>=minimumFrequency:
                                        distances[distance] = distances.get(distance,0) + number
                                        relevantConnectionNumber+=number
                                        direction = dataEntry[toLink][toDirection]["direction"]
                                        if direction in ["f","b"]:
                                            forwardBase = dataEntry[toLink][toDirection]["forwardBase"]
                                            forwardBases[forwardBase] = forwardBases.get(forwardBase,0) + number
                                        if direction in ["r","b"]:
                                            reverseBase = dataEntry[toLink][toDirection]["reverseBase"]
                                            reverseBases[reverseBase] = reverseBases.get(reverseBase,0) + number
                        possibleReverseBases = [base for base in reverseBases.keys() 
                                                if reverseBases[base]==relevantConnectionNumber]
                        possibleForwardBases = [base for base in forwardBases.keys() 
                                                if forwardBases[base]==relevantConnectionNumber]
                        
                        def set_problematic(distanceData,value):
                            for distance in distanceData.keys():
                                for toLink in distanceData[distance].keys():
                                    for toDirection in distanceData[distance][toLink].keys():
                                        distanceData[distance][toLink][toDirection]["problematic"] = value
                            return distanceData
                                                
                        if relevantConnectionNumber>0:
                            
                            #try fixing multiple distances if this is potential problematic
                            if len(possibleReverseBases)+len(possibleForwardBases)==0:
                                if len(distances)>1:
                                    maximumDistanceFrequency = max(max(distances.values()),minimumFrequency)
                                    boundary = math.floor(math.sqrt(maximumDistanceFrequency/minimumFrequency))
                                    incorrectDistances = [d for d in distances.keys() if distances[d]<=boundary]
                                    if len(incorrectDistances)>0:
                                        #remove incorrect distances
                                        newValues = {}
                                        for distance in distanceData.keys():
                                            if distance in incorrectDistances:
                                                dataEntry = distanceData[distance] 
                                                for toLink in dataEntry.keys():
                                                    for toDirection in dataEntry[toLink].keys():
                                                        delete_connection(toLink,"l" if toDirection==0 else "r",
                                                              fromLink,"l" if fromDirection==0 else "r")
                                            else:
                                                newValues[distance] = distanceData[distance]
                                        distanceData = newValues
                                        #recompute possible splitting bases
                                        reverseBases = {}
                                        forwardBases = {}
                                        relevantConnectionNumber = 0
                                        distances = {}
                                        for distance in distanceData.keys():
                                            dataEntry = distanceData[distance] 
                                            for toLink in dataEntry.keys():
                                                for toDirection in dataEntry[toLink].keys():
                                                    number = dataEntry[toLink][toDirection]["number"]
                                                    if number>=minimumFrequency:
                                                        distances[distance] = distances.get(distance,0) + number
                                                        relevantConnectionNumber+=number
                                                        direction = dataEntry[toLink][toDirection]["direction"]
                                                        if direction in ["f","b"]:
                                                            forwardBase = dataEntry[toLink][toDirection]["forwardBase"]
                                                            forwardBases[forwardBase] = (forwardBases.get(forwardBase,0) 
                                                                                         + number)
                                                        if direction in ["r","b"]:
                                                            reverseBase = dataEntry[toLink][toDirection]["reverseBase"]
                                                            reverseBases[reverseBase] = (reverseBases.get(reverseBase,0) 
                                                                                         + number)
                                        possibleReverseBases = [base for base in reverseBases.keys() 
                                                                if reverseBases[base]==relevantConnectionNumber]
                                        possibleForwardBases = [base for base in forwardBases.keys() 
                                                                if forwardBases[base]==relevantConnectionNumber]
                                                                                
                            
                            #no single base covering everything
                            if len(possibleReverseBases)+len(possibleForwardBases)==0:
                                distanceData = set_problematic(distanceData,2)                                
                            #found forward solution    
                            elif len(possibleForwardBases)>0:
                                forwardBase = possibleForwardBases[0]  
                                #try to reduce distances
                                if len(distances)>1:
                                    maximumDistanceFrequency = max(max(distances.values()),maximumFrequency)
                                    boundary = math.floor(math.sqrt(maximumDistanceFrequency/minimumFrequency))
                                    incorrectDistances = [d for d in distances.keys() if distances[d]<=boundary]
                                else:
                                    incorrectDistances = []
                                #reduce (if necessary) to single distance and forwardBase
                                if (len(distances) - len(incorrectDistances)) == 1:
                                    newValues = {}
                                    for distance in distanceData.keys():
                                        if distance in incorrectDistances:
                                            dataEntry = distanceData[distance] 
                                            for toLink in dataEntry.keys():
                                                for toDirection in dataEntry[toLink].keys():
                                                    delete_connection(toLink,"l" if toDirection==0 else "r",
                                                          fromLink,"l" if fromDirection==0 else "r")                           
                                        else:
                                            distanceDataEntry = {}
                                            dataEntry = distanceData[distance] 
                                            for toLink in dataEntry.keys():
                                                for toDirection in dataEntry[toLink].keys():
                                                    direction = dataEntry[toLink][toDirection]["direction"]
                                                    if ((direction in ["f","b"]) and 
                                                        (dataEntry[toLink][toDirection]["forwardBase"] == forwardBase)):
                                                        distanceDataEntry[toLink] = distanceDataEntry.get(toLink,{})
                                                        distanceDataEntry[toLink][toDirection] = dataEntry[toLink].get(
                                                            toDirection)
                                                    else:
                                                        delete_connection(toLink,"l" if toDirection==0 else "r",
                                                          fromLink,"l" if fromDirection==0 else "r")  
                                            assert len(distanceDataEntry)>0
                                            newValues[distance] = distanceDataEntry
                                    distanceData = newValues
                                else:
                                    #can't solve this, too many distances
                                    distanceData = set_problematic(distanceData,3) 
                            #found reverse solution    
                            elif len(possibleReverseBases)>0:
                                reverseBase = possibleReverseBases[0]  
                                #try to reduce distances
                                if len(distances)>1:
                                    maximumDistanceFrequency = max(distances.values())
                                    boundary = math.floor(math.sqrt(maximumDistanceFrequency/minimumFrequency))
                                    incorrectDistances = [d for d in distances.keys() if distances[d]<=boundary]
                                else:
                                    incorrectDistances = []
                                #reduce (if necessary) to single distance and reverseBase
                                if (len(distances) - len(incorrectDistances)) == 1:
                                    newValues = {}
                                    numberOfConnections = 0
                                    for distance in distanceData.keys():
                                        if distance in incorrectDistances:
                                            dataEntry = distanceData[distance] 
                                            for toLink in dataEntry.keys():
                                                for toDirection in dataEntry[toLink].keys():
                                                    delete_connection(toLink,"l" if toDirection==0 else "r",
                                                          fromLink,"l" if fromDirection==0 else "r")                           
                                        else:
                                            distanceDataEntry = {}
                                            dataEntry = distanceData[distance] 
                                            for toLink in dataEntry.keys():
                                                for toDirection in dataEntry[toLink].keys():
                                                    direction = dataEntry[toLink][toDirection]["direction"]
                                                    number = dataEntry[toLink][toDirection]["number"]
                                                    #expect single connection, therefore also restriction number
                                                    if ((direction in ["r","b"]) and 
                                                        (dataEntry[toLink][toDirection]["reverseBase"] == reverseBase) and
                                                        number>=minimumFrequency):
                                                        distanceDataEntry[toLink] = distanceDataEntry.get(toLink,{})
                                                        distanceDataEntry[toLink][toDirection] = dataEntry[toLink].get(
                                                            toDirection)
                                                        numberOfConnections+=1
                                                    else:
                                                        delete_connection(toLink,"l" if toDirection==0 else "r",
                                                          fromLink,"l" if fromDirection==0 else "r")
                                            assert len(distanceDataEntry)>0
                                            newValues[distance] = distanceDataEntry
                                    distanceData = newValues  
                                    #expect only single connection
                                    if numberOfConnections>1:
                                        distanceData = set_problematic(distanceData,5)
                                else:
                                    #can't solve this, too many distances
                                    distanceData = set_problematic(distanceData,4)
                            else:
                                raise Exception("this should not happen")
                        else:
                            #not enough relevant connections
                            distanceData = set_problematic(distanceData,1)
                        #return
                        return distanceData
                    
                    stepSizeStorage = 10000
                    for i in range(mergeStart,mergeEnd+1,stepSizeStorage):
                        #initialise
                        cycleEntries = {}
                        reversalEntries = {}
                        directData = [{} for j in range(min(mergeEnd+1-i,stepSizeStorage,(numberOfKmers-i)))]
                        for storageHandler in storageHandlers:
                            rowData = storageHandler["handler"].root.direct[i:min(mergeEnd+1,i+stepSizeStorage)]
                            for index, item in enumerate(rowData):
                                fromLink = i + index
                                #get and process regular data
                                for k in range(len(item[2])):
                                    if item[2][k][0]==0:
                                        break
                                    else:
                                        fromDirection = (item[2][k][0] - 1) & 1
                                        toDirection = int(((item[2][k][0] - 1) & 2)/2)
                                        toLink = item[2][k][1]
                                        number = int(item[2][k][3])
                                        distance = item[2][k][2]
                                        directData[index] = add_directData(directData[index],
                                                                       fromDirection,toLink,toDirection,distance,number)
                                        
                                #now check other data
                                if storageHandler["otherRow"]:
                                    if storageHandler["otherRow"]["fromLink"]<fromLink:
                                        storageHandler["otherRow"] = None
                                        for row in storageHandler["otherIterator"]:
                                            if row["fromLink"]>=fromLink:
                                                if row["fromLink"]<=mergeEnd:
                                                    storageHandler["otherRow"] = row
                                                break
                                    if storageHandler["otherRow"] and storageHandler["otherRow"]["fromLink"]==fromLink:
                                        fromDirection = storageHandler["otherRow"]["fromDirection"]
                                        toDirection = storageHandler["otherRow"]["toDirection"]
                                        toLink = storageHandler["otherRow"]["toLink"]
                                        number = int(storageHandler["otherRow"]["number"])
                                        distance = storageHandler["otherRow"]["distance"]
                                        directData[index] = add_directData(directData[index],
                                                                       fromDirection,toLink,toDirection,distance,number)
                                        for row in storageHandler["otherIterator"]:
                                            if row["fromLink"]==fromLink:
                                                storageHandler["otherRow"] = row
                                                fromDirection = storageHandler["otherRow"]["fromDirection"]
                                                toDirection = storageHandler["otherRow"]["toDirection"]
                                                toLink = storageHandler["otherRow"]["toLink"]
                                                number = int(storageHandler["otherRow"]["number"])
                                                distance = storageHandler["otherRow"]["distance"]
                                                directData[index] = add_directData(directData[index],
                                                               fromDirection,toLink,toDirection,distance,number)
                                                
                                            elif row["fromLink"]<=mergeEnd:
                                                storageHandler["otherRow"] = row
                                                break
                                            else:
                                                storageHandler["otherRow"] = None
                                                break
                                
                            for index, item in enumerate(rowData):
                                if item[0][1]>0:
                                    if index in cycleEntries.keys():
                                        cycleEntries[index] = (cycleEntries[index][0]+int(item[0][0]),
                                                               min(cycleEntries[index][1],item[0][1]),)
                                    else:
                                        cycleEntries[index] = (int(item[0][0]),item[0][1],)
                            for index, item in enumerate(rowData):
                                if item[1][1]>0:
                                    if index in reversalEntries.keys():
                                        reversalEntries[index] = (reversalEntries[index][0]+int(item[1][0]),
                                                               min(reversalEntries[index][1],item[1][1]),)
                                    else:
                                        reversalEntries[index] = (int(item[1][0]),item[1][1],)
                        for ce in cycleEntries.keys():
                            tableCycleRow = tableCycle.row
                            tableCycleRow["ckmerLink"] = i+ce
                            tableCycleRow["number"] = min(maximumFrequency,cycleEntries[ce][0])
                            tableCycleRow["minimumLength"] = cycleEntries[ce][1]
                            tableCycleRow.append()
                        for ce in reversalEntries.keys():
                            tableReversalRow = tableReversal.row
                            tableReversalRow["ckmerLink"] = i+ce
                            tableReversalRow["number"] = min(maximumFrequency,reversalEntries[ce][0])
                            tableReversalRow["minimumLength"] = reversalEntries[ce][1]
                            tableReversalRow.append()
                        
                        for j in range(len(directData)):
                            fromLink = i + j
                            for fromDirection in directData[j].keys():
                                directData[j][fromDirection]= reduceDistances(fromLink,fromDirection, 
                                                                              directData[j][fromDirection])
                                #check for multiple distances (inconsistent)
                                multipleDistances = 1 if not len(directData[j][fromDirection])==1 else 0
                                for distance in directData[j][fromDirection].keys():
                                    dataEntry = directData[j][fromDirection][distance]
                                    for toLink in dataEntry.keys():
                                        for toDirection in dataEntry[toLink].keys():
                                            number = dataEntry[toLink][toDirection]["number"]
                                            problematic = dataEntry[toLink][toDirection]["problematic"]
                                            direction = dataEntry[toLink][toDirection]["direction"]
                                            if direction in ["r","b"]:
                                                reverseBase = dataEntry[toLink][toDirection]["reverseBase"]
                                            else:
                                                reverseBase = 0
                                            if direction in ["f","b"]:
                                                forwardBase = dataEntry[toLink][toDirection]["forwardBase"]
                                            else:
                                                forwardBase = 0
                                            tableDirectRow = tableDirect.row
                                            tableDirectRow["fromLink"] = fromLink
                                            tableDirectRow["fromDirection"] = "l" if fromDirection==0 else "r"
                                            tableDirectRow["toLink"] = toLink
                                            tableDirectRow["toDirection"] = "l" if toDirection==0 else "r"
                                            tableDirectRow["number"] = min(maximumFrequency,number)
                                            tableDirectRow["distance"] = distance
                                            tableDirectRow["splitDirection"] = direction
                                            tableDirectRow["reverseBase"] = reverseBase
                                            tableDirectRow["forwardBase"] = forwardBase
                                            tableDirectRow["problematic"] = multipleDistances + (problematic<<1)
                                            tableDirectRow.append()

                #finished
                pytablesStorageWorker.flush()
                
        def merge_connection_storage(pytablesFileWorker, storageConnectionFiles, mergeStart, mergeNumber, 
                                 numberOfKmers, maximumConnectionLength):

            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                tableConnections = pytablesStorageWorker.root.connections
                tableData = pytablesStorageWorker.root.data
                tablePaired = pytablesStorageWorker.root.paired
                tableDeletePaired = pytablesStorageWorker.root.deletePaired
                tableDataCounter = 0
                
                with ExitStack() as stack:
                    #get handlers
                    storageHandlers = [{"handler": stack.enter_context(tables.open_file(fname, mode="r"))} 
                                               for fname in storageConnectionFiles]
                    
                    #get and position sorted iterators for other entries
                    for i in range(len(storageHandlers)):
                        storageHandlers[i]["otherRow"] = None
                        if storageHandlers[i]["handler"].root.pairedOther.shape[0]>0:
                            storageHandlers[i]["otherIterator"] = storageHandlers[i]["handler"].root.pairedOther.itersorted(
                                "fromLink",checkCSI=True)
                            for row in storageHandlers[i]["otherIterator"]:
                                if row["fromLink"]>=mergeStart:
                                    if row["fromLink"]<=mergeEnd:
                                        storageHandlers[i]["otherRow"] = row
                                    break
                        else:
                            storageHandlers[i]["otherIterator"] = None
                    
                    #add connections
                    def add_connectionData(connectionDataEntry,length,direct,hashValue,number):
                        if not hashValue in connectionDataEntry[1].keys():
                            connectionDataEntry[1][hashValue] = {"status": True, "length": length, 
                                                                 "direct": direct, "number": number}
                            connectionDataEntry[0]+=1
                            storeData = True
                        elif not connectionDataEntry[1][hashValue]["length"]==length:
                            connectionDataEntry[1][hashValue]["status"] = False
                            storeData = False
                        else:
                            connectionDataEntry[1][hashValue]["number"] += number  
                            connectionDataEntry[1][hashValue]["direct"] = max(direct,
                                                      connectionDataEntry[1][hashValue]["direct"])
                            storeData = False
                        return (connectionDataEntry,storeData,)
                    
                    stepSizeStorage = 10000
                    for i in range(mergeStart,mergeEnd+1,stepSizeStorage):
                        #initialise
                        connectionData = [[0,{}] for j in range(min(mergeEnd+1-i,stepSizeStorage,(numberOfKmers-i)))]
                        for storageHandler in storageHandlers:
                            rowData = storageHandler["handler"].root.connections[i:min(mergeEnd+1,i+stepSizeStorage)]
                            for index, item in enumerate(rowData):
                                fromLink = i + index
                                if item[0]>maximumConnectionLength:
                                    break
                                else:
                                    for k in range(min(item[0],len(item[1]))):
                                        lengthValue = item[1][k][1]
                                        directValue = item[1][k][2]
                                        hashValue = item[1][k][3]
                                        numberValue = item[1][k][4]
                                        (connectionData[index],storeData) = add_connectionData(connectionData[index],
                                                                            lengthValue,directValue,
                                                                            hashValue,numberValue)
                                        if storeData:
                                            linkValue = item[1][k][0]
                                            linkData = storageHandler["handler"].root.data[linkValue:linkValue+lengthValue]
                                            connectionData[index][1][hashValue]["data"] = linkData
                        for j in range(len(connectionData)):
                            fromLink = i + j
                            if connectionData[j][0]<=maximumConnectionLength:
                                for hashValue in connectionData[j][1].keys():
                                    entry = connectionData[j][1][hashValue]
                                    if entry["status"]:
                                        tableConnectionsRow = tableConnections.row
                                        tableConnectionsRow["ckmerLink"] = fromLink
                                        tableConnectionsRow["dataLink"] = tableDataCounter
                                        tableConnectionsRow["length"] = entry["length"]
                                        tableConnectionsRow["number"] = entry["number"]
                                        tableConnectionsRow["direct"] = entry["direct"]
                                        tableConnectionsRow.append()
                                        tableData.append(entry["data"])
                                        tableDataCounter+=len(entry["data"])
                                        
                        #initialise
                        pairedData = [{} for j in range(min(mergeEnd+1-i,stepSizeStorage,(numberOfKmers-i)))]
                        for storageHandler in storageHandlers:
                            rowData = storageHandler["handler"].root.paired[i:min(mergeEnd+1,i+stepSizeStorage)]
                            for index, item in enumerate(rowData):
                                fromLink = i + index
                                for j in range(item[0]):
                                    toLink = item[1][j][0]
                                    number = item[1][j][1]
                                    pairedData[index][toLink] = pairedData[index].get(toLink,0)+number
                                    
                                #now check other data
                                if storageHandler["otherRow"]:
                                    if storageHandler["otherRow"]["fromLink"]<fromLink:
                                        storageHandler["otherRow"] = None
                                        for row in storageHandler["otherIterator"]:
                                            if row["fromLink"]>=fromLink:
                                                if row["fromLink"]<=mergeEnd:
                                                    storageHandler["otherRow"] = row
                                                break
                                    if storageHandler["otherRow"] and storageHandler["otherRow"]["fromLink"]==fromLink:
                                        toLink = storageHandler["otherRow"]["toLink"]
                                        pairedData[index][toLink] = pairedData[index].get(toLink,0)+1
                                        for row in storageHandler["otherIterator"]:
                                            if row["fromLink"]==fromLink:
                                                storageHandler["otherRow"] = row
                                                toLink = storageHandler["otherRow"]["toLink"]
                                                pairedData[index][toLink] = pairedData[index].get(toLink,0)+1
                                            elif row["fromLink"]<=mergeEnd:
                                                storageHandler["otherRow"] = row
                                                break
                                            else:
                                                storageHandler["otherRow"] = None
                                                break
                                
                        for j in range(len(pairedData)):
                            fromLink = i + j
                            if len(pairedData[j])>maximumConnectionLength:
                                for toLink in pairedData[j].keys():
                                    tableDeletePairedRow = tableDeletePaired.row
                                    tableDeletePairedRow["fromLink"] = toLink
                                    tableDeletePairedRow["toLink"] = fromLink
                                    tableDeletePairedRow.append()
                            else:
                                for toLink in pairedData[j].keys():
                                    number = pairedData[j][toLink]
                                    tablePairedRow = tablePaired.row
                                    tablePairedRow["fromLink"] = fromLink
                                    tablePairedRow["toLink"] = toLink
                                    tablePairedRow["number"] = number
                                    tablePairedRow.append()          
                                    
                #finished
                pytablesStorageWorker.flush()
                
                
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("merges: shared memory of {} MB used".format(round(shm.size/1048576)))
        try:
            #get shared memory
            shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
            shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maximumFrequency)).type
            kmer_properties = np.ndarray((numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                               ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm.buf)
            #handle queue
            while not shutdown_event.is_set():
                try:
                    item = queue_ranges.get(block=True, timeout=1)
                    if item==None:
                        break
                    elif isinstance(item,tuple):
                        mergeStart = item[0]
                        mergeNumber = item[1]
                        #create storage
                        numberLength = len(str(numberOfKmers))
                        mergeEnd = min(numberOfKmers,mergeStart+mergeNumber)-1
                        pytablesFileRange = (filenameBase+"_tmp_connections_merge_"+
                                             str(mergeStart).zfill(numberLength)+"_"+
                                        str(mergeEnd).zfill(numberLength)+".process.h5")
                        if os.path.exists(pytablesFileRange):
                            os.remove(pytablesFileRange)
                        with tables.open_file(pytablesFileRange, mode="a") as pytablesStorageRange:
                            Storage.create_merge_storage(pytablesStorageRange, numberOfKmers, 
                                                         maximumFrequency, maximumConnectionLength)
                        #merge
                        merge_direct_storage(pytablesFileRange, storageDirectFiles, mergeStart, 
                                             mergeNumber, numberOfKmers, kmer_properties)
                        merge_connection_storage(pytablesFileRange, storageConnectionFiles, mergeStart, 
                                             mergeNumber, numberOfKmers, maximumConnectionLength)
                        #now the file can be released for final merge
                        queue_merges.put(pytablesFileRange) 
                except Empty:
                    continue
        finally:
            shm.close()
        logger.debug("index: shared memory released")
            

    
    """
    Combine merged files: this should be relatively easy, handled by single process
    """
    def combine_merges(mergeFiles,pytablesStorage,numberOfKmers,maximumFrequency,maximumConnectionLength):
        
        nCycle = 0
        nReversal = 0
        nDirect = 0
        nConnections = 0
        nData = 0
        nPaired = 0

        #get dimensions
        sortedMergeFiles = [{"filename": mergeFile} for mergeFile in sorted(mergeFiles)]
        for i in range(len(sortedMergeFiles)):
            with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                sortedMergeFiles[i]["ncycle"] = pytables_merge.root.cycle.shape[0]
                sortedMergeFiles[i]["nreversal"] = pytables_merge.root.reversal.shape[0]
                sortedMergeFiles[i]["ndirect"] = pytables_merge.root.direct.shape[0]
                sortedMergeFiles[i]["ndeleteDirect"] = pytables_merge.root.deleteDirect.shape[0]
                sortedMergeFiles[i]["nconnections"] = pytables_merge.root.connections.shape[0]
                sortedMergeFiles[i]["ndata"] = pytables_merge.root.data.shape[0]
                sortedMergeFiles[i]["npaired"] = pytables_merge.root.paired.shape[0]
                sortedMergeFiles[i]["ndeletePaired"] = pytables_merge.root.deletePaired.shape[0]
                nCycle += sortedMergeFiles[i]["ncycle"]
                nReversal += sortedMergeFiles[i]["nreversal"]
                nDirect += sortedMergeFiles[i]["ndirect"]
                nConnections += sortedMergeFiles[i]["nconnections"]
                nData = max(nData,sortedMergeFiles[i]["ndata"])
                nPaired  += sortedMergeFiles[i]["npaired"]

        #create temporary storage with dimensions
        Storage.create_merge_storage(pytablesStorage, numberOfKmers, maximumFrequency, maximumConnectionLength, 
                                     nCycle, nReversal, nDirect, nConnections, nPaired)
        tableCycle = pytablesStorage.root.cycle
        tableReversal = pytablesStorage.root.reversal
        tableDirect = pytablesStorage.root.direct
        tableDeleteDirect = pytablesStorage.root.deleteDirect
        tableConnections = pytablesStorage.root.connections
        tableData = pytablesStorage.root.data
        tablePaired = pytablesStorage.root.paired
        tableDeletePaired = pytablesStorage.root.deletePaired
        
        stepSizeStorage = 100 #FIX THIS BACK TO 100000
        maximumCycleLength = 0
        maximumCycleNumber = 0
        maximumReversalLength = 0
        maximumReversalNumber = 0
        maximumDirectDistance = 0
        maximumDirectNumber = 0
        maximumConnectionsNumber = 0
        maximumConnectionsLength = 0
        maximumPairedNumber = 0
        
        connectionsDataCounter = 0
        
        #merge delete entries
        for i in range(len(sortedMergeFiles)):
            with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                tableDeleteDirect.append(pytables_merge.root.deleteDirect[:])
                tableDeletePaired.append(pytables_merge.root.deletePaired[:])
        pytablesStorage.flush()
        tableDeleteDirect.cols.fromLink.create_csindex()
        tableDeletePaired.cols.fromLink.create_csindex()
        pytablesStorage.flush()
                
        for i in range(len(sortedMergeFiles)):
            with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                for j in range(0,sortedMergeFiles[i]["ncycle"],stepSizeStorage):
                    reducedCycleDataBlock = pytables_merge.root.cycle[j:j+stepSizeStorage]
                    maximumCycleNumber = max(maximumCycleNumber,max(reducedCycleDataBlock["number"]))
                    maximumCycleLength = max(maximumCycleLength,max(reducedCycleDataBlock["minimumLength"]))
                    tableCycle.append(reducedCycleDataBlock)
                for j in range(0,sortedMergeFiles[i]["nreversal"],stepSizeStorage):
                    reducedReversalDataBlock = pytables_merge.root.reversal[j:j+stepSizeStorage]
                    maximumReversalNumber = max(maximumReversalNumber,max(reducedReversalDataBlock["number"]))
                    maximumReversalLength = max(maximumReversalLength,max(reducedReversalDataBlock["minimumLength"]))
                    tableReversal.append(reducedReversalDataBlock)
                for j in range(0,sortedMergeFiles[i]["ndirect"],stepSizeStorage):
                    directDataBlock = pytables_merge.root.direct[j:j+stepSizeStorage]
                    deleteDataBlock = tableDeleteDirect.read_where(
                        "(fromLink>={}) & (fromLink<={})".format(
                        directDataBlock[0]["fromLink"],directDataBlock[-1]["fromLink"]))
                    deletableIndices = np.where(
                        np.in1d(directDataBlock[["fromLink","fromDirection","toLink","toDirection"]],
                                deleteDataBlock[["fromLink","fromDirection","toLink","toDirection"]]))[0]
                    reducedDirectDataBlock = np.delete(directDataBlock,deletableIndices,0)
                    maximumDirectNumber = max(maximumDirectNumber,max(reducedDirectDataBlock["number"]))
                    maximumDirectDistance = max(maximumDirectDistance,max(reducedDirectDataBlock["distance"]))
                    #add reduced direct data
                    tableDirect.append(reducedDirectDataBlock)
                for j in range(0,sortedMergeFiles[i]["nconnections"],stepSizeStorage):
                    reducedConnectionsBlock = pytables_merge.root.connections[j:j+stepSizeStorage]
                    for k in range(len(reducedConnectionsBlock)):
                        reducedConnectionsBlock[k][1]+=connectionsDataCounter
                    maximumConnectionsNumber = max(maximumConnectionsNumber,max(reducedConnectionsBlock["number"]))
                    maximumConnectionsLength = max(maximumConnectionsLength,max(reducedConnectionsBlock["length"]))
                    tableConnections.append(reducedConnectionsBlock)
                for j in range(0,sortedMergeFiles[i]["ndata"],stepSizeStorage):
                    reducedDataBlock = pytables_merge.root.data[j:j+stepSizeStorage]
                    tableData.append(reducedDataBlock)
                    connectionsDataCounter+=len(reducedDataBlock)
                for j in range(0,sortedMergeFiles[i]["npaired"],stepSizeStorage):
                    pairedDataBlock = pytables_merge.root.paired[j:j+stepSizeStorage]
                    deleteDataBlock = tableDeletePaired.read_where(
                        "(fromLink>={}) & (fromLink<={})".format(
                        pairedDataBlock[0]["fromLink"],pairedDataBlock[-1]["fromLink"]))
                    deletableIndices = np.where(
                        np.in1d(pairedDataBlock[["fromLink","toLink"]],
                                deleteDataBlock[["fromLink","toLink"]]))[0]
                    reducedPairedDataBlock = np.delete(pairedDataBlock,deletableIndices,0)
                    maximumPairedNumber = max(maximumPairedNumber,max(reducedPairedDataBlock["number"]))
                    #add reduced paired data
                    tablePaired.append(reducedPairedDataBlock)

        tableCycle.attrs["maximumNumber"] = maximumCycleNumber
        tableCycle.attrs["maximumLength"] = maximumCycleLength
        tableReversal.attrs["maximumNumber"] = maximumReversalNumber
        tableReversal.attrs["maximumLength"] = maximumReversalLength
        tableDirect.attrs["maximumDistance"] = maximumDirectDistance
        tableDirect.attrs["maximumNumber"] = maximumDirectNumber
        tableConnections.attrs["maximumNumber"] = maximumConnectionsNumber
        tableConnections.attrs["maximumLength"] = maximumConnectionsLength
        pytablesStorage.flush()
            
    """
    Store everything in the final database, handled by single process
    """    
    def store_merged_connections(h5file,pytablesStorage,numberOfKmers,minimumFrequency):
        logger = logging.getLogger(__name__)     
        frequencyHistogram = {"distance": {}}
        
        def classifyDirectProblematic(dataBlock,directGood,directCoverage,directProblematic):
            dataBlock = np.array(dataBlock, dtype=object)                       
            if max(dataBlock[:,4])==0:
                directGood+=len(dataBlock)
            elif len(dataBlock)==1:
                dataBlock[:,4] = 0
                directGood+=len(dataBlock)
            elif max(dataBlock[:,2])<minimumFrequency:
                dataBlock[:,4] = 1
                directCoverage+=len(dataBlock)
            else:
                dataBlock[:,4] = 2
                directProblematic+=len(dataBlock)
            number = sum(dataBlock[:,2])
            distinct = len(dataBlock)
            dataBlock = tuple([tuple(e) for e in dataBlock])
            return (dataBlock,directGood,directCoverage,directProblematic,distinct,number)
        
        stepSizeStorage = 100000
        dsCkmer = h5file["/split/ckmer"]
        
        #cycles
        numberOfCycles = pytablesStorage.root.cycle.shape[0]
        maximumNumber = pytablesStorage.root.cycle.attrs.maximumNumber
        maximumLength = pytablesStorage.root.cycle.attrs.maximumLength
        dtypeCycleList=[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                         ("minimumLength",haplotyping.index.Database.getUint(maximumLength)),
                         ("number",haplotyping.index.Database.getUint(maximumNumber))]
        dtCycle=np.dtype(dtypeCycleList)
        dsCycle=h5file["/relations/"].create_dataset("cycle",(numberOfCycles,), 
                                                     dtype=dtCycle, chunks=None)
        maximumCycleLength = 0
        maximumCycleNumber = 0
        for i in range(0,numberOfCycles,stepSizeStorage):
            stepData = pytablesStorage.root.cycle[i:i+stepSizeStorage]                    
            dsCycle[i:i+stepSizeStorage] = stepData
            for row in stepData:
                ckmerRow = dsCkmer[row[0]]
                ckmerRow[7] = row[2]
                maximumCycleLength = max(maximumCycleLength,row[1])
                maximumCycleNumber = max(maximumCycleNumber,row[2])
                dsCkmer[row[0]] = ckmerRow    
        h5file["/config/"].attrs["numberOfCycles"]=numberOfCycles        
        h5file["/config/"].attrs["maximumCycleLength"]=maximumCycleLength
        h5file["/config/"].attrs["maximumCycleNumber"]=maximumCycleNumber
        logger.info("store "+str(numberOfCycles)+" cycles")
        #reversals
        numberOfReversals = pytablesStorage.root.reversal.shape[0]
        maximumNumber = pytablesStorage.root.reversal.attrs.maximumNumber
        maximumLength = pytablesStorage.root.reversal.attrs.maximumLength
        dtypeReversalList=[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                         ("minimumLength",haplotyping.index.Database.getUint(maximumLength)),
                         ("number",haplotyping.index.Database.getUint(maximumNumber))]
        dtReversal=np.dtype(dtypeReversalList)
        dsReversal=h5file["/relations/"].create_dataset("reversal",(numberOfReversals,), 
                                                        dtype=dtReversal, chunks=None)
        maximumReversalLength = 0
        maximumReversalNumber = 0
        for i in range(0,numberOfReversals,stepSizeStorage):
            stepData = pytablesStorage.root.reversal[i:i+stepSizeStorage]
            dsReversal[i:i+stepSizeStorage] = stepData
            for row in stepData:
                ckmerRow = dsCkmer[row[0]]
                ckmerRow[8] = row[2]
                maximumReversalLength = max(maximumReversalLength,row[1])
                maximumReversalNumber = max(maximumReversalNumber,row[2])
                dsCkmer[row[0]] = ckmerRow     
        h5file["/config/"].attrs["numberOfReversals"]=numberOfReversals
        h5file["/config/"].attrs["maximumReversalLength"]=maximumReversalLength
        h5file["/config/"].attrs["maximumReversalNumber"]=maximumReversalNumber
        logger.info("store "+str(numberOfReversals)+" reversals")
        #direct relations
        numberOfDirectRelations = pytablesStorage.root.direct.shape[0]
        maximumNumber = pytablesStorage.root.direct.attrs.maximumNumber
        maximumDistance = pytablesStorage.root.direct.attrs.maximumDistance
        dtypeDirectList=[("from",[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                                  ("direction","S1")]),
                         ("to",[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                                  ("direction","S1")]),
                         ("number",haplotyping.index.Database.getUint(maximumNumber)),
                         ("distance",haplotyping.index.Database.getUint(maximumDistance)),
                         ("problematic","uint8")]
        dtDirect=np.dtype(dtypeDirectList)
        dsDirect=h5file["/relations/"].create_dataset("direct",(numberOfDirectRelations,), 
                                                      dtype=dtDirect, chunks=None)
        #connections
        numberOfConnections = pytablesStorage.root.connections.shape[0]
        numberOfData = pytablesStorage.root.data.shape[0]
        numberOfPaired = pytablesStorage.root.paired.shape[0]
        maximumNumber = pytablesStorage.root.connections.attrs.maximumNumber
        maximumLength = pytablesStorage.root.connections.attrs.maximumLength
        dtypeConnectionsList=[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                              ("dataLink",haplotyping.index.Database.getUint(numberOfData)),
                              ("length",haplotyping.index.Database.getUint(maximumLength)),
                              ("direct","uint8"),
                              ("number",haplotyping.index.Database.getUint(maximumNumber))
                             ]
        dtConnections=np.dtype(dtypeConnectionsList)
        dsConnections=h5file["/connections/"].create_dataset("index",(numberOfConnections,), 
                                                      dtype=dtConnections, chunks=None)
        dsData=h5file["/connections/"].create_dataset("data",(numberOfData,), 
                          dtype=haplotyping.index.Database.getUint(numberOfKmers), chunks=None)
        dtypePairedList=[("fromLink",haplotyping.index.Database.getUint(numberOfKmers)),
                         ("toLink",haplotyping.index.Database.getUint(numberOfKmers)),
                         ("number",haplotyping.index.Database.getUint(maximumNumber))
                        ]
        dtPaired=np.dtype(dtypePairedList)
        dsPaired=h5file["/connections/"].create_dataset("paired",(numberOfPaired,), 
                                                      dtype=dtPaired, chunks=None)
        #process direct relations
        directCounter = 0
        previousStepData = []
        directGood = 0
        directCoverage = 0
        directProblematic = 0
        for i in range(0,numberOfDirectRelations,stepSizeStorage):
            stepData = pytablesStorage.root.direct[i:i+stepSizeStorage]
            #move a bit to get all fromLinks in the selection
            if len(previousStepData)>0:
                stepData = np.concatenate((previousStepData, stepData))
            if (i+stepSizeStorage-1)<numberOfDirectRelations:
                lastFromLink = stepData[-1]["fromLink"]
                previousStepData = np.take(stepData, np.where(stepData["fromLink"] == lastFromLink))[0]
                stepData = np.delete(stepData, np.where(stepData["fromLink"] == lastFromLink))
            else:
                previousStepData = [] 
            #get relevant ckmer data
            firstFromLink = stepData[0]["fromLink"]
            lastFromLink = stepData[-1]["fromLink"]
            ckmerStepData = dsCkmer[firstFromLink:(lastFromLink+1)]
            #now create appropiate final selection
            newStepData = []
            newStepDataBlock = []
            previousFromLink = None
            previousFromDirection = None
            previousUpdateLink = None
            updateLink = None
            for row in stepData:
                if not ((previousFromLink==row["fromLink"]) and (previousFromDirection==row["fromDirection"])):
                    if not (previousFromLink == None):
                        updateLink = previousFromLink-firstFromLink
                        (newStepDataBlock,directGood,directCoverage,
                         directProblematic,distinct,number) = classifyDirectProblematic(
                            newStepDataBlock,directGood,directCoverage,directProblematic)
                        #update ckmer-data
                        if not updateLink == previousUpdateLink:
                            ckmerStepData[updateLink][4][0] = directCounter + len(newStepData)
                        if previousFromDirection == b"l":
                            ckmerStepData[updateLink][4][1] = (distinct,number,)
                        elif previousFromDirection == b"r":
                            ckmerStepData[updateLink][4][2] = (distinct,number,)
                        else:
                            raise Exception("unexpected direction")
                        newStepData.extend(newStepDataBlock)
                        newStepDataBlock = []
                    previousFromLink = row["fromLink"]
                    previousFromDirection = row["fromDirection"]
                    previousUpdateLink = updateLink
                newStepDataBlock.append([(row["fromLink"],row["fromDirection"],),
                                    (row["toLink"],row["toDirection"],),
                                    row["number"],row["distance"],row["problematic"]])   
                frequencyHistogram["distance"][row["distance"]] = (
                    frequencyHistogram["distance"].get(row["distance"],0)+row["number"])
            updateLink = previousFromLink-firstFromLink
            (newStepDataBlock,directGood,directCoverage,
             directProblematic,distinct,number) = classifyDirectProblematic(
                        newStepDataBlock,directGood,directCoverage,directProblematic)
            #update ckmer-data
            if not updateLink == previousUpdateLink:
                ckmerStepData[updateLink][4][0] = directCounter + len(newStepData)
            if previousFromDirection == b"l":
                ckmerStepData[updateLink][4][1] = (distinct,number,)
            elif previousFromDirection == b"r":
                ckmerStepData[updateLink][4][2] = (distinct,number,)
            else:
                raise Exception("unexpected direction")
            newStepData.extend(newStepDataBlock)
            #store all
            dsDirect[directCounter:directCounter+len(newStepData)] = newStepData
            dsCkmer[firstFromLink:(lastFromLink+1)] = ckmerStepData
            directCounter+=len(newStepData)

        logger.info("store {} direct connections, {} low coverage, {} problematic".format(
            numberOfDirectRelations,directCoverage,directProblematic))

        #process connections
        connectionsCounter = 0
        dataCounter = 0
        pairedCounter = 0
        ckmerConnectionsIndexLink = np.zeros(dsCkmer.len(), dtype=int)
        ckmerConnectionsIndexCounter = np.zeros(dsCkmer.len(), dtype=int)
        ckmerConnectionsCounter = np.zeros(dsCkmer.len(), dtype=int)
        ckmerPairedLink = np.zeros(dsCkmer.len(), dtype=int)
        ckmerPairedCounter = np.zeros(dsCkmer.len(), dtype=int)
        for i in range(0,numberOfConnections,stepSizeStorage):
            stepData = pytablesStorage.root.connections[i:i+stepSizeStorage]  
            for j in range(i,i+len(stepData)):    
                if ckmerConnectionsIndexCounter[stepData[j-i][0]]==0:
                    ckmerConnectionsIndexLink[stepData[j-i][0]] = j
                    ckmerConnectionsIndexCounter[stepData[j-i][0]] = 1
                else:
                    ckmerConnectionsIndexCounter[stepData[j-i][0]]+=1
            dsConnections[connectionsCounter:connectionsCounter+len(stepData)] = stepData
            connectionsCounter+=len(stepData)
        for i in range(0,numberOfData,stepSizeStorage):
            stepData = pytablesStorage.root.data[i:i+stepSizeStorage]   
            for row in stepData:
                ckmerConnectionsCounter[row]+=1
            dsData[dataCounter:dataCounter+len(stepData)] = stepData
            dataCounter+=len(stepData)
        for i in range(0,numberOfPaired,stepSizeStorage):
            stepData = pytablesStorage.root.paired[i:i+stepSizeStorage]   
            for j in range(i,i+len(stepData)):
                if ckmerPairedCounter[stepData[j-i][0]]==0:
                    ckmerPairedLink[stepData[j-i][0]] = j
                    ckmerPairedCounter[stepData[j-i][0]] = 1
                else:
                    ckmerPairedCounter[stepData[j-i][0]]+=1
            dsPaired[pairedCounter:pairedCounter+len(stepData)] = stepData
            pairedCounter+=len(stepData)
        for i in range(0,dsCkmer.len(),stepSizeStorage):
            stepData = dsCkmer[i:i+stepSizeStorage]  
            for j in range(i,i+len(stepData)):                
                stepData[j-i][5] = (ckmerConnectionsIndexLink[j],ckmerConnectionsIndexCounter[j],ckmerConnectionsCounter[j],)
                stepData[j-i][6] = (ckmerPairedLink[j],ckmerPairedCounter[j],)
            dsCkmer[i:i+len(stepData)] = stepData
        
        logger.info("store {} indirect connections with {} data entries and {} paired".format(
            numberOfConnections,numberOfData,numberOfPaired))

        # store histogram direct distances
        maximumDistance = max(frequencyHistogram["distance"].keys())
        maximumNumber = max(frequencyHistogram["distance"].values())
        dtypeFrequencyHistogramDistanceList=[
                    ("distance",haplotyping.index.Database.getUint(maximumDistance)),
                    ("number",haplotyping.index.Database.getUint(maximumNumber)),]
        dtFrequencyHistogramDistance=np.dtype(dtypeFrequencyHistogramDistanceList)
        dsFrequencyHistogramDistance=h5file["/histogram/"].create_dataset("distance",
              (len(frequencyHistogram["distance"]),), dtype=dtFrequencyHistogramDistance, chunks=None)
        logger.info("store {} entries direct distances histogram".format(
            len(frequencyHistogram["distance"])))
        #store histogram
        dsFrequencyHistogramDistance[0:len(frequencyHistogram["distance"])] = list(
            sorted(frequencyHistogram["distance"].items()))



