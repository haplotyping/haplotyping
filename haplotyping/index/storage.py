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
    
    def create_merge_storage(pytables_storage, numberOfKmers, nCycle=None, nReversal=None, nDirect=None, deleteDirect=True):
        tableCycle = pytables_storage.create_table(pytables_storage.root, 
                            "cycle",{
                            "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                            "minimumLength": tables.UInt8Col(pos=1),
                            "number": tables.UInt8Col(pos=2),
                        }, "Cycles", expectedrows=nCycle)
        tableReversal = pytables_storage.create_table(pytables_storage.root, 
                            "reversal",{
                            "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                            "minimumLength": tables.UInt8Col(pos=1),
                            "number": tables.UInt8Col(pos=2),
                        }, "Reversals", expectedrows=nReversal)
        tableDirect = pytables_storage.create_table(pytables_storage.root, 
                            "direct",{
                            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                            "fromDirection": tables.StringCol(itemsize=1,pos=1),
                            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,2),
                            "toDirection": tables.StringCol(itemsize=1,pos=3),
                            "number": tables.UInt16Col(pos=4),
                            "distance": tables.UInt16Col(pos=5),
                            "splitDirection": tables.StringCol(itemsize=1,pos=6),
                            "reverseBase": haplotyping.index.Database.getTablesUint(numberOfKmers,7),
                            "forwardBase": haplotyping.index.Database.getTablesUint(numberOfKmers,8),
                            "problematic": tables.UInt8Col(pos=9),
                        }, "Direct relations", expectedrows=nDirect)
        if deleteDirect:
            tableDeleteDirect = pytables_storage.create_table(pytables_storage.root, 
                            "deleteDirect",{
                            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                            "fromDirection": tables.StringCol(itemsize=1,pos=1),
                            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,2),
                            "toDirection": tables.StringCol(itemsize=1,pos=3),
                        }, "Incorrect direct relations", expectedrows=nDirect)
    
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
        logger.debug("automaton: fsm loaded from {} bytes".format(os.stat(automatonFile).st_size))
            
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
            return matchesList
        
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("index: shared memory of {} bytes used".format(shm.size))
        try:
            while not shutdown_event.is_set():
                try:
                    item = queue_index.get(block=True, timeout=1)
                    if item==None:
                        break
                    elif isinstance(item,tuple):
                        if len(item)==1:
                            matches = compute_matches(item[0][0],item[0][1])
                            if len(matches)==0 or (len(matches)==1 and len(matches[0])==1):
                                pass
                            else:
                                queue_matches.put((matches,))
                        elif len(item)==2:
                            sequence0 = item[0][0]
                            sequence1 = item[1][0]
                            matches0 = compute_matches(item[0][0],item[0][1])
                            matches1 = compute_matches(item[1][0],item[1][1])
                            if len(matches0)==0 and len(matches1)==0:
                                pass
                            elif len(matches0)==0:
                                if len(matches1)==1 and len(matches1[0])==1:
                                    pass
                                else:
                                    queue_matches.put((matches1,))
                            elif len(matches1)==0:
                                if len(matches0)==1 and len(matches0[0])==1:
                                    pass
                                else:
                                    queue_matches.put((matches0,))
                            else:
                                queue_matches.put((matches0,matches1,))
                except Empty:
                    continue
        finally:
            shm.close()
        logger.debug("index: shared memory released")
        queue_finished.put("index")
            
    def worker_matches(shutdown_event,queue_matches,queue_connections,queue_storage,queue_finished,
                       filenameBase,numberOfKmers,arrayNumber):
        
        logger = logging.getLogger(__name__)
        try:
            curr_proc = current_process()
            pytablesFile = filenameBase+"_tmp_direct_"+str(curr_proc.name)+".process.h5"
            if os.path.exists(pytablesFile):
                os.remove(pytablesFile)

            with tables.open_file(pytablesFile, mode="a") as pytables_storage:
                
                logger.debug("matches: store direct connections")

                dtype = [("cycle", [("number", np.uint16),("minimum", "uint16")]),
                         ("reversal", [("number", np.uint16),("minimum", "uint16")]),
                         ("direct", [("d"+str(i),[("type","uint8"),
                                                  ("link",haplotyping.index.Database.getUint(numberOfKmers)),
                                                  ("distance","uint16"),
                                                  ("number","uint16")]) 
                                     for i in range(arrayNumber)]),]
                connections = np.ndarray((numberOfKmers,), dtype=dtype, order="C")
                connections.fill(((0,0,),(0,0,),tuple((0,0,0,0,) for i in range(arrayNumber))))
                
                logger.debug("created memory storage for direct connections: {} bytes".format(connections.nbytes))
                
                tableDirectOther = pytables_storage.create_table(pytables_storage.root, 
                    "directOther",{
                    "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                    "fromDirection": tables.UInt8Col(pos=1),
                    "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,2),
                    "toDirection": tables.UInt8Col(pos=3),
                    "number": tables.UInt16Col(pos=4),
                    "distance": tables.UInt16Col(pos=5),
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
                    numberFilled = arrayNumber
                    connectionType = 1 + fromDirection + (2*toDirection)
                    directRow = connections[fromLink][2]
                    for i in range(arrayNumber):
                        if directRow[i][0]==0:
                            numberFilled = i
                            break
                        elif (directRow[i][1]==toLink and directRow[i][0]==connectionType 
                              and directRow[i][2]==distance):
                            directRow[i][3]+=1 
                            return
                    if numberFilled<arrayNumber:
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
                        if 2*len(sameFromDirectionDistances)>=arrayNumber:
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
                        elif 2*len(otherFromDirectionDistances)>arrayNumber:
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
                    return (links,max(0,endPos-startPos),len(matches)==1,)

                while not shutdown_event.is_set():
                    try:
                        item = queue_matches.get(block=True, timeout=1)
                        if item==None:
                            break
                        elif isinstance(item,tuple):
                            if len(item)==1:
                                matchesList = item[0]
                                (links,length,directConnected,) = process_matches(matchesList)
                                if len(links)>2:
                                    queue_connections.put(((links,length,directConnected,)))
                            elif len(item)==2:
                                matchesList0 = item[0]
                                matchesList1 = item[1]
                                (links0,length0,directConnected0,) = process_matches(matchesList0)
                                (links1,length1,directConnected1,) = process_matches(matchesList1)
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
                pytables_storage.flush()
                tableDirectOther.cols.fromLink.create_csindex()
                pytables_storage.flush()
                tableDirect = pytables_storage.create_table(pytables_storage.root, 
                                              name="direct", obj=connections, expectedrows=numberOfKmers)
                logger.debug("matches: memory storage saved, other connections indexed".format(pytablesFile))
                pytables_storage.flush()
                
            #now the file can be released for merge
            queue_storage.put(pytablesFile) 
        finally:
            pass
        queue_finished.put("matches")
            
            
    def worker_connections(shutdown_event,queue_connections,queue_storage,queue_finished,
                       filenameBase,numberOfKmers,arrayNumber,maxFrequency,shm_name):
        
        logger = logging.getLogger(__name__)
        
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("connections: shared memory of {} bytes used".format(shm.size))
            
        try:
            curr_proc = current_process()
            pytablesFile = filenameBase+"_tmp_connections_"+str(curr_proc.name)+".process.h5"
            if os.path.exists(pytablesFile):
                os.remove(pytablesFile)

            with tables.open_file(pytablesFile, mode="a") as pytables_storage:
                
                logger.debug("connections: store additional connections")

                dtype = [("number", haplotyping.index.Database.getUint(arrayNumber+1),
                         ("connections", [("c"+str(i),[
                                            ("length","uint16"),
                                            ("direct","uint8"),
                                            ("link",haplotyping.index.Database.getUint(numberOfKmers)),
                                            ("hash","int64")]) 
                                          for i in range(arrayNumber)]),]
                connections = np.ndarray((numberOfKmers,), dtype=dtype, order="C")
                connections.fill((0,tuple((0,0,0,0,) for i in range(arrayNumber))))
                
                logger.debug("created memory storage for additional connections: {} bytes".format(connections.nbytes))
                
                tableConnectionsData = pytables_storage.create_vlarray(pytables_storage.root, "data", 
                               haplotyping.index.Database.getTablesUintAtom(numberOfKmers), "Data")

                mapSplitDirection = [{b"l":"l", b"r":"r", b"b":"b"},{b"l":"r", b"r":"l", b"b":"b"}]
                
                def store_connections(connectionNumber,connectionsList,length,direct,paired=False,pairLink=0):
                    fullConnectionsList = []
                    for (link,direction) in connectionsList:
                        props = kmer_properties[link] 
                        newEntry = (link,mapSplitDirection[direction][props[0]],props[1],)
                        #contract pairs with equal dosage (rightsplitter -> leftsplitter)
                        if newEntry[1]=="l" and (len(fullConnectionsList)>0) and fullConnectionsList[-1][1]=="r":
                            fullConnectionsList[-1] = (fullConnectionsList[-1],newEntry,)
                        else:
                            fullConnectionsList.append(newEntry)
                    #select entries with (potential) minimal dosage
                    potentialMinimum = None
                    potentialDirection = None
                    finalConnectionsList = []       
                    for entry in fullConnectionsList:
                        newDirection = entry[0][1] if len(entry)==2 else entry[1]
                        if potentialMinimum==None:
                            potentialMinimum = entry
                            potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                        else:
                            if potentialDirection=="l":
                                if newDirection=="l":
                                    pass
                                else:
                                    finalConnectionsList.append(entry)
                                    potentialMinimum = entry
                                    potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                            elif potentialDirection=="r":
                                if newDirection=="l":
                                    raise Exception("can't happen")
                                else:
                                    potentialMinimum = entry
                                    potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                            else:
                                if newDirection=="l":
                                    pass
                                else:
                                    finalConnectionsList.append(entry)
                                    potentialMinimum = entry
                                    potentialDirection = entry[-1][1] if len(entry)==2 else entry[1]
                            
                    if not potentialMinimum==None:
                        finalConnectionsList.append(entry)
                    print("{}\t{}\t{}%".format(len(connectionsList),len(finalConnectionsList),
                                              int(100*len(finalConnectionsList)/len(connectionsList))))
                    return connectionNumber+1
                    
                #get shared memory
                shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
                shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maxFrequency)).type
                kmer_properties = np.ndarray((numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                                   ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm.buf)
            
                connectionNumber = 0
                while not shutdown_event.is_set():
                    try:
                        item = queue_connections.get(block=True, timeout=1)
                        if item==None:
                            break
                        elif isinstance(item,tuple):
                            if 1==1:
                                #for now, just ignore
                                pass
                            elif len(item)==1:
                                if len(item[0][0])>2:
                                    connectionNumber = store_connections(connectionNumber,
                                                         item[0][0],item[0][1],item[0][2])
                            elif len(item)==2:
                                if len(item[0][0])==0:
                                    if len(item[1][0])>0:
                                        connectionNumber = store_connections(connectionNumber,
                                                         item[1][0],item[1][1],item[1][2])
                                elif len(item[1][0])==0:
                                    if len(item[0][0])>0:
                                        connectionNumber = store_connections(connectionNumber,
                                                         item[0][0],item[0][1],item[0][2])
                                else:
                                    pairedLink0 = connectionNumber
                                    pairedLink1 = connectionNumber + 1
                                    connectionNumber = store_connections(connectionNumber,
                                                         item[0][0],item[0][1],item[0][2],True,pairedLink1)
                                    connectionNumber = store_connections(connectionNumber,
                                                         item[1][0],item[1][1],item[1][2],True,pairedLink0)
                    except Empty:
                        continue   
                
                pytables_storage.flush()
                
                tableConnection = pytables_storage.create_table(pytables_storage.root, 
                                              name="direct", obj=connections, expectedrows=numberOfKmers)
                logger.debug("matches: memory storage connections saved".format(pytablesFile))
                pytables_storage.flush()
                
            #now the file can be released for merge
            queue_storage.put(pytablesFile) 
        finally:
            shm.close()
            logger.debug("index: shared memory released")

        queue_finished.put("connections")
                
                
    def worker_merges(shutdown_event,queue_ranges,queue_merges,storageFiles,
                      filenameBase,numberOfKmers,arrayNumber,maxFrequency,minFrequency,shm_name):
        
        logger = logging.getLogger(__name__)
        
        def merge_storage(filenameBase, storageFiles, mergeStart, mergeNumber, numberOfKmers, arrayNumber, kmer_properties):
            
            numberLength = len(str(numberOfKmers))
            mergeEnd = min(numberOfKmers,mergeStart+mergeNumber)-1
            pytablesFile = (filenameBase+"_tmp_direct_merge_"+str(mergeStart).zfill(numberLength)+"_"+
                            str(mergeEnd).zfill(numberLength)+".process.h5")
            if os.path.exists(pytablesFile):
                os.remove(pytablesFile)
                
            with tables.open_file(pytablesFile, mode="a") as pytables_storage:
                #create temporary storage
                Storage.create_merge_storage(pytables_storage, numberOfKmers)
                tableCycle = pytables_storage.root.cycle
                tableReversal = pytables_storage.root.reversal
                tableDirect = pytables_storage.root.direct
                tableDeleteDirect = pytables_storage.root.deleteDirect

                with ExitStack() as stack:
                    #get handlers
                    storageHandlers = [{"handler": stack.enter_context(tables.open_file(fname, mode="r"))} 
                                               for fname in storageFiles]
                    
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
                            directDataEntry[fromDirection][distance][toLink][toDirection]["number"] += number
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
                                    if number>=minFrequency:
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
                                    if number>=minFrequency:
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
                                    maxDistanceFrequency = max(max(distances.values()),minFrequency)
                                    boundary = math.floor(math.sqrt(maxDistanceFrequency/minFrequency))
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
                                                    if number>=minFrequency:
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
                                    maxDistanceFrequency = max(max(distances.values()),maxFrequency)
                                    boundary = math.floor(math.sqrt(maxDistanceFrequency/minFrequency))
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
                                    maxDistanceFrequency = max(distances.values())
                                    boundary = math.floor(math.sqrt(maxDistanceFrequency/minFrequency))
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
                                                        number>=minFrequency):
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
                                        number = item[2][k][3]
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
                                        number = storageHandler["otherRow"]["number"]
                                        distance = storageHandler["otherRow"]["distance"]
                                        directData[index] = add_directData(directData[index],
                                                                       fromDirection,toLink,toDirection,distance,number)
                                        for row in storageHandler["otherIterator"]:
                                            if row["fromLink"]==fromLink:
                                                storageHandler["otherRow"] = row
                                                fromDirection = storageHandler["otherRow"]["fromDirection"]
                                                toDirection = storageHandler["otherRow"]["toDirection"]
                                                toLink = storageHandler["otherRow"]["toLink"]
                                                number = storageHandler["otherRow"]["number"]
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
                                        cycleEntries[index] = (cycleEntries[index][0]+item[0][0],
                                                               min(cycleEntries[index][1],item[0][1]),)
                                    else:
                                        cycleEntries[index] = (item[0][0],item[0][1],)
                            for index, item in enumerate(rowData):
                                if item[1][1]>0:
                                    if index in reversalEntries.keys():
                                        reversalEntries[index] = (reversalEntries[index][0]+item[1][0],
                                                               min(reversalEntries[index][1],item[1][1]),)
                                    else:
                                        reversalEntries[index] = (item[1][0],item[1][1],)
                        for ce in cycleEntries.keys():
                            minimumLength = 0
                            tableCycleRow = tableCycle.row
                            tableCycleRow["ckmerLink"] = i+ce
                            tableCycleRow["number"] = cycleEntries[ce][0]
                            tableCycleRow["minimumLength"] = cycleEntries[ce][1]
                            tableCycleRow.append()
                        for ce in reversalEntries.keys():
                            minimumLength = 0
                            tableReversalRow = tableReversal.row
                            tableReversalRow["ckmerLink"] = i+ce
                            tableReversalRow["number"] = reversalEntries[ce][0]
                            tableReversalRow["minimumLength"] = reversalEntries[ce][1]
                            tableReversalRow.append()
                        
                        #TODO: apply some filtering, combining and detection of problematic entries
                        # problematic if
                        # - different distances in single direction
                        # - too much connections
                        # then try to filter by (using data available in kmer_properties array from shared memory)
                        # - check shared bases for each direction
                        # - check split direction
                        # otherwise classify problem
                        # Filtering has to be symmetric!

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
                                            tableDirectRow["number"] = number
                                            tableDirectRow["distance"] = distance
                                            tableDirectRow["splitDirection"] = direction
                                            tableDirectRow["reverseBase"] = reverseBase
                                            tableDirectRow["forwardBase"] = forwardBase
                                            tableDirectRow["problematic"] = multipleDistances + (problematic<<1)
                                            tableDirectRow.append()
                                            
                        
                
                #finished
                pytables_storage.flush()
                tableDeleteDirect.cols.fromLink.create_csindex()
                pytables_storage.flush()
                return pytablesFile
                
                
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("merges: shared memory of {} bytes used".format(shm.size))
        try:
            #get shared memory
            shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
            shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maxFrequency)).type
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
                        pytablesFile = merge_storage(filenameBase, storageFiles, mergeStart, 
                                                     mergeNumber, numberOfKmers, arrayNumber, kmer_properties)
                        #now the file can be released for final merge
                        queue_merges.put(pytablesFile) 
                except Empty:
                    continue
        finally:
            shm.close()
        logger.debug("index: shared memory released")
            

    
    def combine_merges(mergeFiles,filenameBase,numberOfKmers):
        pytablesFile = (filenameBase+"_tmp_direct_merge.h5")
        if os.path.exists(pytablesFile):
            os.remove(pytablesFile)
            
        with tables.open_file(pytablesFile, mode="w") as pytables_storage:
            
            nCycle = 0
            nReversal = 0
            nDirect = 0
            nDeleteDirect = 0
            
            #get dimensions
            sortedMergeFiles = [{"filename": mergeFile} for mergeFile in sorted(mergeFiles)]
            for i in range(len(sortedMergeFiles)):
                with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                    sortedMergeFiles[i]["ncycle"] = pytables_merge.root.cycle.shape[0]
                    sortedMergeFiles[i]["nreversal"] = pytables_merge.root.reversal.shape[0]
                    sortedMergeFiles[i]["ndirect"] = pytables_merge.root.direct.shape[0]
                    sortedMergeFiles[i]["ndeleteDirect"] = pytables_merge.root.deleteDirect.shape[0]
                    nCycle += sortedMergeFiles[i]["ncycle"]
                    nReversal += sortedMergeFiles[i]["nreversal"]
                    nDirect += sortedMergeFiles[i]["ndirect"]
                    nDeleteDirect += sortedMergeFiles[i]["ndeleteDirect"]
                    
            #create temporary storage with dimensions
            Storage.create_merge_storage(pytables_storage, numberOfKmers, nCycle, nReversal, nDirect, False)
            tableCycle = pytables_storage.root.cycle
            tableReversal = pytables_storage.root.reversal
            tableDirect = pytables_storage.root.direct
                    
            stepSizeStorage = 100000
            maxCycleLength = 0
            maxCycleNumber = 0
            maxReversalLength = 0
            maxReversalNumber = 0
            maxDirectDistance = 0
            maxDirectNumber = 0
            for i in range(len(sortedMergeFiles)):
                with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                    for j in range(0,sortedMergeFiles[i]["ncycle"],stepSizeStorage):
                        reducedCycleDataBlock = pytables_merge.root.cycle[j:j+stepSizeStorage]
                        maxCycleNumber = max(maxCycleNumber,max(reducedCycleDataBlock["number"]))
                        maxCycleLength = max(maxCycleLength,max(reducedCycleDataBlock["minimumLength"]))
                        tableCycle.append(reducedCycleDataBlock)
                    for j in range(0,sortedMergeFiles[i]["nreversal"],stepSizeStorage):
                        reducedReversalDataBlock = pytables_merge.root.reversal[j:j+stepSizeStorage]
                        maxReversalNumber = max(maxReversalNumber,max(reducedReversalDataBlock["number"]))
                        maxReversalLength = max(maxReversalLength,max(reducedReversalDataBlock["minimumLength"]))
                        tableReversal.append(reducedReversalDataBlock)
                    for j in range(0,sortedMergeFiles[i]["ndirect"],stepSizeStorage):
                        directDataBlock = pytables_merge.root.direct[j:j+stepSizeStorage]
                        deleteDataBlock = pytables_merge.root.deleteDirect.read_where(
                            "(fromLink>={}) & (fromLink<={})".format(
                            directDataBlock[0]["fromLink"],directDataBlock[-1]["fromLink"]))
                        deletableIndices = np.where(
                            np.in1d(directDataBlock[["fromLink","fromDirection","toLink","toDirection"]],
                                    deleteDataBlock[["fromLink","fromDirection","toLink","toDirection"]]))[0]
                        reducedDirectDataBlock = np.delete(directDataBlock,deletableIndices,0)
                        maxDirectNumber = max(maxDirectNumber,max(reducedDirectDataBlock["number"]))
                        maxDirectDistance = max(maxDirectDistance,max(reducedDirectDataBlock["distance"]))
                        #add reduced direct data
                        tableDirect.append(reducedDirectDataBlock)
            tableCycle.attrs.maxNumber = maxCycleNumber
            tableCycle.attrs.maxLength = maxCycleLength
            tableReversal.attrs.maxNumber = maxReversalNumber
            tableReversal.attrs.maxLength = maxReversalLength
            tableDirect.attrs.maxDistance = maxDirectDistance
            tableDirect.attrs.maxNumber = maxDirectNumber
            pytables_storage.flush()
            
    def store_merged_direct(h5file,filenameBase,numberOfKmers,minFrequency):
        pytablesFile = (filenameBase+"_tmp_direct_merge.h5")
        logger = logging.getLogger(__name__)                
        
        def classifyDirectProblematic(dataBlock,directGood,directCoverage,directProblematic):
            dataBlock = np.array(dataBlock, dtype=object)                       
            if max(dataBlock[:,4])==0:
                directGood+=len(dataBlock)
            elif len(dataBlock)==1:
                dataBlock[:,4] = 0
                directGood+=len(dataBlock)
            elif max(dataBlock[:,2])<minFrequency:
                dataBlock[:,4] = 1
                directCoverage+=len(dataBlock)
            else:
                dataBlock[:,4] = 2
                directProblematic+=len(dataBlock)
            number = sum(dataBlock[:,2])
            distinct = len(dataBlock)
            dataBlock = tuple([tuple(e) for e in dataBlock])
            return (dataBlock,directGood,directCoverage,directProblematic,distinct,number)
        
        if not os.path.exists(pytablesFile):
            logger.error("couldn't find {}".format(pytablesFile))
        else:
            stepSizeStorage = 100000
            dsCkmer = h5file["/split/ckmer"]
            with tables.open_file(pytablesFile, mode="r") as pytables:
                #cycles
                numberOfCycles = pytables.root.cycle.shape[0]
                maxNumber = pytables.root.cycle.attrs.maxNumber
                maxLength = pytables.root.cycle.attrs.maxLength
                dtypeCycleList=[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                                 ("minimumLength",haplotyping.index.Database.getUint(maxLength)),
                                 ("number",haplotyping.index.Database.getUint(maxNumber))]
                dtCycle=np.dtype(dtypeCycleList)
                dsCycle=h5file["/relations/"].create_dataset("cycle",(numberOfCycles,), 
                                                             dtype=dtCycle, chunks=None)
                maxCycleLength = 0
                maxCycleNumber = 0
                for i in range(0,numberOfCycles,stepSizeStorage):
                    stepData = pytables.root.cycle[i:i+stepSizeStorage]                    
                    dsCycle[i:i+stepSizeStorage] = stepData
                    for row in stepData:
                        ckmerRow = dsCkmer[row[0]]
                        ckmerRow[6] = row[2]
                        maxCycleLength = max(maxCycleLength,row[1])
                        maxCycleNumber = max(maxCycleNumber,row[2])
                        dsCkmer[row[0]] = ckmerRow    
                h5file["/config/"].attrs["numberOfCycles"]=numberOfCycles        
                h5file["/config/"].attrs["maxCycleLength"]=maxCycleLength
                h5file["/config/"].attrs["maxCycleNumber"]=maxCycleNumber
                logger.info("store "+str(numberOfCycles)+" cycles")
                #reversals
                numberOfReversals = pytables.root.reversal.shape[0]
                maxNumber = pytables.root.reversal.attrs.maxNumber
                maxLength = pytables.root.reversal.attrs.maxLength
                dtypeReversalList=[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                                 ("minimumLength",haplotyping.index.Database.getUint(maxLength)),
                                 ("number",haplotyping.index.Database.getUint(maxNumber))]
                dtReversal=np.dtype(dtypeReversalList)
                dsReversal=h5file["/relations/"].create_dataset("reversal",(numberOfReversals,), 
                                                                dtype=dtReversal, chunks=None)
                maxReversalLength = 0
                maxReversalNumber = 0
                for i in range(0,numberOfReversals,stepSizeStorage):
                    stepData = pytables.root.reversal[i:i+stepSizeStorage]
                    dsReversal[i:i+stepSizeStorage] = stepData
                    for row in stepData:
                        ckmerRow = dsCkmer[row[0]]
                        ckmerRow[7] = row[2]
                        maxReversalLength = max(maxReversalLength,row[1])
                        maxReversalNumber = max(maxReversalNumber,row[2])
                        dsCkmer[row[0]] = ckmerRow     
                h5file["/config/"].attrs["numberOfReversals"]=numberOfReversals
                h5file["/config/"].attrs["maxReversalLength"]=maxReversalLength
                h5file["/config/"].attrs["maxReversalNumber"]=maxReversalNumber
                logger.info("store "+str(numberOfReversals)+" reversals")
                #direct relations
                numberOfDirectRelations = pytables.root.direct.shape[0]
                maxNumber = pytables.root.direct.attrs.maxNumber
                maxDistance = pytables.root.direct.attrs.maxDistance
                dtypeDirectList=[("from",[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                                          ("direction","S1")]),
                                 ("to",[("ckmerLink",haplotyping.index.Database.getUint(numberOfKmers)),
                                          ("direction","S1")]),
                                 ("number",haplotyping.index.Database.getUint(maxNumber)),
                                 ("distance",haplotyping.index.Database.getUint(maxDistance)),
                                 ("problematic","uint8")]
                dtDirect=np.dtype(dtypeDirectList)
                dsDirect=h5file["/relations/"].create_dataset("direct",(numberOfDirectRelations,), 
                                                              dtype=dtDirect, chunks=None)
                directCounter = 0
                previousStepData = []
                directGood = 0
                directCoverage = 0
                directProblematic = 0
                for i in range(0,numberOfDirectRelations,stepSizeStorage):
                    stepData = pytables.root.direct[i:i+stepSizeStorage]
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

            