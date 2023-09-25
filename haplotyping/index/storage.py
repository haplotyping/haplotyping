import haplotyping.index.database
from multiprocessing import shared_memory, current_process
from statistics import multimode
from queue import Empty
import ahocorasick, metis, networkit as nk
import os, re, pickle, tables, statistics, logging, time, psutil
import numpy as np, math
from contextlib import ExitStack

class Storage:
    
    """
    Internal use, storage and processing
    """
    
    stepSizeStorage = 1000000
    
    def create_mergeDirect_storage(pytablesStorage, numberOfKmers, maximumFrequency,
                             nCycle=None, nReversal=None, nDirect=None, nPaired=None):
        pytablesStorage.create_table(pytablesStorage.root, 
            "tmpPaired",{
            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
        }, "Temporary Paired", expectedrows=nPaired)
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
        
    def create_filteredReads_storage(pytablesStorage,numberOfReads,numberOfData,numberOfKmers,
                                  numberOfPartitions,maximumReadLength,numberOfPaired):
        filters = tables.Filters(complevel=9, complib="blosc")
        pytablesStorage.create_earray(
                      pytablesStorage.root, "readPartitionData",
                      haplotyping.index.Database.getTablesUintAtom(numberOfKmers), 
                      shape=(0,), filters=filters, 
                      expectedrows=numberOfData)
        pytablesStorage.create_table(
                      pytablesStorage.root, "readPartitionInfo", {
                        "length": haplotyping.index.Database.getTablesUint(maximumReadLength,0),
                        "partition": haplotyping.index.Database.getTablesUint(numberOfPartitions,1)
                      }, "Read size and partition",
                      expectedrows=numberOfReads)
        pytablesStorage.create_table(
                      pytablesStorage.root, "readPaired", {
                        "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                        "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
                      }, "Paired Reads",
                      expectedrows=numberOfPaired)
        
    def create_filteredMergedReads_storage(pytablesStorage,numberOfReads,numberOfData,numberOfKmers,
                                            numberOfPartitions,maximumReadLength):
        filters = tables.Filters(complevel=9, complib="blosc")
        pytablesStorage.create_earray(
                      pytablesStorage.root, "readPartitionData",
                      haplotyping.index.Database.getTablesUintAtom(numberOfKmers), 
                      shape=(0,), filters=filters, 
                      expectedrows=numberOfData)
        pytablesStorage.create_table(
                      pytablesStorage.root, "readPartitionInfo", {
                        "length": haplotyping.index.Database.getTablesUint(maximumReadLength,0),
                        "number": haplotyping.index.Database.getTablesUint(numberOfReads,1)
                      }, "Read size and number",
                      expectedrows=numberOfReads)
        pytablesStorage.create_table(
                      pytablesStorage.root, "readPartition", {
                        "readPartitionDataLink": haplotyping.index.Database.getTablesUint(numberOfData,0),
                        "readPartitionDataNumber": haplotyping.index.Database.getTablesUint(numberOfData,1),
                        "readPartitionInfoLink": haplotyping.index.Database.getTablesUint(numberOfReads,2),
                        "readPartitionInfoNumber": haplotyping.index.Database.getTablesUint(numberOfReads,3)
                      }, "Index data and info readPartitions",
                      expectedrows=numberOfPartitions)
        
    def create_mergeReads_storage(pytablesStorage,numberOfReads,numberOfData,numberOfKmers,
                                  numberOfPartitions,maximumReadLength):
        filters = tables.Filters(complevel=9, complib="blosc")
        pytablesStorage.create_carray(
                      pytablesStorage.root, "readPartitionData",
                      haplotyping.index.Database.getTablesUintAtom(numberOfKmers), 
                      shape=(numberOfData,), filters=filters)
        pytablesStorage.create_carray(
                      pytablesStorage.root, "readPartitionLength",
                      haplotyping.index.Database.getTablesUintAtom(maximumReadLength), 
                      shape=(numberOfReads,), filters=filters)
        pytablesStorage.create_carray(
                      pytablesStorage.root, "readPartitionInfo",
                      haplotyping.index.Database.getTablesUintAtom(numberOfPartitions), 
                      shape=(numberOfReads,), filters=filters)        
    
    def workerAutomaton(shutdown_event,queue_start,queue_automaton,queue_index,
                         queue_finished,k,automatonKmerSize,automatonFile):
        
        logger = logging.getLogger("{}.worker.automaton".format(__name__))
                
        def compute_matches(sequence,automatonSplits):
            rsequence = haplotyping.General.reverse_complement(sequence)
            boundary = len(sequence)-k+automatonKmerSize
            correction_reverse = boundary - 1
            correction_forward = automatonKmerSize - 1
            #use automaton to check reverse complement sequence
            rdict = {(correction_reverse - end_index): (number, startLinks,)
                     for (end_index, (number,startLinks)) in automatonSplits.iter(rsequence,0,boundary)}
            if len(rdict)>0:
                flist = [(end_index - correction_forward, number, startLinks,)
                         for (end_index, (number,startLinks)) in automatonSplits.iter(sequence,0,boundary)]
                #combine results
                clist = [(pos, (number, startLinks,) ,rdict[pos]) 
                         for (pos, number, startLinks) in flist 
                         if pos in rdict.keys()]
            else:
                clist = []
            return clist
        
        try:
            
            #wait for permission to start loading automaton
            while True:
                try:
                    item = queue_start.get(block=True, timeout=1)
                    if item=="automaton":
                        break
                    else:
                        logger.error("automaton ({}): unexpected start value in queue: {}".format(os.getpid(),item))
                except:
                    pass
                time.sleep(1)
            
            logger.debug("automaton ({}): load automaton from {}".format(os.getpid(),automatonFile))
            
            #load automaton into memory
            automatonSplits = ahocorasick.load(automatonFile,pickle.loads)
            logger.debug("automaton ({}): fsm loaded from {} MB".format(
                os.getpid(),math.ceil(os.stat(automatonFile).st_size/1048576)))

            process = psutil.Process(os.getpid())
            logger.debug("automaton ({}): used memory {} MB".format(
                os.getpid(),math.ceil(process.memory_info().rss/1048576)))
            
            #finished loading automaton
            queue_finished.put("automaton:started")

            while not shutdown_event.is_set():
                try:
                    item = queue_automaton.get(block=True, timeout=1)
                    if item==None:
                        logger.debug("autmaton ({}): none item".format(os.getpid()))
                        break
                    elif isinstance(item,tuple) and len(item)==2:
                        queue_index.put((
                            (item[0],compute_matches(item[0],automatonSplits),),
                            (item[1],compute_matches(item[1],automatonSplits),),
                        ))
                    elif isinstance(item,str):
                        queue_index.put((
                            (item,compute_matches(item,automatonSplits),),
                        ))
                except Empty:
                    logger.debug("automaton ({}): empty".format(os.getpid()))
                    time.sleep(5)
                    continue
        except Exception as ex:
            logger.error("automaton ({}): problem with worker: {}".format(os.getpid(),ex))
        finally:
            del automatonSplits
            logger.debug("automaton ({}): fsm released".format(os.getpid()))
        queue_finished.put("automaton:ended")
            
                
    
    def workerIndex(shutdown_event,queue_index,queue_matches,queue_storage,queue_finished,
                     filenameBase,numberOfKmers,k,indexType,shm_name):

        logger = logging.getLogger("{}.worker.index".format(__name__))
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
            totalChecks=len(clist)
            totalMatches = 0
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
                                nkmer = shm.buf[indexLink:indexLink+k].tobytes().decode()
                                if rkmer==nkmer:
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
                                nkmer = shm.buf[indexLink:indexLink+k].tobytes().decode()
                                if kmer==nkmer:
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
                                nkmer = shm.buf[indexLink:indexLink+k].tobytes().decode()
                                if kmer==nkmer:
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
                                nkmer = shm.buf[indexLink:indexLink+k].tobytes().decode()
                                if rkmer==nkmer:
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
                    totalMatches+=len(matches)
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
                totalMatches+=len(matches)
            return (matchesList, len(matchesList)==1, totalChecks, totalMatches)
        
        def normalised_matchesList(matches):
            matchesList = []
            for item in matches:
                matchesList.extend([s[1] for s in item])
            if len(matchesList)==0:
                return matchesList
            elif matchesList[0]<matchesList[-1]:
                return matchesList
            elif matchesList[0]>matchesList[-1]:
                return matchesList[::-1]
            elif hash(tuple(matchesList))<hash(tuple(matchesList[::-1])):
                return matchesList
            else:
                return matchesList[::-1]
            
        def store_matches(matches,readData,readInfo):
            matchesList = normalised_matchesList(matches)
            if len(matchesList)>0:
                readData.append(tuple(matchesList))
                readInfo.append([(len(matchesList),0,)])
            
        def store_paired_matches(matches0,matches1,readData,readInfo):
            matches0List = normalised_matchesList(matches0)
            matches1List = normalised_matchesList(matches1)
            if len(matches0List)==0 and len(matches1List)==0:
                pass
            elif len(matches0List)==0:
                store_matches(matches1,readData,readInfo)
            elif len(matches1List)==0:
                store_matches(matches0,readData,readInfo)
            else:
                readData.append(tuple(matches0List))
                readData.append(tuple(matches1List))
                readInfo.append([(len(matches0List),1,)])
                readInfo.append([(len(matches1List),2,)])
            
                    
        #stats
        totalChecks = 0
        totalMatches = 0
        
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("index ({}): shared memory of {} MB used".format(
            os.getpid(),math.ceil(shm.size/1048576)))
        
        process = psutil.Process(os.getpid())
        logger.debug("index ({}): used memory {} MB".format(
            os.getpid(),math.ceil(process.memory_info().rss/1048576)))
        
        try:
            curr_proc = current_process()
            pytablesFileWorker = filenameBase+"_tmp_reads_{}.process.h5".format(curr_proc.name)
            if os.path.exists(pytablesFileWorker):
                os.remove(pytablesFileWorker)

            with (open(os.devnull,"w") if indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS
                      else tables.open_file(pytablesFileWorker, mode="w")) as pytablesStorageWorker:
                
                if not indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS:
                    filters = tables.Filters(complevel=9, complib="blosc")
                    readData = pytablesStorageWorker.create_earray(pytablesStorageWorker.root, "readRawData",
                                      haplotyping.index.Database.getTablesUintAtom(numberOfKmers), 
                                      shape=(0,), filters=filters)
                    readInfo = pytablesStorageWorker.create_table(
                                      pytablesStorageWorker.root, "readRawInfo", {
                                        "length": tables.UInt32Col(pos=0),
                                        "paired": tables.UInt8Col(pos=1)
                                      }, "Read size and paired")
                    
                while not shutdown_event.is_set():
                    try:
                        item = queue_index.get(block=True, timeout=1)
                        if item==None:
                            logger.debug("index ({}): none item".format(os.getpid()))
                            break
                        elif isinstance(item,tuple):
                            if len(item)==1:
                                (matches,direct,tmpTotalChecks,tmpTotalMatches,) = compute_matches(item[0][0],item[0][1])
                                totalChecks+=tmpTotalChecks
                                totalMatches+=tmpTotalMatches
                                if tmpTotalMatches<=1:
                                    pass
                                else:
                                    queue_matches.put(((matches, direct,),))
                                    if (not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS)
                                        and tmpTotalMatches>2):
                                        store_matches(matches,readData,readInfo)
                            elif len(item)==2:
                                (matches0,direct0,tmpTotalChecks0,tmpTotalMatches0,) = compute_matches(item[0][0],item[0][1])
                                (matches1,direct1,tmpTotalChecks1,tmpTotalMatches1,) = compute_matches(item[1][0],item[1][1])
                                totalChecks+=tmpTotalChecks0+tmpTotalChecks1
                                totalMatches+=tmpTotalMatches0+tmpTotalMatches1
                                if tmpTotalMatches0==0 and tmpTotalMatches1==0:
                                    pass
                                elif tmpTotalMatches0==0:
                                    if tmpTotalMatches1>1:
                                        queue_matches.put(((matches1, direct1,),))
                                        if (not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS) 
                                            and tmpTotalMatches1>2):
                                            store_matches(matches1,readData,readInfo)
                                elif tmpTotalMatches1==0:
                                    if tmpTotalMatches0>1:
                                        queue_matches.put(((matches0, direct0,),))
                                        if (not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS)
                                            and tmpTotalMatches0>2):
                                            store_matches(matches0,readData,readInfo)
                                else:
                                    queue_matches.put(((matches0, direct0, ),
                                                       (matches1, direct1, )))                                    
                                    if (not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS)):
                                        if (tmpTotalMatches0>1) and (tmpTotalMatches1>1):
                                            store_paired_matches(matches0,matches1,readData,readInfo)
                                        else:
                                            if tmpTotalMatches0>1:
                                                store_matches(matches0,readData,readInfo)
                                            if tmpTotalMatches1>1:
                                                store_matches(matches1,readData,readInfo)
                    except Empty:
                        logger.debug("index ({}): empty".format(os.getpid()))
                        time.sleep(5)
                        continue
            #now the file can be released for later processing
            if not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS):
                queue_storage.put(pytablesFileWorker) 
        except Exception as ex:
           logger.error("index  ({}): problem with worker: {}".format(os.getpid(),ex))
        finally:
            shm.close()
        logger.debug("index ({}): found {} matches in {} checks".format(os.getpid(),totalMatches,totalChecks))
        queue_finished.put("index:ended:{}:{}".format(totalChecks,totalMatches))
            
    
    """
    Process matches into direct connections and queue full list of matches for indirect connections
    """
    def worker_matches_dtype(numberOfKmers,maximumFrequency,estimatedMaximumReadLength,numberDirectArray):
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
                                     for i in range(numberDirectArray)]),]
        return dtype
    
    def workerMatches(shutdown_event,queue_matches,queue_storage,queue_finished,
                       filenameBase,numberOfKmers,maximumFrequency,estimatedMaximumReadLength,
                       numberDirectArray,indexType,shm_name):
        
        logger = logging.getLogger("{}.worker.matches".format(__name__))
        try:
            curr_proc = current_process()
            pytablesFileWorker = filenameBase+"_tmp_direct_{}.process.h5".format(curr_proc.name)
            if os.path.exists(pytablesFileWorker):
                os.remove(pytablesFileWorker)

            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                
                logger.debug("matches ({}): store direct connections (dimension: {})".format(
                    os.getpid(),numberDirectArray))
                
                #stats
                totalChecks = 0
                totalDirect = 0
                totalCycle = 0
                totalReversal = 0
                
                #define correct maxValues based on previous results             
                dtype = haplotyping.index.storage.Storage.worker_matches_dtype(
                    numberOfKmers,maximumFrequency,estimatedMaximumReadLength,numberDirectArray)
                connections = np.ndarray((numberOfKmers,), dtype=dtype, order="C")
                connections.fill(((0,0,),(0,0,),tuple((0,0,0,0,) for i in range(numberDirectArray))))
                
                logger.debug("matches ({}): created memory storage for direct connections: {} MB".format(
                    os.getpid(), math.ceil(connections.nbytes/1048576)))
                
                #store paired connections for partition graph
                tablePaired = pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                    "tmpPaired",{
                    "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
                    "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
                }, "Temporary to dump paired relations", track_times=False)
                
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

                def store_paired(linkFrom, linkTo):
                    if not linkFrom==linkTo:
                        tablePaired.append([(linkFrom,linkTo,)])

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
                    numberFilled = numberDirectArray
                    connectionType = 1 + fromDirection + (2*toDirection)
                    directRow = connections[fromLink][2]
                    for i in range(numberDirectArray):
                        if directRow[i][0]==0:
                            numberFilled = i
                            break
                        elif (directRow[i][1]==toLink and directRow[i][0]==connectionType 
                              and directRow[i][2]==distance):
                            directRow[i][3]+=1 
                            return
                    if numberFilled<numberDirectArray:
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
                        if 2*len(sameFromDirectionDistances)>=numberDirectArray:
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
                        elif 2*len(otherFromDirectionDistances)>numberDirectArray:
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

                def process_matches(matchesList,totalDirect,totalReversal,totalCycle,totalChecks):
                    history = {}
                    links = []
                    startPos = 0
                    endPos = 0
                    for matches in matchesList:
                        totalChecks+=len(matches)
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
                                totalDirect+=1
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
                                        totalCycle+=1
                                        store_cycle(currentLink,currentPos-history[currentLink][1])
                                    else:
                                        totalReversal+=1
                                        store_reversal(currentLink,currentPos-history[currentLink][1])
                                history[previousLink]=(previousDirection,previousPos,)
                    return (links,max(0,endPos-startPos),totalDirect,totalReversal,totalCycle,totalChecks,)

                process = psutil.Process(os.getpid())
                logger.debug("matches ({}): used memory {} MB".format(
                    os.getpid(),math.ceil(process.memory_info().rss/1048576)))
                
                shm = shared_memory.SharedMemory(shm_name)
                logger.debug("matches ({}): shared memory of {} MB used".format(os.getpid(),math.ceil(shm.size/1048576)))
                
                #get shared memory
                shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
                shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maximumFrequency)).type
                kmer_properties = np.ndarray((numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                                       ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm.buf)
        
                while not shutdown_event.is_set():
                    try:
                        item = queue_matches.get(block=True, timeout=1)
                        if item==None:
                            logger.debug("matches ({}): none item".format(os.getpid()))
                            break
                        elif isinstance(item,tuple):
                            if len(item)==1:
                                matchesList = item[0][0]
                                directConnected = item[0][1]
                                (links,length,totalDirect,totalReversal,totalCycle,totalChecks,) = process_matches(
                                    matchesList,totalDirect,totalReversal,totalCycle,totalChecks)
                            elif len(item)==2:
                                matchesList0 = item[0][0]
                                directConnected0 = item[0][1]
                                matchesList1 = item[1][0]
                                directConnected1 = item[1][1]
                                (links0,length0,totalDirect,totalReversal,totalCycle,totalChecks,) = process_matches(
                                    matchesList0,totalDirect,totalReversal,totalCycle,totalChecks)
                                (links1,length1,totalDirect,totalReversal,totalCycle,totalChecks,) = process_matches(
                                    matchesList1,totalDirect,totalReversal,totalCycle,totalChecks)
                                #register pair data for single matches
                                if (len(links0)==1) and (len(links1)>0):
                                    pairFrom = links0[0][0]
                                    pairTo = None
                                    minFreq = 0
                                    for link in links1:
                                        freq = kmer_properties[link[0]][1]
                                        #not informative
                                        if link[0]==pairFrom:
                                            pairFrom = None
                                            break
                                        #connect to first k-mer with minimum frequency
                                        if freq==minFreq:                                            
                                            pairTo = link[0] if (pairTo==None) else min(pairTo,link[0])
                                        elif minFreq==0 or freq<minFreq:
                                            minFreq = freq
                                            pairTo = link[0]
                                    if not pairFrom is None:
                                        if not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS):
                                            store_paired(pairFrom,pairTo)
                                if (len(links1)==1) and (len(links0)>0):
                                    pairFrom = links1[0][0]
                                    pairTo = None
                                    minFreq = 0
                                    for link in links0:
                                        freq = kmer_properties[link[0]][1]
                                        #not informative
                                        if link[0]==pairTo:
                                            pairTo = None
                                            break
                                        #connect to first k-mer with minimum frequency
                                        if freq==minFreq:
                                            pairTo = link[0] if (pairTo==None) else min(pairTo,link[0])
                                        elif minFreq==0 or freq<minFreq:
                                            minFreq = freq
                                            pairTo = link[0]
                                    if not pairTo is None:
                                        if not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS):
                                            store_paired(pairFrom,pairTo)
                    except Empty:
                        logger.debug("matches ({}): empty".format(os.getpid()))
                        time.sleep(5)
                        continue   
                pytablesStorageWorker.flush()
                logger.debug("matches ({}): create indices temporary tables".format(os.getpid()))
                tableDirectOther.cols.fromLink.create_csindex()
                pytablesStorageWorker.flush()
                if not (indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS):
                    tablePaired.cols.fromLink.create_csindex()
                    pytablesStorageWorker.flush()
                pytablesStorageWorker.create_table(pytablesStorageWorker.root, 
                                              name="direct", obj=connections, expectedrows=numberOfKmers)
                logger.debug("matches ({}): found {} direct, {} cycle and {} reversal in {} matches".format(
                    os.getpid(), totalDirect,totalCycle,totalReversal,totalChecks))
                pytablesStorageWorker.flush()
                
            #now the file can be released for merge
            queue_storage.put(pytablesFileWorker) 
        except Exception as ex:
            logger.error("matches ({}): problem with worker: {}".format(os.getpid(), ex))
        finally:
            shm.close()
        queue_finished.put("matches:ended")
            
                
    """
    Merge stored direct and indirect connections
    """                
    def workerMergeDirect(shutdown_event,queue_ranges,queue_merges,storageDirectFiles,
                      filenameBase,numberOfKmers,maximumFrequency,minimumFrequency,shm_name):
        
        logger = logging.getLogger("{}.worker.merges".format(__name__))

        def merge_direct_storage(pytablesFileWorker, storageDirectFiles, mergeStart, mergeNumber, 
                                 numberOfKmers, kmer_properties):

            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                tableCycle = pytablesStorageWorker.root.cycle
                tableReversal = pytablesStorageWorker.root.reversal
                tableDirect = pytablesStorageWorker.root.direct
                tableDeleteDirect = pytablesStorageWorker.root.deleteDirect                
                tablePaired = pytablesStorageWorker.root.tmpPaired
                
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
                    #get and position sorted iterators for paired entries
                    for i in range(len(storageHandlers)):
                        storageHandlers[i]["pairedRow"] = None
                        if storageHandlers[i]["handler"].root.tmpPaired.shape[0]>0:
                            storageHandlers[i]["pairedIterator"] = storageHandlers[i]["handler"].root.tmpPaired.itersorted(
                                "fromLink",checkCSI=True)
                            for row in storageHandlers[i]["pairedIterator"]:
                                if row["fromLink"]>=mergeStart:
                                    if row["fromLink"]<=mergeEnd:
                                        storageHandlers[i]["pairedRow"] = row
                                    break
                        else:
                            storageHandlers[i]["pairedIterator"] = None
                    
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
                            #first
                            tableDeleteDirectRow = tableDeleteDirect.row
                            tableDeleteDirectRow["fromLink"] = fromLink
                            tableDeleteDirectRow["fromDirection"] = fromDirection
                            tableDeleteDirectRow["toLink"] = toLink
                            tableDeleteDirectRow["toDirection"] = toDirection
                            tableDeleteDirectRow.append()
                            #second
                            tableDeleteDirectRow = tableDeleteDirect.row
                            tableDeleteDirectRow["fromLink"] = toLink
                            tableDeleteDirectRow["fromDirection"] = toDirection
                            tableDeleteDirectRow["toLink"] = fromLink
                            tableDeleteDirectRow["toDirection"] = fromDirection
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
                    
                    for i in range(mergeStart,mergeEnd+1,Storage.stepSizeStorage):
                        #initialise
                        cycleEntries = {}
                        reversalEntries = {}
                        directData = [{} for j in range(min(mergeEnd+1-i,Storage.stepSizeStorage,(numberOfKmers-i)))]
                        for storageHandler in storageHandlers:
                            rowData = storageHandler["handler"].root.direct[i:min(mergeEnd+1,i+Storage.stepSizeStorage)]
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
                        
                        unlinkedEntries = []
                        for j in range(len(directData)):
                            fromLink = i + j
                            if len(directData[j])==0:
                                unlinkedEntries.append(fromLink)
                            else:
                                for fromDirection in directData[j].keys():
                                    directData[j][fromDirection] = reduceDistances(fromLink,fromDirection, 
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

                        #store paired data for not directly connected entries
                        if (len(unlinkedEntries)>0):
                            for storageHandler in storageHandlers:
                                if not storageHandler["pairedRow"]==None:
                                    try:
                                        unlinkedIterator = iter(unlinkedEntries)
                                        paired = storageHandler["pairedRow"]
                                        unlinked = next(unlinkedIterator)
                                        while True:
                                            if paired[0] > unlinked:
                                                unlinked = next(unlinkedIterator)
                                            else:
                                                if paired[0] == unlinked:     
                                                    tablePaired.append([(paired[0],paired[1],)])
                                                paired = next(storageHandler["pairedIterator"])
                                    except StopIteration:
                                        pass                           

                #finished
                pytablesStorageWorker.flush()
                
                        
                
        shm = shared_memory.SharedMemory(shm_name)
        logger.debug("merges ({}): shared memory of {} MB used".format(os.getpid(),math.ceil(shm.size/1048576)))
        try:
            #get shared memory
            shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
            shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maximumFrequency)).type
            kmer_properties = np.ndarray((numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                               ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm.buf)
            
            process = psutil.Process(os.getpid())
            logger.debug("merges ({}): used memory {} MB".format(
                os.getpid(),math.ceil(process.memory_info().rss/1048576)))
            
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
                        pytablesFileRange = (filenameBase+"_tmp_direct_merge_"+
                                             str(mergeStart).zfill(numberLength)+"_"+
                                        str(mergeEnd).zfill(numberLength)+".process.h5")
                        if os.path.exists(pytablesFileRange):
                            os.remove(pytablesFileRange)
                        with tables.open_file(pytablesFileRange, mode="a") as pytablesStorageRange:
                            Storage.create_mergeDirect_storage(pytablesStorageRange, numberOfKmers, 
                                                         maximumFrequency)
                        #merge
                        merge_direct_storage(pytablesFileRange, storageDirectFiles, mergeStart, 
                                             mergeNumber, numberOfKmers, kmer_properties)
                        #now the file can be released for final merge
                        queue_merges.put(pytablesFileRange) 
                except Empty:
                    time.sleep(1)
                    continue
        except Exception as ex:
            logger.error("merges ({}): problem with worker: {}".format(os.getpid(),ex))
        finally:
            shm.close()
            
            
    def workerProcessReads(shutdown_event,queue_rawReads,queue_filteredReads,queue_finished,filenameBase,numberOfKmers,
                     numberOfPartitions,numberOfDirect,maximumFrequency,maximumReadLength,shm_kmer_name,shm_direct_name):

        logger = logging.getLogger("{}.worker.index".format(__name__))
        
        try:
            shm_kmer = shared_memory.SharedMemory(shm_kmer_name)
            shm_kmer_partition = np.dtype(haplotyping.index.Database.getUint(numberOfPartitions)).type
            shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(maximumFrequency)).type
            shm_kmer_reference = np.dtype(haplotyping.index.Database.getUint(numberOfDirect)).type
            shm_kmer_numberLeft = np.dtype("uint8").type
            shm_kmer_numberRight = np.dtype("uint8").type

            shm_direct = shared_memory.SharedMemory(shm_direct_name)
            shm_direct_kmer = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
            kmer_properties = np.ndarray((numberOfKmers,), 
                                         dtype=[("partition",shm_kmer_partition),
                                                ("number",shm_kmer_number),
                                                ("reference",shm_kmer_reference),
                                                ("numberLeft",shm_kmer_numberLeft),
                                                ("numberRight",shm_kmer_numberRight)], 
                                         buffer=shm_kmer.buf)
            direct_properties = np.ndarray((numberOfDirect,), 
                                         dtype=shm_direct_kmer, 
                                         buffer=shm_direct.buf)

            def getRead(counters, structureIterator, stepData,pytablesStorageWorker):
                row = next(structureIterator)
                counters[0]=counters[1]
                counters[1]+=row[0]
                rowData = stepData[counters[0]-counters[2]:counters[1]-counters[2]]
                if counters[1]>counters[3]:
                    counters[2] = counters[3]
                    counters[3] += counters[4]
                    stepData = pytablesStorageWorker.root.readRawData[counters[2]:counters[3]]
                    rowData = np.append(rowData,stepData[0:counters[1]-counters[2]])
                return row, rowData, counters, structureIterator, stepData

            def _getNeighbours(nodeId):
                left  = direct_properties[kmer_properties[nodeId][2]:
                                          kmer_properties[nodeId][2]+kmer_properties[nodeId][3]]
                right = direct_properties[kmer_properties[nodeId][2]+kmer_properties[nodeId][3]:
                                          kmer_properties[nodeId][2]+kmer_properties[nodeId][3]+kmer_properties[nodeId][4]]
                return left,right

            def _getConnections(nodeId,left,right,other=False):
                if other:
                    return (len(left) if not nodeId in left else 0) + (len(right) if not nodeId in right else 0)
                else:
                    return (len(left) if nodeId in left else 0) + (len(right) if nodeId in right else 0)

            def _filterByType(nodes,types,filtered):
                # storage rules:
                # - first only if right splitting (2)
                # - last only if left splitting (1)
                # - pairs only both if first right (2) and second left (1)
                # - type-filtered: don't store pair direct connected
                # - type-filtered: don't store single
                if len(nodes)>2:
                    filtering = [0] * len(nodes)
                    filtering[0] = 1 if types[0]|2==2 else 0
                    filtering[-1] = 1 if types[-1]|1==1 else 0
                    for i in range(1,len(nodes)):
                        if types[i-1]|2==2 and types[i]|1==1:
                            filtering[i-1] = 1
                            filtering[i] = 1
                    s = sum(filtering)
                    if s>2 or (s==2 and len(np.trim_zeros(filtering))>2):
                        filteredNodes = [n for n,f in zip(nodes,filtering) if f==1]
                        filtered.append(filteredNodes)
                return filtered

            def filterReadData(rowData):
                filtered = []
                filteredTypes = []
                #types: 0 - single connected, 1 - splitting to the left, 2 - splitting to the right , 3 - both
                if len(rowData)>1:
                    nodeId = rowData[0]
                    left,right = _getNeighbours(nodeId)
                    newRow = [nodeId]
                    newTypes = [0]
                    for i in range(1,len(rowData)):
                        nNodeId = rowData[i]
                        nLeft,nRight = _getNeighbours(nNodeId)
                        forwardNumber = _getConnections(nNodeId,left,right)
                        backwardNumber = _getConnections(nodeId,nLeft,nRight)
                        #no direct connections (probably read error)
                        if forwardNumber==0 or backwardNumber==0:
                            filtered = _filterByType(newRow,newTypes,filtered)
                            newRow = [nNodeId]
                            newTypes = [0]
                        else:
                            #direct connected
                            newRow.append(nNodeId)
                            #not trivial (backward)
                            if backwardNumber>1:
                                newTypes.append(1)
                            else:
                                newTypes.append(0)
                            #not trivial (forward)
                            if forwardNumber>1:
                                newTypes[-2] |= 2
                        left = nLeft
                        right = nRight
                        nodeId = nNodeId
                    filtered = _filterByType(newRow,newTypes,filtered)
                return filtered


            while not shutdown_event.is_set():
                try:
                    item = queue_rawReads.get(block=True, timeout=1)
                    if item==None:
                        logger.debug("reads ({}): none item".format(os.getpid()))
                        break
                    else: 
                        curr_proc = current_process()
                        pytablesFileWorker = filenameBase+"_tmp_processed_reads_{}.process.h5".format(curr_proc.name)
                        if os.path.exists(pytablesFileWorker):
                            os.remove(pytablesFileWorker)                    
                        with (tables.open_file(item, mode="r") as pytablesStorageWorkerRaw, 
                              tables.open_file(pytablesFileWorker, mode="w") as pytablesStorageWorkerFiltered):
                            #create partition table
                            Storage.create_filteredReads_storage(pytablesStorageWorkerFiltered,
                                                      pytablesStorageWorkerRaw.root.readRawInfo.shape[0],
                                                      pytablesStorageWorkerRaw.root.readRawInfo.shape[0],
                                                      numberOfKmers, numberOfPartitions, maximumReadLength,
                                                      pytablesStorageWorkerRaw.root.readRawInfo.shape[0])
                            readPartitionData = pytablesStorageWorkerFiltered.root.readPartitionData
                            readPartitionInfo = pytablesStorageWorkerFiltered.root.readPartitionInfo
                            readPaired = pytablesStorageWorkerFiltered.root.readPaired                        

                            #compute and store read partition
                            counters = [0,0,0,Storage.stepSizeStorage,Storage.stepSizeStorage] #n0,n1,m0,m1,stepSize
                            #buffer
                            computedPartitionData = []
                            computedPartitionInfo = []
                            computedPaired = []
                            #initialize
                            previousPaired = None
                            try:
                                structureIterator = pytablesStorageWorkerRaw.root.readRawInfo.iterrows()
                                stepData = pytablesStorageWorkerRaw.root.readRawData[counters[2]:counters[3]]
                                while True:
                                    row, rowData, counters, structureIterator, stepData = getRead(
                                        counters, structureIterator, stepData, pytablesStorageWorkerRaw)
                                    assert row[0]==len(rowData)
                                    #filter
                                    filteredRowData = filterReadData(rowData)
                                    #only if filtered set non-empty
                                    if len(filteredRowData)>0:   
                                        rowNodes = []
                                        for filteredEntry in filteredRowData:
                                            rowNodes.extend(filteredEntry)
                                            partitions = [kmer_properties[nodeId][0] for nodeId in filteredEntry]
                                            #store filtered read data
                                            computedPartitionData.extend(filteredEntry)
                                            readPartitions = multimode(partitions)
                                            if len(readPartitions)==1:
                                                computedPartitionInfo.append((len(filteredEntry),readPartitions[0],))
                                            else:
                                                sizes = [kmer_properties[nodeId][1] for nodeId in filteredEntry]
                                                readPartitionSizes = {}
                                                for p in readPartitions:
                                                    selection = np.where(partitions==p)[0]
                                                    readPartitionSizes[p] = max(np.array(sizes)[selection])
                                                computedPartitionInfo.append((
                                                    len(filteredEntry),min(readPartitionSizes, key=readPartitionSizes.get),))
                                        #handle paired data
                                        #TODO: additional filter direct connected???
                                        if row[1]==1:
                                            previousRead = [(nodeId,kmer_properties[nodeId][1]) for nodeId in rowNodes]
                                            previousPaired = min(previousRead,key=lambda x:(x[1],x[0]))[0]
                                        elif row[1]==2 and not previousPaired==None:
                                            currentRead = [(nodeId,kmer_properties[nodeId][1]) for nodeId in rowNodes]
                                            currentPaired = min(currentRead,key=lambda x:(x[1],x[0]))[0]
                                            if not (previousPaired==currentPaired or previousPaired in rowNodes):
                                                computedPaired.append((previousPaired,currentPaired,))
                                                computedPaired.append((currentPaired,previousPaired,))
                                        else:
                                            previousPaired = None
                                    else:
                                        previousPaired = None
                                    #update to storage
                                    if len(computedPartitionData)>Storage.stepSizeStorage:
                                        readPartitionData.append(tuple(computedPartitionData))
                                        readPartitionInfo.append(computedPartitionInfo)
                                        computedPartitionData = []
                                        computedPartitionInfo = []
                                        if len(computedPaired)>0:
                                            computedPaired = sorted(computedPaired,key=lambda x:(x[0],x[1]))
                                            readPaired.append(computedPaired)
                                            computedPaired = []
                            except StopIteration:
                                pass
                            finally:
                                #update last entries to storage
                                if len(computedPartitionData)>0:
                                    readPartitionData.append(tuple(computedPartitionData))
                                    readPartitionInfo.append(computedPartitionInfo)
                                    computedPartitionData = []
                                    computedPartitionInfo = []
                                    if len(computedPaired)>0:
                                        computedPaired = sorted(computedPaired,key=lambda x:(x[0],x[1]))
                                        readPaired.append(computedPaired)
                                        computedPaired = [] 
                                #flush and finish        
                                pytablesStorageWorkerFiltered.flush()
                                logger.debug("reads ({}): filtered {}/{} reads to {}/{}".format(
                                  os.getpid(),
                                  pytablesStorageWorkerRaw.root.readRawInfo.shape[0],
                                  pytablesStorageWorkerRaw.root.readRawData.shape[0],
                                  pytablesStorageWorkerFiltered.root.readPartitionInfo.shape[0],
                                  pytablesStorageWorkerFiltered.root.readPartitionData.shape[0]))
                                #store filtered read files
                                queue_filteredReads.put(pytablesFileWorker)
                except Empty:
                    logger.debug("reads ({}): empty".format(os.getpid()))
                    time.sleep(5)
                    continue      
        finally:
            shm_kmer.close()
            shm_direct.close()

    def workerMergeReads(shutdown_event, queue_ranges, queue_merges, storageReadFiles, 
                          partitionSizes, filenameBase, numberOfKmers, numberOfPartitions, maximumReadLength):

        logger = logging.getLogger("{}.worker.index".format(__name__))
        
        def merge_reads_storage(pytablesFileWorker, storageReadFiles, mergeStart, mergeEnd, partitionPosition):

            #initialise buffers
            buffer = []
            for i in range(mergeStart,mergeEnd+1):
                buffer.append([[],[]])
                
            def getRead(counters, infoIterator, stepData, pytablesMergeSource):
                row = next(infoIterator)
                counters[0]=counters[1]
                counters[1]+=row[0]
                rowData = stepData[counters[0]-counters[2]:counters[1]-counters[2]]
                while counters[1]>counters[3]:                    
                    counters[2] = counters[3]
                    counters[3] += counters[4]
                    stepData = pytablesMergeSource.root.readPartitionData[counters[2]:counters[3]]
                    rowData = np.append(rowData,stepData[0:counters[1]-counters[2]])
                return row, rowData, counters, infoIterator, stepData
            
            with tables.open_file(pytablesFileWorker, mode="a") as pytablesStorageWorker:
                assert len(partitionPosition) == 1 + mergeEnd - mergeStart
                #set partitions
                for i in range(mergeStart,mergeEnd+1):
                    values = [i] * partitionPosition[i-mergeStart][1]
                    startPosition = partitionPosition[i-mergeStart][0]
                    pytablesStorageWorker.root.readPartitionInfo[startPosition:startPosition+len(values)] = values
                #set data and length
                for item in storageReadFiles:
                    with tables.open_file(item, mode="r") as pytablesMergeSource:
                        #compute and store read partition
                        stepSize = 10
                        counters = [0,0,0,stepSize,stepSize] #n0,n1,m0,m1,stepSize
                        #loop
                        try:
                            infoIterator = pytablesMergeSource.root.readPartitionInfo.iterrows()
                            stepData = pytablesMergeSource.root.readPartitionData[counters[2]:counters[3]]
                            while True:
                                row, rowData, counters, infoIterator, stepData = getRead(
                                    counters, infoIterator, stepData, pytablesMergeSource)
                                if row[1]<mergeStart or row[1]>mergeEnd:
                                    continue
                                assert row[0]==len(rowData)
                                partitionId = row[1]-mergeStart
                                buffer[partitionId][0].append(row[0])
                                buffer[partitionId][1].extend(rowData)
                                if len(buffer[partitionId][1])>100000:
                                    startPositionLength = partitionPosition[partitionId][0]
                                    startPositionData = partitionPosition[partitionId][2]
                                    pytablesStorageWorker.root.readPartitionLength[
                                        startPositionLength:
                                        startPositionLength+len(buffer[partitionId][0])] = buffer[partitionId][0]
                                    pytablesStorageWorker.root.readPartitionData[
                                        startPositionData:
                                        startPositionData+len(buffer[partitionId][1])] = buffer[partitionId][1]
                                    partitionPosition[partitionId][0] += len(buffer[partitionId][0])
                                    partitionPosition[partitionId][2] += len(buffer[partitionId][1])
                                    buffer[partitionId] = [[],[]]
                                
                        except StopIteration:
                            pass
                        finally:
                            for i in range(mergeStart,mergeEnd+1):
                                partitionId = i-mergeStart
                                if len(buffer[partitionId][1])>0:
                                    startPositionLength = partitionPosition[partitionId][0]
                                    startPositionData = partitionPosition[partitionId][2]
                                    pytablesStorageWorker.root.readPartitionLength[
                                        startPositionLength:
                                        startPositionLength+len(buffer[partitionId][0])] = buffer[partitionId][0]
                                    pytablesStorageWorker.root.readPartitionData[
                                        startPositionData:
                                        startPositionData+len(buffer[partitionId][1])] = buffer[partitionId][1]
                                    partitionPosition[partitionId][0] += len(buffer[partitionId][0])
                                    partitionPosition[partitionId][2] += len(buffer[partitionId][1])
                                    buffer[partitionId] = [[],[]]
   
        try:
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
                        numberLength = len(str(numberOfPartitions))
                        mergeEnd = min(numberOfPartitions,mergeStart+mergeNumber)-1
                        pytablesFileRange = (filenameBase+"_tmp_read_merge_"+
                                             str(mergeStart).zfill(numberLength)+"_"+
                                        str(mergeEnd).zfill(numberLength)+".process.h5")
                        if os.path.exists(pytablesFileRange):
                            os.remove(pytablesFileRange)
                            
                        #prepare structure to store sorted data
                        partitionPosition = []
                        numberOfReads = 0
                        numberOfData = 0
                        for i in range(mergeStart,mergeEnd+1):
                            partitionPosition.append([numberOfReads,partitionSizes[i][0],numberOfData,partitionSizes[i][1]])
                            numberOfReads+=partitionSizes[i][0]
                            numberOfData+=partitionSizes[i][1]
                        #only merge if reads available
                        if (numberOfReads>0) and (numberOfData>0):
                            with tables.open_file(pytablesFileRange, mode="a") as pytablesStorageRange:
                                Storage.create_mergeReads_storage(pytablesStorageRange,
                                                      numberOfReads, numberOfData, numberOfKmers, 
                                                      numberOfPartitions, maximumReadLength)
                            #merge
                            merge_reads_storage(pytablesFileRange, storageReadFiles, 
                                                mergeStart, mergeEnd, partitionPosition)
                            #now the file can be released for final merge
                            queue_merges.put(pytablesFileRange) 
                except Empty:
                    time.sleep(1)
                    continue
        except Exception as ex:
            logger.error("merges ({}): problem with worker: {}".format(os.getpid(),ex))
        finally:
            pass
    
    
    """
    Combine merged files: this should be relatively easy, handled by single process
    """
    def combineDirectMerges(mergeFiles,pytablesStorage,numberOfKmers,maximumFrequency):
        
        nCycle = 0
        nReversal = 0
        nDirect = 0
        nPaired = 0
        
        #get dimensions
        sortedMergeFiles = [{"filename": mergeFile} for mergeFile in sorted(mergeFiles)]
        for i in range(len(sortedMergeFiles)):
            with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                sortedMergeFiles[i]["ncycle"] = pytables_merge.root.cycle.shape[0]
                sortedMergeFiles[i]["nreversal"] = pytables_merge.root.reversal.shape[0]
                sortedMergeFiles[i]["ndirect"] = pytables_merge.root.direct.shape[0]
                sortedMergeFiles[i]["ndeleteDirect"] = pytables_merge.root.deleteDirect.shape[0]
                sortedMergeFiles[i]["npaired"] = pytables_merge.root.tmpPaired.shape[0]                
                nCycle += sortedMergeFiles[i]["ncycle"]
                nReversal += sortedMergeFiles[i]["nreversal"]
                nDirect += sortedMergeFiles[i]["ndirect"]
                nPaired += sortedMergeFiles[i]["npaired"]
                                
        #create temporary storage with dimensions
        Storage.create_mergeDirect_storage(pytablesStorage, numberOfKmers, maximumFrequency, 
                                     nCycle, nReversal, nDirect, nPaired)
        tableCycle = pytablesStorage.root.cycle
        tableReversal = pytablesStorage.root.reversal
        tableDirect = pytablesStorage.root.direct
        tableDeleteDirect = pytablesStorage.root.deleteDirect
        tablePaired = pytablesStorage.root.tmpPaired
        
        maximumCycleLength = 0
        maximumCycleNumber = 0
        maximumReversalLength = 0
        maximumReversalNumber = 0
        maximumDirectDistance = 0
        maximumDirectNumber = 0
        
        #merge delete entries
        for i in range(len(sortedMergeFiles)):
            with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                tableDeleteDirect.append(pytables_merge.root.deleteDirect[:])
        pytablesStorage.flush()
        tableDeleteDirect.cols.fromLink.create_csindex()
        pytablesStorage.flush()
                
        for i in range(len(sortedMergeFiles)):
            with tables.open_file(sortedMergeFiles[i]["filename"], mode="r") as pytables_merge:
                for j in range(0,sortedMergeFiles[i]["ncycle"],Storage.stepSizeStorage):
                    reducedCycleDataBlock = pytables_merge.root.cycle[j:j+Storage.stepSizeStorage]
                    maximumCycleNumber = max(maximumCycleNumber,max(reducedCycleDataBlock["number"]))
                    maximumCycleLength = max(maximumCycleLength,max(reducedCycleDataBlock["minimumLength"]))
                    tableCycle.append(reducedCycleDataBlock)
                for j in range(0,sortedMergeFiles[i]["nreversal"],Storage.stepSizeStorage):
                    reducedReversalDataBlock = pytables_merge.root.reversal[j:j+Storage.stepSizeStorage]
                    maximumReversalNumber = max(maximumReversalNumber,max(reducedReversalDataBlock["number"]))
                    maximumReversalLength = max(maximumReversalLength,max(reducedReversalDataBlock["minimumLength"]))
                    tableReversal.append(reducedReversalDataBlock)
                for j in range(0,sortedMergeFiles[i]["ndirect"],Storage.stepSizeStorage):
                    directDataBlock = pytables_merge.root.direct[j:j+Storage.stepSizeStorage]
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
                for j in range(0,sortedMergeFiles[i]["npaired"],Storage.stepSizeStorage):
                    reducedPairedDataBlock = pytables_merge.root.tmpPaired[j:j+Storage.stepSizeStorage]
                    tablePaired.append(reducedPairedDataBlock)

        tableCycle.attrs["maximumNumber"] = maximumCycleNumber
        tableCycle.attrs["maximumLength"] = maximumCycleLength
        tableReversal.attrs["maximumNumber"] = maximumReversalNumber
        tableReversal.attrs["maximumLength"] = maximumReversalLength
        tableDirect.attrs["maximumDistance"] = maximumDirectDistance
        tableDirect.attrs["maximumNumber"] = maximumDirectNumber
        pytablesStorage.flush()
        #sort paired data
        tablePaired.cols.fromLink.create_csindex()
        pytablesStorage.flush()
        
    def combineFilteredPairs(mergeFiles,pytablesStorage,numberOfKmers):
        logger = logging.getLogger(__name__) 
        #collect, sort and combine data
        dataList = []
        for item in mergeFiles:
            with tables.open_file(item, mode="r") as pytablesSource:
                size = pytablesSource.root.readPaired.shape[0]
                itemData = pytablesSource.root.readPaired[0:size]
                dataList.append(np.unique(itemData))
                del itemData
        data = np.unique(np.concatenate(dataList))
        del dataList
        #store data
        numberOfPaired = len(data)
        logger.info("store "+str(numberOfPaired)+" paired")
        pairedTable = pytablesStorage.create_table(pytablesStorage.root, 
            "readPaired",{
            "fromLink": haplotyping.index.Database.getTablesUint(numberOfKmers,0),
            "toLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
        }, "Paired reads", expectedrows=numberOfPaired)
        pairedTable.append(data)
        pairedTable.flush()
        
        
    def combineReadMerges(mergeFiles,pytablesStorage,numberOfKmers,
                                  numberOfPartitions, maximumOriginalReadLength):
        numberOfReads = 0
        numberOfData = 0
        partitionIndex = [0] * numberOfPartitions
        #compute number of reads for each partition
        for item in mergeFiles:
            with tables.open_file(item, mode="r") as pytablesSource:
                numberOfReads+=pytablesSource.root.readPartitionInfo.shape[0]
                numberOfData+=pytablesSource.root.readPartitionData.shape[0]
                for i in range(0,numberOfReads,Storage.stepSizeStorage):
                    stepData = pytablesSource.root.readPartitionInfo[i:i+Storage.stepSizeStorage]
                    for p in stepData:
                        partitionIndex[p]+=1
        #only if reads found
        if numberOfReads>0:
            Storage.create_filteredMergedReads_storage(pytablesStorage,numberOfReads,numberOfData,
                      numberOfKmers, numberOfPartitions, maximumOriginalReadLength)
            #variables for final position reads and data
            tReads = 0
            tData = 0
            #final storage
            readPartition=pytablesStorage.root.readPartition
            readPartitionInfo=pytablesStorage.root.readPartitionInfo
            readPartitionData=pytablesStorage.root.readPartitionData
            #statistics
            maximumReadLength = 0
            maximumTotalReadLength = 0
            maximumReadNumber = 0
            maximumTotalReadNumber = 0
            #processed partition
            partition = 0            
            for item in mergeFiles:
                with tables.open_file(item, mode="r") as pytablesSource:
                    #get sizes source files
                    sInfo=pytablesSource.root.readPartitionInfo.shape[0]
                    sLength=pytablesSource.root.readPartitionLength.shape[0]
                    sData=pytablesSource.root.readPartitionData.shape[0]
                    #variables for source position reads and data
                    pReads=0
                    pData=0
                    while pReads<sInfo:
                        #collect data, can be empty
                        nReads = partitionIndex[partition]
                        partitionLength = pytablesSource.root.readPartitionLength[pReads:pReads+nReads]
                        nData = sum(partitionLength)
                        partitionData = pytablesSource.root.readPartitionData[pData:pData+nData]
                        #filter
                        partitionDataIndex = {}
                        partitionDataInfo = []
                        fData = 0
                        mData = 0
                        for i in range(nReads):
                            rowData = tuple(partitionData[fData:fData+partitionLength[i]])
                            rowKey = hash(rowData)
                            if rowKey in partitionDataIndex:
                                partitionDataInfo[partitionDataIndex[rowKey]][1]+=1
                            else:
                                partitionDataIndex[rowKey] = len(partitionDataInfo)
                                partitionDataInfo.append([partitionLength[i],1])
                                readPartitionData.append(rowData)
                                mData+=partitionLength[i]
                            fData+=partitionLength[i]
                        if len(partitionDataInfo)>0:
                            readPartitionInfo.append([tuple(item) for item in partitionDataInfo])
                            maximumReadLength = max(maximumReadLength,max([item[0] 
                                                   for item in partitionDataInfo]))
                            maximumTotalReadLength = max(maximumTotalReadLength,sum([item[0] 
                                                   for item in partitionDataInfo]))
                            maximumReadNumber = max(maximumReadNumber,max([item[1] 
                                                   for item in partitionDataInfo]))
                            maximumTotalReadNumber = max(maximumTotalReadNumber,len(partitionDataInfo))
                        #always store information for partition
                        readPartition.append([(tData,mData,tReads,len(partitionDataInfo))])
                        pReads+=nReads
                        pData+=nData
                        tReads+=len(partitionDataInfo)
                        tData+=mData
                        partition+=1
            while partition<numberOfPartitions:
                readPartition.append((tData,0,tReads,0))
                partition+=1
            readPartitionInfo.attrs["maximumReadLength"] = maximumReadLength
            readPartitionInfo.attrs["maximumTotalReadLength"] = maximumTotalReadLength
            readPartitionInfo.attrs["maximumReadNumber"] = maximumReadNumber
            readPartitionInfo.attrs["maximumTotalReadNumber"] = maximumTotalReadNumber

    """
    Store direct data in the final database, handled by single process
    """    
    def storeMergedDirect(h5file,pytablesStorage,numberOfKmers,minimumFrequency):
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
                                                     dtype=dtCycle, chunks=None, 
                                                     compression="gzip", compression_opts=9)
        maximumCycleLength = 0
        maximumCycleNumber = 0
        for i in range(0,numberOfCycles,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.cycle[i:i+Storage.stepSizeStorage]                    
            dsCycle[i:i+Storage.stepSizeStorage] = stepData
            for row in stepData:
                ckmerRow = dsCkmer[row[0]]
                ckmerRow[6] = row[2]
                maximumCycleLength = max(maximumCycleLength,row[1])
                maximumCycleNumber = max(maximumCycleNumber,row[2])
                dsCkmer[row[0]] = ckmerRow    
        h5file["/config/"].attrs["numberCycles"]=numberOfCycles        
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
                                                        dtype=dtReversal, chunks=None, 
                                                        compression="gzip", compression_opts=9)
        maximumReversalLength = 0
        maximumReversalNumber = 0
        for i in range(0,numberOfReversals,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.reversal[i:i+Storage.stepSizeStorage]
            dsReversal[i:i+Storage.stepSizeStorage] = stepData
            for row in stepData:
                ckmerRow = dsCkmer[row[0]]
                ckmerRow[7] = row[2]
                maximumReversalLength = max(maximumReversalLength,row[1])
                maximumReversalNumber = max(maximumReversalNumber,row[2])
                dsCkmer[row[0]] = ckmerRow     
        h5file["/config/"].attrs["numberReversals"]=numberOfReversals
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
                                                      dtype=dtDirect, chunks=None, 
                                                      compression="gzip", compression_opts=9)
        #process direct relations
        directCounter = 0
        previousStepData = []
        directGood = 0
        directCoverage = 0
        directProblematic = 0
        for i in range(0,numberOfDirectRelations,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.direct[i:i+Storage.stepSizeStorage]
            #move a bit to get all fromLinks in the selection
            if len(previousStepData)>0:
                stepData = np.concatenate((previousStepData, stepData))
            if (i+Storage.stepSizeStorage-1)<numberOfDirectRelations:
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

        logger.info("store {} direct connections".format(numberOfDirectRelations))
        logger.info("{} low coverage and {} problematic".format(directCoverage,directProblematic))

        # store histogram direct distances
        if len(frequencyHistogram["distance"])>0:
            maximumDistance = max(frequencyHistogram["distance"].keys())
            maximumNumber = max(frequencyHistogram["distance"].values())
            dtypeFrequencyHistogramDistanceList=[
                        ("distance",haplotyping.index.Database.getUint(maximumDistance)),
                        ("number",haplotyping.index.Database.getUint(maximumNumber)),]
            dtFrequencyHistogramDistance=np.dtype(dtypeFrequencyHistogramDistanceList)
            dsFrequencyHistogramDistance=h5file["/histogram/"].create_dataset("distance",
                  (len(frequencyHistogram["distance"]),), dtype=dtFrequencyHistogramDistance, chunks=None, 
                  compression="gzip", compression_opts=9)
            logger.info("store {} entries direct distances histogram".format(
                len(frequencyHistogram["distance"])))
            #store histogram
            dsFrequencyHistogramDistance[0:len(frequencyHistogram["distance"])] = list(
                sorted(frequencyHistogram["distance"].items()))
            
    def storeMergedReads(h5file,pytablesStorage,numberOfKmers,numberOfPartitions):
        logger = logging.getLogger(__name__)    
        
        #paired
        ckmerIndex = [[0,0] for i in range(numberOfKmers)]
        numberOfPaired = pytablesStorage.root.readPaired.shape[0]
        dtypePairedList=[("fromLink",haplotyping.index.Database.getUint(numberOfKmers)),
                         ("toLink",haplotyping.index.Database.getUint(numberOfKmers))]
        dtPaired=np.dtype(dtypePairedList)
        dsPaired=h5file["/relations/"].create_dataset("paired",(numberOfPaired,), 
                                                      dtype=dtPaired, chunks=None, 
                                                      compression="gzip", compression_opts=9)
        for i in range(0,numberOfPaired,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.readPaired[i:i+Storage.stepSizeStorage]
            for j in range(len(stepData)):
                row = stepData[j]
                if ckmerIndex[row[0]][1]>0:
                    ckmerIndex[row[0]][1] += 1
                else:
                    ckmerIndex[row[0]] = [i+j,1]
            dsPaired[i:i+len(stepData)] = stepData
        logger.info("store {} paired connections".format(numberOfPaired))
        
        #update kmer properties
        dsCkmer = h5file["/split/ckmer"]
        for i in range(0,numberOfKmers,Storage.stepSizeStorage):
            stepData = dsCkmer[i:i+Storage.stepSizeStorage]
            stepData["paired"] = [tuple(item) for item in ckmerIndex[i:i+len(stepData)]]
            dsCkmer[i:i+Storage.stepSizeStorage] = stepData
        
        #reads
        numberOfReadPartition = pytablesStorage.root.readPartition.shape[0]
        numberOfReadPartitionData = pytablesStorage.root.readPartitionData.shape[0]
        numberOfReadPartitionInfo = pytablesStorage.root.readPartitionInfo.shape[0]
        
        maxReadLength = pytablesStorage.root.readPartitionInfo.attrs["maximumReadLength"]
        maxTotalReadLength = pytablesStorage.root.readPartitionInfo.attrs["maximumTotalReadLength"]
        maxReadNumber = pytablesStorage.root.readPartitionInfo.attrs["maximumReadNumber"]
        maxTotalReadNumber = pytablesStorage.root.readPartitionInfo.attrs["maximumTotalReadNumber"]
        
        dsReadData=h5file["/relations/"].create_dataset("readData",(numberOfReadPartitionData,), 
                                                      dtype=haplotyping.index.Database.getUint(numberOfKmers), 
                                                      chunks=None, compression="gzip", compression_opts=9)
        for i in range(0,numberOfReadPartitionData,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.readPartitionData[i:i+Storage.stepSizeStorage]
            dsReadData[i:i+len(stepData)] = stepData
        logger.info("store {} read data points".format(numberOfReadPartitionData))
        dtypeReadInfoList=[("length",haplotyping.index.Database.getUint(maxReadLength)),
                           ("number",haplotyping.index.Database.getUint(maxReadNumber))]
        dtReadInfo=np.dtype(dtypeReadInfoList)
        dsReadInfo=h5file["/relations/"].create_dataset("readInfo",(numberOfReadPartitionInfo,), 
                                                      dtype=dtReadInfo, chunks=None, 
                                                      compression="gzip", compression_opts=9)
        for i in range(0,numberOfReadPartitionInfo,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.readPartitionInfo[i:i+Storage.stepSizeStorage]
            dsReadInfo[i:i+len(stepData)] = stepData
        logger.info("store {} read info".format(numberOfReadPartitionInfo))
        dtypeReadPartitionList=[("readData",[
                                    ("link",haplotyping.index.Database.getUint(numberOfReadPartitionData)),
                                    ("number",haplotyping.index.Database.getUint(maxTotalReadLength))]),
                                ("readInfo",[
                                    ("link",haplotyping.index.Database.getUint(numberOfReadPartitionInfo)),
                                    ("number",haplotyping.index.Database.getUint(maxTotalReadNumber))])]
        dtReadPartition=np.dtype(dtypeReadPartitionList)
        dsReadPartition=h5file["/relations/"].create_dataset("readPartition",(numberOfReadPartition,), 
                                                      dtype=dtReadPartition, chunks=None, 
                                                      compression="gzip", compression_opts=9)
        for i in range(0,numberOfReadPartition,Storage.stepSizeStorage):
            stepData = pytablesStorage.root.readPartition[i:i+Storage.stepSizeStorage]
            dsReadPartition[i:i+len(stepData)] = [((item[0],item[1]),(item[2],item[3])) for item in stepData]
        logger.info("store {} read partitions".format(numberOfReadPartition))
        
        
        
    """
    Partition k-mers
    """    
    def partitionKmers(h5file,pytablesStorage,maxNumberOfPartitions):
        #initialize
        logger = logging.getLogger(__name__)  
        edges = []
        
        #compute edges for partition graph
        numberOfKmers = h5file["split"]["ckmer"].shape[0]
        numberOfBases = h5file["split"]["base"].shape[0]
        numberOfDirect = h5file["relations"]["direct"].shape[0]
        maxNumberOfPartitions = int(numberOfKmers ** (2/3))
        previousVertex = 0
        previousNeighbours = []
        logger.debug("create graph for partitioning")
        for i in range(0,numberOfDirect,Storage.stepSizeStorage):
            block = h5file["relations"]["direct"][i:i+Storage.stepSizeStorage]
            logger.debug("add edges: {}-{} of {} ({}%)".format(
                i,i+len(block),numberOfDirect,int(100*(i+len(block))/numberOfDirect)))
            for j in range(len(block)):
                while block[j][0][0]>previousVertex:
                    edges.append(tuple(previousNeighbours))
                    previousVertex+=1
                    previousNeighbours = []
                if not block[j][1][0] in previousNeighbours:
                    previousNeighbours.append(block[j][1][0])
        while numberOfKmers>previousVertex:
            edges.append(tuple(previousNeighbours))
            previousVertex+=1
            previousNeighbours = []

        #check symmetry graph
        n_self = 0
        n_asymmetric = 0
        for i in range(len(edges)):
            if i in edges[i]:
                n_self+=1
            for j in edges[i]:
                if not i in edges[j]:
                    n_asymmetric+=1
        logger.debug("detected {} self connected nodes and {} asymmetric connections".format(
            n_self,n_asymmetric))
            
        #update disconnected or only self connected
        disconnected = [i for i in range(len(edges)) if len(edges[i])==0]
        selfconnected = [i for i in range(len(edges)) if len(edges[i])==1 and i in edges[i]]
        logger.debug("connect {} disconnected and {} only self-connected nodes to base connections".format(
            len(disconnected),len(selfconnected)))

        #combine
        disconnected = disconnected + selfconnected
        disconnected.sort()

        #create memory for links with bases
        bases = h5file["split/base"]
        ckmers = h5file["split/ckmer"]
        disconnectedRows = ckmers[disconnected]
        for ckmerId,ckmerRow in zip(disconnected,disconnectedRows):
            if ckmerRow[1]==b"l":
                baseIds = [ckmerRow[3][0]]
            elif ckmerRow[1]==b"r":
                baseIds = [ckmerRow[3][1]]
            elif ckmerRow[1]==b"b":
                baseIds = [ckmerRow[3][0],ckmerRow[3][1]]
            else:
                baseIds = []
            ckmerLinks = set()
            for baseId in baseIds:
                ckmerLinks.update([baseEntry[1] for baseEntry in bases[baseId][2] if baseEntry[0]>0])
                ckmerLinks.remove(ckmerId) 
                #keep it symmetric
                for ckmerLink in ckmerLinks:
                    if not ckmerLink in edges[ckmerId]:
                        edges[ckmerId] = edges[ckmerId] + (ckmerLink,)
                    if not ckmerId in edges[ckmerLink]:
                        edges[ckmerLink] = edges[ckmerLink] + (ckmerId,)

        #connected components
        g = nk.graph.Graph(weighted=False)
        g.addNodes(len(edges))
        for i in range(len(edges)):
            for j in edges[i]:
                g.addEdge(i,j)
        cc = nk.components.ConnectedComponents(g)
        cc.run()
        ccList = cc.getComponents()
        logger.debug("found {} connected components".format(len(ccList)))

        #analyse connected components
        nc_trivial = 0
        nc_small = 0
        trivialPartitions = []
        smallPartitions = []
        selectedNodes = []
        partitionSize = int(len(edges)/maxNumberOfPartitions)
        for component in ccList:
            if len(component)<=2:
                nc_trivial+=1
                if (len(trivialPartitions)==0) or len(trivialPartitions[-1])>=partitionSize:
                    trivialPartitions.append([])
                trivialPartitions[-1].extend(component)
            elif len(component)<=partitionSize:
                nc_small+=1
                smallPartitions.append(component)
            else:
                selectedNodes.extend(component)
        logger.debug("distribute {} trivial components with at most 2 nodes over {} partitions".format(
            nc_trivial,len(trivialPartitions)))
        logger.debug("assign {} small components with at most {} nodes to separate partitions".format(
            nc_small,partitionSize))
        
        #compute final partitions
        partitions = [-1] * len(edges)
        counts = []
        for partition in trivialPartitions:
            for node in partition:
                partitions[node] = len(counts)
            counts.append(len(partition))
        for partition in smallPartitions:
            for node in partition:
                partitions[node] = len(counts)
            counts.append(len(partition))
        numberOfPartitions = len(counts)

        #compute normal partitions
        if len(selectedNodes)>0:
            #recompute number of partitions
            newMaxNumberOfPartitions = int(len(selectedNodes)/partitionSize)
            logger.debug("compute metis graph partitioning with at most {} partitions".format(
                newMaxNumberOfPartitions))
            #compute selected graph edges
            selectedNodes.sort()
            selectedMap = {selectedNodes[i]: i for i in range(len(selectedNodes))}
            selectedEdges = [tuple([selectedMap[j] for j in edges[selectedNodes[i]]]) 
                             for i in range(len(selectedNodes))]

            m = metis.adjlist_to_metis(selectedEdges)
            del selectedEdges
            (objval,selectedPartitions) = metis.part_graph(m,newMaxNumberOfPartitions)
            (selectedValues, selectedPartitions, selectedCounts) = np.unique(selectedPartitions, 
                                                 return_inverse=True, return_counts=True) 
            numberOfSelectedPartitions = len(selectedCounts)
            logger.debug("computed {} regular partitions".format(numberOfSelectedPartitions))
            #update final partitions
            for i in range(len(selectedPartitions)):
                partitions[selectedNodes[i]] = selectedPartitions[i] + numberOfPartitions
            counts = counts + ([0] * numberOfSelectedPartitions)
            for i in range(len(selectedCounts)):
                counts[i + numberOfPartitions] = selectedCounts[i]
            numberOfPartitions = len(counts)
            
        #some checks
        assert len([p for p in partitions if p<0])==0
        assert sum(counts) == len(edges)

        #store partition in ckmer properties
        dsCkmer = h5file["/split/ckmer"]
        assert numberOfKmers==len(partitions)
        for i in range(0,numberOfKmers,Storage.stepSizeStorage):
            stepData = dsCkmer[i:i+Storage.stepSizeStorage] 
            stepData["partition"] = partitions[i:i+Storage.stepSizeStorage] 
            dsCkmer[i:i+Storage.stepSizeStorage] = stepData
        
        #store partition size
        h5file["/config/"].attrs["numberPartitions"]=numberOfPartitions  
        dtypePartitionList=[("ckmer",haplotyping.index.Database.getUint(numberOfKmers))]
        dtPartition=np.dtype(dtypePartitionList)
        dsPartition=h5file["/histogram/"].create_dataset("partition",(numberOfPartitions,), 
                                                        dtype=dtPartition, chunks=None, 
                                                        compression="gzip", compression_opts=9)
        dsPartition[0:numberOfPartitions] = counts
        
        return numberOfPartitions 
        
        



