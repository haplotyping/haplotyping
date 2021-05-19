import logging, h5py, tables, gzip, time
import os, sys, tempfile, ahocorasick, pickle
import re, haplotyping
import numpy as np

class Relations:
    
    """
    Internal use, parse read files and store results in database
    """
    
    def __init__(self, readFiles,pairedReadFiles, h5file, filenameBase, debug):
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("parse "+str(len(readFiles))+" readfiles and "+str(2*len(pairedReadFiles))+" paired readfiles")
        
        self.stepSizeStorage=10000
        self.debug = debug
        self.filenameBase = filenameBase
        
        #set variables
        self.k = h5file["/config"].attrs["k"]
        self.h5file = h5file
        self.readFiles = readFiles
        self.pairedReadFiles = pairedReadFiles
        self.numberOfKmers = self.h5file["/split/ckmer"].shape[0]
        self.numberOfDirectRelations = 0
        self.numberOfCycles = 0
        self.numberOfReversals = 0
        self.maxNumberofLinkDigits = int(np.ceil(np.log10(self.numberOfKmers)))
        
        self.readLengthMinimum=None
        self.readLengthMaximum=None
        self.readPairedTotal=0
        self.readUnpairedTotal=0
        self.readTotal=0
        self.processReadsTime=0
        
        self.maxDirectDistance = 0
        self.maxDirectNumber = 0  
        self.maxConnectedIndexSize = 0
        self.maxConnectedIndexLength = 0
        self.maxConnectedIndexNumber = 0
        self.maxConnectedPairNumber = 0
        self.maxCycleLength = 0
        self.maxCycleNumber = 0
        self.maxReversalLength = 0
        self.maxReversalNumber = 0
        
        #get config
        self.minimumFrequency = self.h5file["/config"].attrs["minimumFrequency"]
        
        #check existence group
        if not "/relations" in h5file:
            h5file.create_group("/relations")
        
        #create direct relations dataset
        if "/relations/direct" in h5file:
            self._logger.warning("direct relation dataset already exists in hdf5 storage")
        elif "/relations/connected" in h5file:
            self._logger.warning("connected relation dataset already exists in hdf5 storage")
        elif "/config/pairedReads" in h5file or "/config/unpairedReads" in h5file:
            self._logger.warning("readfiles already exist in hdf5 storage")
        else:
            
            try:
                
                if self.debug:
                    self.tmpDirectory = None
                    pytablesMainFile = self.filenameBase+"_tmp_relations_main.h5"
                    pytablesDumpFile = self.filenameBase+"_tmp_relations_dump.h5"
                    if os.path.exists(pytablesMainFile):
                        os.remove(pytablesMainFile)
                    if os.path.exists(pytablesDumpFile):
                        os.remove(pytablesDumpFile)
                    self._logger.debug("store temporary in "+pytablesMainFile) 
                    self._logger.debug("store temporary in "+pytablesDumpFile) 
                else:
                    self.tmpDirectory = tempfile.TemporaryDirectory()
                    pytablesMainFile = self.tmpDirectory.name+"/kmer.data_tmp_relations_main.h5"
                    pytablesDumpFile = self.tmpDirectory.name+"/kmer.data_tmp_relations_dump.h5"
                
                with tables.open_file(pytablesMainFile, mode="w", 
                                      title="Temporary main storage") as self.pytables_main_storage, \
                     tables.open_file(pytablesDumpFile, mode="w", 
                                      title="Temporary dump storage") as self.pytables_dump_storage :
                    #initialisation
                    self._getAutomaton()
                    
                    #process
                    self._processReadFiles()
                    self._group()
                    self._store()
                    
                    #flush
                    self.h5file.flush()
               
            #except:
            #    self._logger.error("problem occurred while constructing relations")
            finally:
                if not self.debug:
                    self.tmpDirectory.cleanup()
                
              
    def _processReadFiles(self):    
        #create temporary storage
        tableDumpDirectDef = {
            "sortKey": tables.StringCol(2+(2*self.maxNumberofLinkDigits),pos=0),
            "fromCkmerLink": self._getTablesUint(self.numberOfKmers,1),
            "fromDirection": tables.StringCol(1,pos=2),
            "toCkmerLink": self._getTablesUint(self.numberOfKmers,3),
            "toDirection": tables.StringCol(1,pos=4),
            "distance": tables.UInt16Col(pos=5),
        }
        tableDumpConnectedDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
        }
        tableDumpConnectedIndexDef = {
            "hashKey": tables.StringCol(21,pos=0),
            "sizeConnected": tables.UInt64Col(pos=1),
            "lengthConnected": tables.UInt64Col(pos=2),
            "firstDumpConnectedLink": tables.UInt64Col(pos=3),
        }
        tableDumpConnectedPairDef = {
            "hashKey0": tables.StringCol(21,pos=0),
            "hashKey1": tables.StringCol(21,pos=1),
        }
        tableDumpCycleDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
            "length": tables.UInt64Col(pos=1),
        }
        tableDumpReversalDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
            "length": tables.UInt64Col(pos=1),
        }
        self.tableDumpDirect = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root, 
                                "dumpDirect", tableDumpDirectDef, 
                                "Temporary to dump direct relations")
        self.tableDumpConnected = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root, 
                                "dumpConnected", tableDumpConnectedDef, 
                                "Temporary to dump connections")
        self.numberTableDumpConnected = 0
        self.tableDumpConnectedIndex = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root,
                                "dumpConnectedIndex", tableDumpConnectedIndexDef, 
                                "Temporary for index dumped connections")
        self.tableDumpConnectedPair = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root,
                                "dumpConnectedPair", tableDumpConnectedPairDef, 
                                "Temporary for dumped connection pairs")
        self.tableDumpCycle = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root,
                                "dumpCycle", tableDumpCycleDef, 
                                "Temporary for dumped cycles")
        self.tableDumpReversal = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root,
                                "dumpReversal", tableDumpReversalDef, 
                                "Temporary for dumped reversals")
        
        #process and register unpaired read files
        dtypeList = [("file","S255"),("readLength","uint64"),
                     ("readNumber","uint64"),("processTime","uint32")]
        ds = self.h5file["/config/"].create_dataset("unpairedReads",(len(self.readFiles),),
                                          dtype=np.dtype(dtypeList),chunks=None)
        for i in range(len(self.readFiles)):
            (readLength,readNumber,processTime) = self._processReadFile(self.readFiles[i])
            ds[i] = (self.readFiles[i],
                     readLength,readNumber,int(processTime))


        #process and register paired read files
        dtypeList = [("file0","S255"),("file1","S255"),("readLength","uint64"),
                     ("readNumber","uint64"),("processTime","uint32")]
        ds = self.h5file["/config/"].create_dataset("pairedReads",(len(self.pairedReadFiles),),
                                          dtype=np.dtype(dtypeList),chunks=None)
        for i in range(len(self.pairedReadFiles)):
            (readLength,readNumber,processTime) = self._processPairedReadFiles(self.pairedReadFiles[i][0],
                                                                         self.pairedReadFiles[i][1])
            ds[i] = (self.pairedReadFiles[i][0],self.pairedReadFiles[i][1],
                     readLength,readNumber,int(processTime))
            
        self.tableDumpDirect.flush()
        self.tableDumpConnected.flush()
        self.tableDumpConnectedIndex.flush()
        self.tableDumpConnectedPair.flush()
        self.tableDumpCycle.flush()
        self.tableDumpReversal.flush()
        
        
            
                
    def _processReadFile(self, filename: str):
        startTime = time.time()
        with gzip.open(filename, "rt") as f:
            readLengthMinimum=None
            readLengthMaximum=None
            readNumber=0
            while True:               
                f.readline()
                sequence = f.readline().rstrip()  
                f.readline()
                f.readline()
                if sequence:
                    readLengthMinimum=(len(sequence) if readLengthMinimum==None 
                                       else min(readLengthMinimum,len(sequence)))
                    readLengthMaximum=(len(sequence) if readLengthMaximum==None 
                                       else max(readLengthMaximum,len(sequence)))
                    readNumber+=1
                    matchesList=self._computeMatchesList(sequence)
                    for matches in matchesList:
                        self._processMatches(matches)
                    self._processConnected(sum(matchesList,[]),False)
                else:
                    break
            endTime = time.time()
            if readLengthMinimum>0:
                self.readLengthMinimum=(readLengthMinimum if self.readLengthMinimum==None 
                                        else min(self.readLengthMinimum,readLengthMinimum))
            self.readLengthMaximum=(readLengthMaximum if self.readLengthMaximum==None 
                                    else max(self.readLengthMaximum,readLengthMaximum))
            self.readUnpairedTotal+=readNumber
            self.readTotal+=readNumber
            self.processReadsTime+=endTime-startTime
            return (readLengthMaximum,readNumber,endTime-startTime)

    def _processPairedReadFiles(self, filename0: str, filename1: str):
        startTime = time.time()
        
        def storePairedConnected(hashKey0,hashKey1):
            if hashKey0 and hashKey1 and not hashKey0==hashKey1:
                dumpConnectedPair = self.tableDumpConnectedPair.row
                dumpConnectedPair["hashKey0"]=hashKey0
                dumpConnectedPair["hashKey1"]=hashKey1
                dumpConnectedPair.append()
                dumpConnectedPair["hashKey0"]=hashKey1
                dumpConnectedPair["hashKey1"]=hashKey0
                dumpConnectedPair.append()
        
        with gzip.open(filename0, "rt") as f0, gzip.open(filename1, "rt") as f1:
            readLengthMinimum=None
            readLengthMaximum=None
            readNumber=0
            while True:               
                #first of pair
                f0.readline()
                sequence0 = f0.readline().rstrip()  
                f0.readline()
                f0.readline()
                #second of pair
                f1.readline()
                sequence1 = haplotyping.General.reverse_complement(f1.readline().rstrip())
                f1.readline()
                f1.readline()
                #process
                if sequence0 and sequence1:
                    readLengthMinimum=(min(len(sequence0),len(sequence1)) if readLengthMinimum==None 
                                           else min(readLengthMinimum,len(sequence0),len(sequence1)))
                    readLengthMaximum=(max(len(sequence0),len(sequence1)) if readLengthMaximum==None 
                                       else max(readLengthMaximum,len(sequence0),len(sequence1)))
                    readNumber+=2      
                    matchesList0=self._computeMatchesList(sequence0)
                    matchesList1=self._computeMatchesList(sequence1)
                    connectedHashKeys0=set()
                    connectedHashKeys1=set()
                    for matches0 in matchesList0:
                        self._processMatches(matches0)
                    hashKey0 = self._processConnected(sum(matchesList0,[]),len(matchesList1)>0)
                    for matches1 in matchesList1:
                        self._processMatches(matches1)
                    hashKey1 = self._processConnected(sum(matchesList1,[]),len(matchesList0)>0)
                    #connecting real pair
                    storePairedConnected(hashKey0,hashKey1) 
                else:
                    break
            endTime = time.time()
            if readLengthMinimum>0:
                self.readLengthMinimum=(readLengthMinimum if self.readLengthMinimum==None 
                                        else min(self.readLengthMinimum,readLengthMinimum))
            self.readLengthMaximum=(readLengthMaximum if self.readLengthMaximum==None 
                                    else max(self.readLengthMaximum,readLengthMaximum))
            self.readPairedTotal+=readNumber
            self.readTotal+=readNumber
            self.processReadsTime+=endTime-startTime
            return (readLengthMaximum,readNumber,endTime-startTime)
        
    def _group(self):
        self._groupDirect()
        self._groupConnected()
        self._indexConnected()
        
    def _groupDirect(self):
        #create temporary storage
        tableGroupedDirectDef = {
            "sortKey": tables.StringCol(2+(2*self.maxNumberofLinkDigits),pos=0),
            "fromCkmerLink": self._getTablesUint(self.numberOfKmers,1),
            "fromDirection": tables.StringCol(1,pos=2),
            "toCkmerLink": self._getTablesUint(self.numberOfKmers,3),
            "toDirection": tables.StringCol(1,pos=4),
            "distance": self._getTablesUint(self.maxDirectDistance,5),
            "number": tables.UInt16Col(pos=6),
        }
        tableGroupedDirectErrorDef = {
            "sortKey": tables.StringCol(2+(2*self.maxNumberofLinkDigits),pos=0),
        }
        tableFilteredDirectDef = {
            "fromCkmerLink": self._getTablesUint(self.numberOfKmers,0),
            "fromDirection": tables.StringCol(1,pos=1),
            "toCkmerLink": self._getTablesUint(self.numberOfKmers,2),
            "toDirection": tables.StringCol(1,pos=3),
            "distance": self._getTablesUint(self.maxDirectDistance,4),
            "number": tables.UInt16Col(pos=5),
        }
        tableGroupedCycleDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
            "minimumLength": tables.UInt64Col(pos=1),
            "number": tables.UInt64Col(pos=2),
        }
        tableGroupedReversalDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
            "minimumLength": tables.UInt64Col(pos=1),
            "number": tables.UInt64Col(pos=2),
        }
        tableDumpDirectErrorDef = {
            "sortKey": tables.StringCol(2+(2*self.maxNumberofLinkDigits),pos=0),
        }
        self.tableGroupedDirect = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedDirect", tableGroupedDirectDef, 
                                "Temporary to store grouped direct relations")
        self.tableGroupedDirectError = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedDirectError", tableGroupedDirectErrorDef, 
                                "Temporary to store grouped direct error relations")
        self.tableFilteredDirect = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "filteredDirect", tableFilteredDirectDef, 
                                "Temporary to store filtered direct relations")
        self.tableGroupedCycle = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedCycle", tableGroupedCycleDef, 
                                "Temporary to store grouped cycles")
        self.tableGroupedReversal = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedReversal", tableGroupedReversalDef, 
                                "Temporary to store grouped cycles")
        self.tableDumpDirectError = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root, 
                                "dumpDirectError", tableDumpDirectErrorDef, 
                                "Temporary to dump direct errors")
        def saveGroupedDirect(sortKey,fromCkmer,fromDirection,toCkmer,toDirection,distance,number):   
            groupedDirectRow = self.tableGroupedDirect.row            
            groupedDirectRow["sortKey"] = sortKey
            groupedDirectRow["fromCkmerLink"] = fromCkmer
            groupedDirectRow["fromDirection"] = fromDirection
            groupedDirectRow["toCkmerLink"] = toCkmer
            groupedDirectRow["toDirection"] = toDirection
            groupedDirectRow["distance"] = distance
            groupedDirectRow["number"] = number
            groupedDirectRow.append()    
        def saveGroupedDirectError(sortKey):   
            groupedDirectErrorRow = self.tableGroupedDirectError.row            
            groupedDirectErrorRow["sortKey"] = sortKey
            groupedDirectErrorRow.append()    
        def saveFilteredDirect(fromCkmer,fromDirection,toCkmer,toDirection,distance,number):   
            filteredDirectRow = self.tableFilteredDirect.row            
            filteredDirectRow["fromCkmerLink"] = fromCkmer
            filteredDirectRow["fromDirection"] = fromDirection
            filteredDirectRow["toCkmerLink"] = toCkmer
            filteredDirectRow["toDirection"] = toDirection
            filteredDirectRow["distance"] = distance
            filteredDirectRow["number"] = number
            filteredDirectRow.append()    
        def saveGroupedCycle(ckmerLink,minimumLength,number):  
            self.maxCycleLength = max(self.maxCycleLength,minimumLength)
            self.maxCycleNumber = max(self.maxCycleNumber,number)
            groupedCycleRow = self.tableGroupedCycle.row            
            groupedCycleRow["ckmerLink"] = ckmerLink
            groupedCycleRow["minimumLength"] = minimumLength
            groupedCycleRow["number"] = number
            groupedCycleRow.append()
        def saveGroupedReversal(ckmerLink,minimumLength,number): 
            self.maxReversalLength = max(self.maxReversalLength,minimumLength)
            self.maxReversalNumber = max(self.maxReversalNumber,number)
            groupedReversalRow = self.tableGroupedReversal.row            
            groupedReversalRow["ckmerLink"] = ckmerLink
            groupedReversalRow["minimumLength"] = minimumLength
            groupedReversalRow["number"] = number
            groupedReversalRow.append()
        def saveDumpDirectErrors(removeList):
            dumpDirectErrorRow = self.tableDumpDirectError.row
            for sortKey in removeList:
                dumpDirectErrorRow["sortKey"] = sortKey
                dumpDirectErrorRow.append()
                #keep symmetry
                n=round(len(sortKey)/2)
                dumpDirectErrorRow["sortKey"] = sortKey[n:]+sortKey[:n]
                dumpDirectErrorRow.append()
            
        def checkIndexCounter(indexCounter):
            #initialise
            removeList = set()            
            
            #distance based check and compute trunks
            ckmerInfo = self.h5file.get("/split/ckmer")
            trunkFreqs = {}
            distanceFreqs = {}
            for toCkmerLink in indexCounter.keys():
                for toDirection in indexCounter[toCkmerLink].keys():
                    toCkmer = ckmerInfo[toCkmerLink][0].decode()
                    number = indexCounter[toCkmerLink][toDirection][0]
                    distance = int(np.median(indexCounter[toCkmerLink][toDirection][1]))
                    if toDirection.decode()=="l":
                        trunk=toCkmer[:-1]
                    else:
                        trunk=haplotyping.General.reverse_complement(toCkmer)[:-1]
                    #update with median distance and trunk
                    indexCounter[toCkmerLink][toDirection][1] = distance
                    indexCounter[toCkmerLink][toDirection][3] = trunk
                    if not trunk in trunkFreqs.keys():
                        trunkFreqs[trunk]=number
                    else:
                        trunkFreqs[trunk]+=number
                    if not distance in distanceFreqs.keys():
                        distanceFreqs[distance]=number
                    else:
                        distanceFreqs[distance]+=number
            if not ((len(trunkFreqs)>1) or (len(distanceFreqs)>1)):
                #(probably) no problem
                return (indexCounter, removeList, False)
            
            #try restricting to the most frequent trunk
            mostFrequentTrunks = [(k,v) for k, v in trunkFreqs.items() if v==max(trunkFreqs.values())]
            if len(mostFrequentTrunks)==1:
                mostFrequentTrunk = mostFrequentTrunks[0][0]
                newIndexCounter = {}
                for toCkmerLink in indexCounter.keys():
                    for toDirection in indexCounter[toCkmerLink].keys():
                        number = indexCounter[toCkmerLink][toDirection][0]
                        distance = indexCounter[toCkmerLink][toDirection][1]
                        sortKey = indexCounter[toCkmerLink][toDirection][2]
                        trunk = indexCounter[toCkmerLink][toDirection][3]
                        if trunk==mostFrequentTrunk:
                            if not toCkmerLink in newIndexCounter.keys():
                                newIndexCounter[toCkmerLink] = {}
                            newIndexCounter[toCkmerLink][toDirection] = [number,distance,sortKey,None]
                        else:
                            removeList.add(sortKey)
                return (newIndexCounter,removeList, False)
            
            #problems, but not directly fixable
            return (indexCounter, removeList, True)

        #create grouped storage direct relations
        self.tableDumpDirect.flush()  
        numberOfDirect=self.tableDumpDirect.shape[0]
        self._logger.debug("group "+str(numberOfDirect)+" direct relations")
        self.tableDumpDirect.cols.sortKey.create_csindex()
        self.tableDumpDirect.flush()  
        previousSortKey=None
        fromLink=None
        fromDirection=None
        indexCounter={}
        self.maxDirectNumber = 0                
        for row in self.tableDumpDirect.itersorted("sortKey"):
            if not row["sortKey"]==previousSortKey:
                if not ((row["fromCkmerLink"]==fromLink) and 
                        (row["fromDirection"]==fromDirection)):
                    if previousSortKey:
                        (indexCounter,removeList,detectedProblem) = checkIndexCounter(indexCounter)
                        saveDumpDirectErrors(removeList)
                        for toCkmerLink in indexCounter.keys():
                            for toDirection in indexCounter[toCkmerLink].keys():
                                number = indexCounter[toCkmerLink][toDirection][0]
                                distance = indexCounter[toCkmerLink][toDirection][1]
                                sortKey = indexCounter[toCkmerLink][toDirection][2]
                                saveGroupedDirect(sortKey,fromLink,fromDirection,
                                                  toCkmerLink,toDirection,distance,number) 
                                self.maxDirectNumber = max(self.maxDirectNumber,number)
                    indexCounter={}
                    fromLink=row["fromCkmerLink"]
                    fromDirection=row["fromDirection"]
                previousSortKey=row["sortKey"]                    
                if not row["toCkmerLink"] in indexCounter.keys():
                    indexCounter[row["toCkmerLink"]]={}
                if not row["toDirection"] in indexCounter[row["toCkmerLink"]].keys():
                    indexCounter[row["toCkmerLink"]][row["toDirection"]] = [0,[],row["sortKey"],None]
            indexCounter[row["toCkmerLink"]][row["toDirection"]][0]+=1
            indexCounter[row["toCkmerLink"]][row["toDirection"]][1].append(row["distance"])
        if previousSortKey:
            (indexCounter,removeList,detectedProblem) = checkIndexCounter(indexCounter)
            saveDumpDirectErrors(removeList)
            for toCkmerLink in indexCounter.keys():
                for toDirection in indexCounter[toCkmerLink].keys():
                    number = indexCounter[toCkmerLink][toDirection][0]
                    distance = indexCounter[toCkmerLink][toDirection][1]
                    sortKey = indexCounter[toCkmerLink][toDirection][2]
                    saveGroupedDirect(sortKey,fromLink,fromDirection,
                                      toCkmerLink,toDirection,distance,number) 
                    self.maxDirectNumber = max(self.maxDirectNumber,number)  
                        
        #create grouped storage direct relations that are wrong
        self.tableDumpDirectError.flush() 
        numberOfDirectError=self.tableDumpDirectError.shape[0]
        self._logger.debug("group "+str(numberOfDirectError)+" direct error relations")
        self.tableDumpDirectError.cols.sortKey.create_csindex()
        self.tableDumpDirectError.flush() 
        previousSortKey=None
        for row in self.tableDumpDirectError.itersorted("sortKey"):
            if not row["sortKey"]==previousSortKey:
                if previousSortKey:
                    saveGroupedDirectError(previousSortKey.decode())
                previousSortKey=row["sortKey"]
        if previousSortKey:
            saveGroupedDirectError(previousSortKey.decode())
        self.tableGroupedDirectError.flush() 
        
        #create filtered storage direct relations
        self.tableGroupedDirect.flush() 
        numberOfDirect=self.tableGroupedDirect.shape[0]
        numberOfDirectError=self.tableGroupedDirectError.shape[0]
        self.tableGroupedDirect.cols.sortKey.create_csindex()
        self.tableGroupedDirect.flush() 
        self.tableGroupedDirectError.cols.sortKey.create_csindex()
        self.tableGroupedDirectError.flush()
        errorIterator = self.tableGroupedDirectError.itersorted("sortKey")
        self._logger.debug("filter "+str(numberOfDirect)+" by checking "+str(numberOfDirectError)+" errors")
        errorRow = None
        errorIteratorFinished = False
        for row in self.tableGroupedDirect.itersorted("sortKey"):
            if not errorIteratorFinished:
                if (errorRow==None) or (errorRow and errorRow[0]<row[0]):
                    errorIteratorFinished = True
                    for errorRow in errorIterator:
                        if errorRow[0]>=row[0]:
                            errorIteratorFinished = False
                            break  
            if errorIteratorFinished or not (errorRow[0]==row[0]):
                saveFilteredDirect(row[1],row[2],row[3],row[4],row[5],row[6])
        self.tableFilteredDirect.flush()
                
               
        #create grouped storage cycles
        self.tableDumpCycle.flush()  
        numberOfCycles=self.tableDumpCycle.shape[0]
        self._logger.debug("group "+str(numberOfCycles)+" cycles")
        self.tableDumpCycle.cols.ckmerLink.create_csindex()
        self.tableDumpCycle.flush()  
        previousLink=None
        previousLength=None
        previousNumber=None    
        for row in self.tableDumpCycle.itersorted("ckmerLink"):
            if not previousLink==row["ckmerLink"]:
                if not previousLink==None:
                    saveGroupedCycle(previousLink,previousLength,previousNumber)
                previousLink=row["ckmerLink"]
                previousLength=row["length"]
                previousNumber=0
            previousLength=min(previousLength,row["length"])
            previousNumber+=1
        if not previousLink==None:
            saveGroupedCycle(previousLink,previousLength,previousNumber)
        self.tableGroupedCycle.flush()  
         
        #create grouped storage reversals
        self.tableDumpReversal.flush()  
        numberOfReversals=self.tableDumpReversal.shape[0]
        self._logger.debug("group "+str(numberOfReversals)+" reversals")
        self.tableDumpReversal.cols.ckmerLink.create_csindex()
        self.tableDumpReversal.flush()  
        previousLink=None
        previousLength=None
        previousNumber=None              
        for row in self.tableDumpReversal.itersorted("ckmerLink"):
            if not previousLink==row["ckmerLink"]:
                if not previousLink==None:
                    saveGroupedReversal(previousLink,previousLength,previousNumber)
                previousLink=row["ckmerLink"]
                previousLength=row["length"]
                previousNumber=0
            previousLength=min(previousLength,row["length"])
            previousNumber+=1
        if not previousLink==None:
            saveGroupedReversal(previousLink,previousLength,previousNumber)
        self.tableGroupedReversal.flush()  
            
    def _groupConnected(self):
        #create temporary storage
        numberDumpConnected = self.tableDumpConnected.shape[0]
        tableGroupedConnectedDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
        }   
        tableGroupedConnectedIndexDef = {
            "hashKey": tables.StringCol(21,pos=0),
            "sizeConnected": self._getTablesUint(self.maxConnectedIndexSize,1),
            "lengthConnected": self._getTablesUint(self.maxConnectedIndexLength,2),
            "firstDumpConnectedLink": self._getTablesUint(numberDumpConnected,3),
            "firstGroupedConnectedLink": self._getTablesUint(numberDumpConnected,4),
            "number": self._getTablesUint(numberDumpConnected,5),
            "processed": tables.UInt8Col(pos=6),
        } 
        self.tableGroupedConnected = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedConnected", tableGroupedConnectedDef, 
                                "Temporary to store grouped connected")
        self.tableGroupedConnectedIndex = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedConnectedIndex", tableGroupedConnectedIndexDef, 
                                "Temporary to store grouped connected index")
        
        def saveGroupedConnectedIndex(hashKey,sizeConnected,lengthConnected,
                                      firstDumpConnectedLink,firstGroupedConnectedLink,number):
            #update number
            self.maxConnectedIndexNumber = max(self.maxConnectedIndexNumber,number)
            #store in index
            groupedConnectedIndexRow = self.tableGroupedConnectedIndex.row
            groupedConnectedIndexRow["hashKey"] = hashKey
            groupedConnectedIndexRow["sizeConnected"] = sizeConnected
            groupedConnectedIndexRow["lengthConnected"] = lengthConnected
            groupedConnectedIndexRow["firstDumpConnectedLink"] = firstDumpConnectedLink
            groupedConnectedIndexRow["firstGroupedConnectedLink"] = firstGroupedConnectedLink
            groupedConnectedIndexRow["number"] = number
            groupedConnectedIndexRow.append()
            #store connected
            groupedConnectedRow = self.tableGroupedConnected.row
            connected = self.tableDumpConnected[firstDumpConnectedLink:firstDumpConnectedLink+sizeConnected]
            for ckmerLink in connected:
                groupedConnectedRow["ckmerLink"] = ckmerLink[0]
                groupedConnectedRow.append()
        
        #create grouped storage connected index
        self.tableDumpConnectedIndex.flush()
        self.tableDumpConnected.flush()
        numberDumpConnectedIndex = self.tableDumpConnectedIndex.shape[0]
        numberDumpConnected = self.tableDumpConnected.shape[0]
        self._logger.debug("group "+str(numberDumpConnectedIndex)+" sets of "+str(numberDumpConnected)+" connected")
        self.tableDumpConnectedIndex.cols.hashKey.create_csindex()
        self.tableDumpConnectedIndex.flush()
        previousHashKey=None
        indexCounter={}        
        firstGroupedConnectedLink=0
        for row in self.tableDumpConnectedIndex.itersorted("hashKey"):
            if not row["hashKey"]==previousHashKey:
                if previousHashKey:
                    for sizeConnected in indexCounter.keys():
                        lengthConnected = int(np.median(indexCounter[sizeConnected][2]))
                        saveGroupedConnectedIndex(previousHashKey,sizeConnected,lengthConnected,
                                                  indexCounter[sizeConnected][0],
                                                  firstGroupedConnectedLink,
                                                  indexCounter[sizeConnected][1]) 
                        firstGroupedConnectedLink+=sizeConnected
                previousHashKey=row["hashKey"]                     
            if not row["sizeConnected"] in indexCounter.keys():
                indexCounter={row["sizeConnected"]:[row["firstDumpConnectedLink"],0,[]]}  
            indexCounter[row["sizeConnected"]][1]+=1
            indexCounter[row["sizeConnected"]][2].append(row["lengthConnected"])
        if previousHashKey:
            for sizeConnected in indexCounter.keys():
                lengthConnected = int(np.median(indexCounter[sizeConnected][2]))
                saveGroupedConnectedIndex(previousHashKey,sizeConnected,lengthConnected,
                                                  indexCounter[sizeConnected][0],
                                                  firstGroupedConnectedLink,
                                                  indexCounter[sizeConnected][1]) 
                firstGroupedConnectedLink+=sizeConnected    
        self.tableGroupedConnectedIndex.flush()
        self.tableGroupedConnected.flush()
                
        
        #create temporary storage
        numberGroupedConnectedIndex = self.tableGroupedConnectedIndex.shape[0]
        numberDumpConnectedPair = self.tableDumpConnectedPair.shape[0]
        tableGroupedConnectedPairDef = {
            "hashKey0": tables.StringCol(21,pos=0),
            "groupedConnectedIndexLink0": self._getTablesUint(numberGroupedConnectedIndex,1),
            "hashKey1": tables.StringCol(21,pos=2),
            "groupedConnectedIndexLink1": self._getTablesUint(numberGroupedConnectedIndex,3),
            "number": self._getTablesUint(numberDumpConnectedPair,4),
        } 
        self.tableGroupedConnectedPair = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "groupedConnectedPair", tableGroupedConnectedPairDef, 
                                "Temporary to store grouped connected pairs")
        
        def saveGroupedConnectedPair(hashKey0,hashKey1,number):
            #update number
            self.maxConnectedPairNumber = max(self.maxConnectedPairNumber,number)
            #store in index
            groupedConnectedPairRow = self.tableGroupedConnectedPair.row
            groupedConnectedPairRow["hashKey0"] = hashKey0
            groupedConnectedPairRow["groupedConnectedIndexLink0"] = 0
            groupedConnectedPairRow["hashKey1"] = hashKey1
            groupedConnectedPairRow["groupedConnectedIndexLink1"] = 0
            groupedConnectedPairRow["number"] = number
            groupedConnectedPairRow.append()            
        
        #create grouped storage connected pairs
        numberDumpConnectedPair = self.tableDumpConnectedPair.shape[0]
        self._logger.debug("group "+str(numberDumpConnectedPair)+" pairs of connected")
        self.tableDumpConnectedPair.cols.hashKey0.create_csindex()
        self.tableDumpConnectedPair.flush()
        previousHashKey0=None
        indexCounter={}  
        for row in self.tableDumpConnectedPair.itersorted("hashKey0"):
            if not row["hashKey0"]==previousHashKey0:
                if previousHashKey0:
                    for previousHashKey1 in indexCounter.keys():
                        saveGroupedConnectedPair(previousHashKey0,previousHashKey1,indexCounter[previousHashKey1])
                previousHashKey0=row["hashKey0"]   
            if not row["hashKey1"] in indexCounter.keys():
                indexCounter={row["hashKey1"]:0}
            indexCounter[row["hashKey1"]]+=1
        if previousHashKey0:
            for previousHashKey1 in indexCounter.keys():
                saveGroupedConnectedPair(previousHashKey0,previousHashKey1,indexCounter[previousHashKey1])
        #replace hashkeys by links
        self.tableGroupedConnectedPair.flush()
        self.tableGroupedConnectedIndex.cols.hashKey.create_csindex()
        self.tableGroupedConnectedIndex.flush()
        self.tableGroupedConnectedPair.cols.hashKey0.create_csindex()
        self.tableGroupedConnectedPair.cols.hashKey1.create_csindex()
        self.tableGroupedConnectedPair.flush()
        #add first link
        iteratorIndex = self.tableGroupedConnectedIndex.itersorted("hashKey")
        iteratorPair = self.tableGroupedConnectedPair.itersorted("hashKey0")
        currentHashKey = None
        currentRow = None
        for rowPair in iteratorPair:
            if rowPair["hashKey0"]==currentHashKey:
                rowPair["groupedConnectedIndexLink0"]=currentRow
                rowPair.update()
            else:
                for rowIndex in iteratorIndex:        
                    if rowIndex["hashKey"]==rowPair["hashKey0"]:
                        currentHashKey = rowIndex["hashKey"]
                        currentRow = rowIndex.nrow
                        rowPair["groupedConnectedIndexLink0"]=currentRow
                        rowPair.update()
                        break
        #add second link
        iteratorIndex = self.tableGroupedConnectedIndex.itersorted("hashKey")
        iteratorPair = self.tableGroupedConnectedPair.itersorted("hashKey1")
        currentHashKey = None
        currentRow = None
        for rowPair in iteratorPair:
            if rowPair["hashKey1"]==currentHashKey:
                rowPair["groupedConnectedIndexLink1"]=currentRow
                rowPair.update()
            else:
                for rowIndex in iteratorIndex:        
                    if rowIndex["hashKey"]==rowPair["hashKey1"]:
                        currentHashKey = rowIndex["hashKey"]
                        currentRow = rowIndex.nrow
                        rowPair["groupedConnectedIndexLink1"]=currentRow
                        rowPair.update()
                        break
        
        
    def _indexConnected(self):
        #create temporary storage
        self.tableGroupedConnectedIndex.flush()
        numberGroupedConnectedIndex = self.tableGroupedConnectedIndex.shape[0]
        tableCkmerConnectedIndexDef = {
            "firstCkmerConnectedIndexListLink": tables.UInt64Col(pos=0),
            "sizeCkmerConnectedIndexList": tables.UInt8Col(dflt=0,pos=1),
        }
        tableCkmerConnectedIndexListDef = {
            "groupedConnectedIndexLink": self._getTablesUint(numberGroupedConnectedIndex,0),
        }
        tableDumpCkmerConnectedIndexDef = {
            "ckmerLink": self._getTablesUint(self.numberOfKmers,0),
            "groupedConnectedIndexLink": self._getTablesUint(numberGroupedConnectedIndex,1),
        }
        self.tableCkmerConnectedIndex = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "ckmerConnectedIndex", tableCkmerConnectedIndexDef, 
                                "Temporary to store ckmer index on the ckmer connected index list",
                                expectedrows=self.numberOfKmers)
        self.tableCkmerConnectedIndexList = self.pytables_main_storage.create_table(self.pytables_main_storage.root, 
                                "ckmerConnectedIndexList", tableCkmerConnectedIndexListDef, 
                                "Temporary to store the ckmer connected index list")
        self.tableDumpCkmerConnectedIndex = self.pytables_dump_storage.create_table(self.pytables_dump_storage.root, 
                                "dumpCkmerConnectedIndex", tableDumpCkmerConnectedIndexDef, 
                                "Temporary to store the ckmer connected index list")
        
        #fill tableGroupedConnectedIndex
        self._logger.debug("dump k-mer index for "+str(numberGroupedConnectedIndex)+" sets of connected k-mers")
        dumpCkmerConnectedIndexRow = self.tableDumpCkmerConnectedIndex.row
        for id in range(numberGroupedConnectedIndex):
            row = self.tableGroupedConnectedIndex[id]
            connected = (self.tableDumpConnected[row["firstDumpConnectedLink"]:
                                                 row["firstDumpConnectedLink"]+row["sizeConnected"]])
            for i in range(len(connected)):
                dumpCkmerConnectedIndexRow["ckmerLink"]=connected[i][0]
                dumpCkmerConnectedIndexRow["groupedConnectedIndexLink"]=id
                dumpCkmerConnectedIndexRow.append()
            
        #fill tableCkmerConnectedIndex and tableCkmerConnectedIndexList
        self.tableDumpCkmerConnectedIndex.flush()
        self.tableDumpCkmerConnectedIndex.cols.ckmerLink.create_csindex()
        self.tableDumpCkmerConnectedIndex.flush()
        ckmerConnectedIndexRow = self.tableCkmerConnectedIndex.row
        ckmerConnectedIndexListRow = self.tableCkmerConnectedIndexList.row
        previousCkmer = None
        currentCkmer = -1
        previousGroupedConnectedIndexLinks=[]
        tableCkmerConnectedIndexListCounter=0
        
        def addCkmerConnectedIndex(currentCkmer, tableCkmerConnectedIndexListCounter, 
                                   ckmerLink, groupedConnectedIndexLinks):
            for i in range(currentCkmer+1,ckmerLink):
                ckmerConnectedIndexRow["firstCkmerConnectedIndexListLink"] = 0
                ckmerConnectedIndexRow["sizeCkmerConnectedIndexList"] = 0
                ckmerConnectedIndexRow.append()
            #only if linked to a limited set
            if len(groupedConnectedIndexLinks)<=255:
                ckmerConnectedIndexRow["firstCkmerConnectedIndexListLink"] = tableCkmerConnectedIndexListCounter
                ckmerConnectedIndexRow["sizeCkmerConnectedIndexList"] = len(groupedConnectedIndexLinks)
                ckmerConnectedIndexRow.append()
                for link in groupedConnectedIndexLinks:
                    ckmerConnectedIndexListRow["groupedConnectedIndexLink"]=link
                    ckmerConnectedIndexListRow.append()
                    tableCkmerConnectedIndexListCounter+=1
            else:
                ckmerConnectedIndexRow["firstCkmerConnectedIndexListLink"] = 0
                ckmerConnectedIndexRow["sizeCkmerConnectedIndexList"] = 0
                ckmerConnectedIndexRow.append()
            return (ckmerLink,tableCkmerConnectedIndexListCounter)

        numberDumpCkmerConnectedIndex = self.tableDumpCkmerConnectedIndex.shape[0]
        self._logger.debug("group "+str(numberDumpCkmerConnectedIndex)+" items from k-mer index")
        for row in self.tableDumpCkmerConnectedIndex.itersorted("ckmerLink"):
            if not previousCkmer==row["ckmerLink"]:
                if previousCkmer:
                    (currentCkmer,tableCkmerConnectedIndexListCounter) = addCkmerConnectedIndex(
                        currentCkmer, tableCkmerConnectedIndexListCounter, 
                        previousCkmer, set(previousGroupedConnectedIndexLinks))
                #reset
                previousCkmer=row["ckmerLink"]
                previousGroupedConnectedIndexLinks=[]
            previousGroupedConnectedIndexLinks.append(row["groupedConnectedIndexLink"])
        if previousCkmer:
            (currentCkmer,tableCkmerConnectedIndexListCounter) = addCkmerConnectedIndex(
                        currentCkmer, tableCkmerConnectedIndexListCounter, 
                        previousCkmer, previousGroupedConnectedIndexLinks)
        for i in range(currentCkmer+1,self.numberOfKmers):
            ckmerConnectedIndexRow["firstCkmerConnectedIndexListLink"] = 0
            ckmerConnectedIndexRow["sizeCkmerConnectedIndexList"] = 0
            ckmerConnectedIndexRow.append()           
               
    def _store(self):
        self._storeCkmer()
        self._storeDirect()
        self._storeConnected()
        self.h5file["/config/"].attrs["readLengthMinimum"]=self.readLengthMinimum
        self.h5file["/config/"].attrs["readLengthMaximum"]=self.readLengthMaximum
        self.h5file["/config/"].attrs["readPairedTotal"]=self.readPairedTotal
        self.h5file["/config/"].attrs["readUnpairedTotal"]=self.readUnpairedTotal
        self.h5file["/config/"].attrs["readTotal"]=self.readTotal
        self.h5file["/config/"].attrs["processReadsTime"]=int(np.ceil(self.processReadsTime))
        self.h5file["/config/"].attrs["maxDirectDistance"]=self.maxDirectDistance
        self.h5file["/config/"].attrs["maxDirectNumber"]=self.maxDirectNumber
        self.h5file["/config/"].attrs["maxConnectedIndexSize"]=self.maxConnectedIndexSize
        self.h5file["/config/"].attrs["maxConnectedIndexLength"]=self.maxConnectedIndexLength
        self.h5file["/config/"].attrs["maxConnectedIndexNumber"]=self.maxConnectedIndexNumber
        self.h5file["/config/"].attrs["maxConnectedPairNumber"]=self.maxConnectedPairNumber  
        self.h5file["/config/"].attrs["numberOfCycles"]=self.numberOfCycles
        self.h5file["/config/"].attrs["maxCycleLength"]=self.maxCycleLength
        self.h5file["/config/"].attrs["maxCycleNumber"]=self.maxCycleNumber
        self.h5file["/config/"].attrs["numberOfCycles"]=self.numberOfCycles
        self.h5file["/config/"].attrs["numberOfReversals"]=self.numberOfReversals
        self.h5file["/config/"].attrs["maxReversalLength"]=self.maxReversalLength
        self.h5file["/config/"].attrs["maxReversalNumber"]=self.maxReversalNumber
        
    def _storeCkmer(self):
        #create storage cycles
        self.tableGroupedCycle.flush()  
        self.numberOfCycles=self.tableGroupedCycle.shape[0]
        dtypeList = [("ckmerLink",self._getUint(self.numberOfKmers)),
                     ("minimumLength",self._getUint(self.maxCycleLength)),
                     ("number",self._getUint(self.maxCycleNumber))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/split/"].create_dataset("cycle",(self.numberOfCycles,),
                                          dtype=dt,chunks=None)
        self._logger.info("store "+str(self.numberOfCycles)+" cycles")
        ckmer = self.h5file.get("/split/ckmer")
        #add stored and grouped cycles to the final unchunked storage
        for i in range(0,self.numberOfCycles,self.stepSizeStorage):
            selection = [x[0] for x in self.tableGroupedCycle[i:i+self.stepSizeStorage]]
            values = [(x[2],) for x in self.tableGroupedCycle[i:i+self.stepSizeStorage]]
            ds[i:i+self.stepSizeStorage] = self.tableGroupedCycle[i:i+self.stepSizeStorage]    
            ckmer[selection,"cycle"] = values
        #create storage reversals
        self.tableGroupedReversal.flush()  
        self.numberOfReversals=self.tableGroupedReversal.shape[0]
        dtypeList = [("ckmerLink",self._getUint(self.numberOfKmers)),
                     ("minimumLength",self._getUint(self.maxReversalLength)),
                     ("number",self._getUint(self.maxReversalNumber))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/split/"].create_dataset("reversal",(self.numberOfReversals,),
                                          dtype=dt,chunks=None)
        self._logger.info("store "+str(self.numberOfReversals)+" reversals")
        #add stored and grouped reversals to the final unchunked storage
        for i in range(0,self.numberOfReversals,self.stepSizeStorage):
            selection = [x[0] for x in self.tableGroupedReversal[i:i+self.stepSizeStorage]]
            values = [(x[2],) for x in self.tableGroupedReversal[i:i+self.stepSizeStorage]]
            ds[i:i+self.stepSizeStorage] = self.tableGroupedReversal[i:i+self.stepSizeStorage]   
            ckmer[selection,"reversal"] = values
        
        
    def _storeDirect(self):
        #create storage  
        self.tableFilteredDirect.flush()  
        self.numberOfDirectRelations=self.tableFilteredDirect.shape[0]
        dtypeList = [("fromCkmerLink",self._getUint(self.numberOfKmers)),("fromDirection","S1"),
                     ("toCkmerLink",self._getUint(self.numberOfKmers)),("toDirection","S1"),
                     ("distance",self._getUint(self.maxDirectDistance)),
                     ("number",self._getUint(self.maxDirectNumber))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/relations/"].create_dataset("direct",(self.numberOfDirectRelations,),
                                          dtype=dt,chunks=None)
        self._logger.info("store "+str(self.numberOfDirectRelations)+" grouped direct relations")
        #add stored and grouped kmers to the final unchunked storage
        for i in range(0,self.numberOfDirectRelations,self.stepSizeStorage):
            ds[i:i+self.stepSizeStorage] = self.tableFilteredDirect[i:i+self.stepSizeStorage] 
        #update ckmers
        ckmer = self.h5file.get("/split/ckmer")
        self.tableFilteredDirect.cols.fromCkmerLink.create_csindex()
        self.tableFilteredDirect.flush()
        previousCkmer = None
        for row in self.tableFilteredDirect.itersorted("fromCkmerLink"):            
            if not row["fromCkmerLink"]==previousCkmer:
                if not previousCkmer==None:
                    ckmer[previousCkmer,"direct"] = ((distinctLeft,numberLeft,),(distinctRight,numberRight,),)
                previousCkmer=row["fromCkmerLink"] 
                numberLeft=0
                numberRight=0
                distinctLeft=0
                distinctRight=0
            if row["fromDirection"].decode()=="l":
                distinctLeft+=1
                numberLeft+=row["number"]
            elif row["fromDirection"].decode()=="r":
                distinctRight+=1    
                numberRight+=row["number"]
        if not previousCkmer==None:
            ckmer[previousCkmer,"direct"] = ((distinctLeft,numberLeft,),(distinctRight,numberRight,),)
                
            
    def _storeConnected(self):
        self.tableGroupedConnected.flush() 
        self.tableGroupedConnectedIndex.flush()
        self.tableCkmerConnectedIndex.flush()
        self.tableCkmerConnectedIndexList.flush()
        #store connected
        numberOfConnected=self.tableGroupedConnected.shape[0]
        dtypeList = [("ckmerLink",self._getUint(self.numberOfKmers))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/relations/"].create_dataset("connected",(numberOfConnected,),
                                          dtype=dt,chunks=None)
        self._logger.info("store "+str(numberOfConnected)+" connected k-mers")
        for i in range(0,numberOfConnected,self.stepSizeStorage):
            ds[i:i+self.stepSizeStorage] = [x for x in self.tableGroupedConnected[i:i+self.stepSizeStorage]]
        #store connectedIndex
        numberOfConnectedIndex=self.tableGroupedConnectedIndex.shape[0]
        dtypeList = [("sizeConnected",self._getUint(self.maxConnectedIndexSize)),
                     ("lengthConnected",self._getUint(self.maxConnectedIndexLength)),
                     ("connectedLink",self._getUint(numberOfConnected)),
                     ("number",self._getUint(self.maxConnectedIndexNumber))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/relations/"].create_dataset("connectedIndex",(numberOfConnectedIndex,),
                                          dtype=dt,chunks=None)
        self._logger.info("store "+str(numberOfConnectedIndex)+" sets of connected k-mers")
        for i in range(0,numberOfConnectedIndex,self.stepSizeStorage):
            ds[i:i+self.stepSizeStorage] = [(x[1],x[2],x[4],x[5],) for x in 
                                            self.tableGroupedConnectedIndex[i:i+self.stepSizeStorage]]    
        #store connectedPair
        numberOfConnectedPair=self.tableGroupedConnectedPair.shape[0]
        dtypeList = [("connectedIndexLink0",self._getUint(numberOfConnectedIndex)),
                     ("connectedIndexLink1",self._getUint(numberOfConnectedIndex)),
                     ("number",self._getUint(self.maxConnectedPairNumber))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/relations/"].create_dataset("connectedPair",(numberOfConnectedPair,),
                                          dtype=dt,chunks=None)
        self._logger.info("store "+str(numberOfConnectedPair)+" pairs of connected sets")
        for i in range(0,numberOfConnectedPair,self.stepSizeStorage):
            ds[i:i+self.stepSizeStorage] = [(x[1],x[3],x[4],) for x in 
                                            self.tableGroupedConnectedPair[i:i+self.stepSizeStorage]]    
        #store ckmer
        numberOfCkmerConnectedIndexList=self.tableCkmerConnectedIndexList.shape[0]
        dtypeList = [("ckmerConnectedLink",self._getUint(numberOfCkmerConnectedIndexList)),
                     ("sizeCkmerConnected","uint8")]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/relations/"].create_dataset("ckmer",(self.numberOfKmers,),
                                          dtype=dt,chunks=None)
        ckmer = self.h5file.get("/split/ckmer")
        for i in range(0,self.numberOfKmers,self.stepSizeStorage):
            values = [(x[1],) for x in self.tableCkmerConnectedIndex[i:i+self.stepSizeStorage]]
            ckmer[i:i+self.stepSizeStorage,"connected"] = values
            ds[i:i+self.stepSizeStorage] = [x for x in self.tableCkmerConnectedIndex[i:i+self.stepSizeStorage]]
        #store ckmerConnected
        numberOfGroupedConnectedIndex=self.tableGroupedConnectedIndex.shape[0]
        numberOfCkmerConnectedIndexList=self.tableCkmerConnectedIndexList.shape[0]
        dtypeList = [("connectedIndexLink",self._getUint(numberOfGroupedConnectedIndex))]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/relations/"].create_dataset("ckmerConnected",(numberOfCkmerConnectedIndexList,),
                                          dtype=dt,chunks=None)
        self._logger.info("store for "+str(self.numberOfKmers)+" k-mers the "+str(numberOfCkmerConnectedIndexList)+
                          " connected sets")
        for i in range(0,numberOfCkmerConnectedIndexList,self.stepSizeStorage):
            ds[i:i+self.stepSizeStorage] = [x for x in self.tableCkmerConnectedIndexList[i:i+self.stepSizeStorage]]
        #update ckmer
         
        
        
    def _hashConnected(self,connected):
        def computeHash(connected,i=0):
            #include escape for symmetry
            if ((2*i)+1)>=len(connected) or (connected[i]<connected[(-1-i)]):
                tupleConnected=tuple(connected)
                h = hash(tupleConnected)
                h += sys.maxsize + 1
                return str(h).zfill(21)
            elif connected[i]>connected[(-1-i)]:
                tupleConnected=tuple(connected)
                h = hash(tupleConnected[::-1])
                h += sys.maxsize + 1
                return str(h).zfill(21)
            else:
                return computeHash(connected,i+1)
        return computeHash(connected)
    
    def _computeProblemStartPositions(self,sequence):
        problemStartPositions = []
        problemPattern = re.compile(r"[ACGT][^ACGT]+")
        for m in re.finditer(problemPattern, sequence):
            problemStartPositions.append(m.span()[0]+1)
        return problemStartPositions
    
    def _computeMatchesList(self,sequence):
        def saveDumpCycle(link,length): 
            #save data
            cycleRow = self.tableDumpCycle.row    
            cycleRow["ckmerLink"] = link
            cycleRow["length"] = length        
            cycleRow.append() 
        def saveDumpReversal(link,length): 
            #save data
            reversalRow = self.tableDumpReversal.row    
            reversalRow["ckmerLink"] = link
            reversalRow["length"] = length        
            reversalRow.append() 
            
        #register history to detect (local) cycles and reversals
        history = {}
        matchesList = []
        problems = self._computeProblemStartPositions(sequence)
        relevantProblem = None if len(problems)==0 else problems[0]
        matches = []
        for end_index, (link,orientation,splitLocation) in self.automatonSplits.iter(sequence):
            pos=1+end_index-self.k
            if link in history.keys():
                if history[link][0]==orientation:
                    saveDumpCycle(link,1+pos-history[link][1])
                else:
                    saveDumpReversal(link,1+pos-history[link][1])
            history[link]=[orientation,pos]
            #check if a problem did occur between last match and current
            if relevantProblem and len(matches)>0 and pos>relevantProblem:
                matchesList.append(matches)
                matches=[]                            
            matches.append([pos,link,orientation,splitLocation])
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
    
    def _processConnected(self,matches,storeSingleConnected):    
        boundary = 0 if storeSingleConnected else 1
        if len(matches)>boundary:   
            startConnected=len(matches)  
            endConnected=0
            for i in range(len(matches)):            
                (pos,link,orientation,splitLocation)=matches[i]
                if not (splitLocation=="r"):
                    startConnected=i
                    break
            for i in range(len(matches))[::-1]:            
                (pos,link,orientation,splitLocation)=matches[i]
                if not (splitLocation=="l"):
                    endConnected=i
                    break  
            if (endConnected-startConnected)>boundary:
                connected = [x[1] for x in matches[startConnected:endConnected+1]]
                sizeConnected=len(connected)
                lengthConnected=1+matches[endConnected][0]-matches[startConnected][0]
                self.maxConnectedIndexSize=max(self.maxConnectedIndexSize,sizeConnected)
                self.maxConnectedIndexLength=max(self.maxConnectedIndexLength,lengthConnected)
                connectedHashKey = self._hashConnected(connected)
                connectedIndexRow = self.tableDumpConnectedIndex.row
                connectedIndexRow["hashKey"] = connectedHashKey
                connectedIndexRow["sizeConnected"] = sizeConnected
                connectedIndexRow["lengthConnected"] = lengthConnected
                connectedIndexRow["firstDumpConnectedLink"] = self.numberTableDumpConnected
                connectedIndexRow.append()
                #store connected
                connectedRow = self.tableDumpConnected.row                        
                for i in range(len(connected)):
                    connectedRow["ckmerLink"] = connected[i]
                    connectedRow.append()
                    self.numberTableDumpConnected+=1    
                return connectedHashKey
        return None
        
    def _processMatches(self,matches):
        def saveDumpDirect(link1,direction1,link2,direction2,distance): 
            #update maximum
            self.maxDirectDistance = max(self.maxDirectDistance,distance)
            #save data
            directRow = self.tableDumpDirect.row    
            #todo: this can be done more space efficient
            directRow["sortKey"] = (str(link1).zfill(self.maxNumberofLinkDigits)+direction1+
                                  str(link2).zfill(self.maxNumberofLinkDigits)+direction2)
            directRow["fromCkmerLink"] = link1
            directRow["fromDirection"] = direction1
            directRow["toCkmerLink"] = link2
            directRow["toDirection"] = direction2
            directRow["distance"] = distance        
            directRow.append()    
    
        for i in range(len(matches)):
            (pos1,link1,orientation1,splitLocation1)=matches[i]
            direction1="l" if orientation1=="r" else "r"
            #loop over right neighbours
            for j in range(i+1,len(matches)):
                (pos2,link2,orientation2,splitLocation2)=matches[j]
                direction2="r" if orientation2=="r" else "l"
                #store direct path relation
                if j==i+1: 
                    saveDumpDirect(link1,direction1,link2,direction2,pos2-pos1)
                    saveDumpDirect(link2,direction2,link1,direction1,pos2-pos1)


    def _getAutomaton(self):        
        if self.debug:
            filenameAutomaton = self.filenameBase+"_tmp_relations_automaton.data"
            if os.path.exists(filenameAutomaton):
                with open(filenameAutomaton, "rb") as f:
                    self.automatonSplits = pickle.load(f)
                    self._logger.warning("automaton loaded from previous run")
                    return
        else:
            filenameAutomaton = None
        
        self._logger.debug("create automaton for "+str(self.numberOfKmers)+" k-mers")
        self.automatonSplits = ahocorasick.Automaton()
        kmers = self.h5file["/split/ckmer"]
        for i in range(self.numberOfKmers):
            kmer = kmers[i][0].decode()
            kmerType = kmers[i][1].decode()
            rkmer = haplotyping.General.reverse_complement(kmer)
            rkmerType = "l" if kmerType=="r" else ("r" if kmerType=="l" else "b")
            self.automatonSplits.add_word(kmer,(i,"c",kmerType))
            #if canonical k-mer equals reverse complement, don't insert for rc
            if not kmer==rkmer:
                self.automatonSplits.add_word(rkmer,(i,"r",rkmerType))
        self.automatonSplits.make_automaton()
        if self.debug:
            with open(filenameAutomaton, "wb") as f:
                pickle.dump(self.automatonSplits, f)

                
    def _getTablesUint(self, maximumValue, position):
        if maximumValue<=np.iinfo(np.uint8).max:
            return tables.UInt8Col(pos=position)
        elif maximumValue<=np.iinfo(np.uint16).max:
            return tables.UInt16Col(pos=position)
        elif maximumValue<=np.iinfo(np.uint32).max:
            return tables.UInt32Col(pos=position)
        else:
            return tables.UInt64Col(pos=position)
        
    def _getUint(self, maximumValue):
        if maximumValue<=np.iinfo(np.uint8).max:
            return "uint8"
        elif maximumValue<=np.iinfo(np.uint16).max:
            return "uint16"
        elif maximumValue<=np.iinfo(np.uint32).max:
            return "uint32"
        else:
            return "uint64"
