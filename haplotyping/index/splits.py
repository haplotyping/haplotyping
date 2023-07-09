import logging, h5py, tables, gzip, time
import os, sys, shutil, psutil, numpy as np
import re, haplotyping, ahocorasick, pickle, signal
import haplotyping.index.database
import multiprocessing as mp
from queue import Empty

class Splits:
    
    """
    Internal use, get splitting k-mers from index
    """
    
    stepSizeStorage = 1000000
    
    def __init__(self, sortedIndexFile: str, h5file, filenameBase, debug=False, keepTemporaryFiles=False):
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("parse k-mers to find right splitting k-mers")
        
        #set variables
        self.k = h5file["/config"].attrs["k"]
        self.automatonKmerSize = h5file["/config"].attrs["automatonKmerSize"]
        self.minimumFrequency = h5file["/config"].attrs["minimumCanonicalSplitFrequency"]
        self.h5file = h5file
        self.maximumNumber = 0
        self.debug = debug
        self.keepTemporaryFiles = keepTemporaryFiles
        
        
        #check existence group
        if not "/split" in h5file:
            h5file.create_group("/split")
        if not "/histogram" in h5file:
            h5file.create_group("/histogram")
        
        #create splitting k-mer dataset
        if "/split/ckmer" in h5file:
            self._logger.warning("canonical splitting k-mer dataset already in hdf5 storage")
        elif "/histogram/kmer" in h5file:
            self._logger.warning("histogram k-mer dataset already in hdf5 storage")
        elif "/histogram/ckmer" in h5file:
            self._logger.warning("histogram splitting k-mer dataset already in hdf5 storage")
        elif "/histogram/base" in h5file:
            self._logger.warning("histogram splitting k-mer bases dataset already in hdf5 storage")
        else:
        
            try:                                                
                pytablesFile = filenameBase+"_tmp_split.h5"
                if os.path.exists(pytablesFile):
                    os.remove(pytablesFile)
                self._logger.debug("store temporary in "+pytablesFile)    
                
                #create datasets
                with tables.open_file(pytablesFile, mode="w", 
                                      title="Temporary storage") as pytablesStorage:
                    #get k-mers
                    self._parseIndex(sortedIndexFile, pytablesStorage)
                    self._sort(pytablesStorage)
                    self._store(pytablesStorage)                    
                    #flush
                    self.h5file.flush()                    
            except Exception as e:
                self._logger.error("problem occurred while constructing splits: "+str(e))
            finally:
                try:
                    if not self.keepTemporaryFiles:
                        os.remove(pytablesFile)
                except:
                    self._logger.error("problem removing "+pytablesFile)    

#----------------
# Main functions
#----------------            
                    
    def _parseIndex(self, filename: str, pytablesStorage):
        
        
        def close_queue(queue_entry):
            time.sleep(0.1)
            queue_entry.close()
            queue_entry.join_thread()
            
        
        #create temporary storage
        tableCkmerDef = {
            "ckmer": tables.StringCol(self.k,pos=0),
            "type": tables.StringCol(1,pos=1),
            "number": tables.UInt64Col(pos=2),
        }
        tableDumpKmers = pytablesStorage.create_table(pytablesStorage.root,
                                    "dumpCkmer",tableCkmerDef, "Temporary to store dump canonical k-mers")
        
        #multiprocessing
        with mp.Manager() as manager:
            shutdown_event = mp.Event()
            qsize = 10000
            queue_splits = mp.Queue(qsize)
            totalNumberOfKmers = mp.Value("L",0)
            totalSumOfKmerFrequencies = mp.Value("L",0)
            minimumAllKmerFrequencies = mp.Value("L",2**32)
            maximumAllKmerFrequencies = mp.Value("L",0)
            histogramKmer = manager.dict()
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            process_list = mp.Process(target=haplotyping.index.splits.Splits.workerList, 
                                      args=(shutdown_event,filename,self.k,self.minimumFrequency,
                                            totalNumberOfKmers,totalSumOfKmerFrequencies,
                                            minimumAllKmerFrequencies,maximumAllKmerFrequencies,histogramKmer,
                                            queue_splits, ))
            signal.signal(signal.SIGINT, original_sigint_handler)

            try:
                #start
                process_list.start()
                #init
                rightSplitBases = 0
                rightSplitKmers = 0
                self.frequencyHistogram = {"kmer": {}, "ckmer": {}, "base": {}}
                #process queue
                while True:
                    try:
                        item = queue_splits.get(block=True, timeout=1)
                        if item==None:
                            break
                        else:
                            rightSplitBases+=1
                            for key,value in item[1].items():
                                kmer = item[0]+key
                                ckmer = haplotyping.General.canonical(kmer)
                                #assume right splitting if k-mer equals rc
                                #todo: can this (in theory) cause problematic situations?
                                ckmerType = "r" if kmer==ckmer else "l"
                                ckmerRow = tableDumpKmers.row            
                                ckmerRow["ckmer"] = ckmer
                                ckmerRow["type"] = ckmerType
                                ckmerRow["number"] = value
                                ckmerRow.append()
                                rightSplitKmers+=1 
                                self.frequencyHistogram["ckmer"][value] = (
                                    self.frequencyHistogram["ckmer"].get(value,0)+1)
                            if (rightSplitBases%1000000)==0:
                                self._logger.debug("processed {} bases and {} k-mers".format(
                                    rightSplitBases,rightSplitKmers))
                    except Empty:
                        continue
                tableDumpKmers.flush()
                self._logger.debug("found {} rightSplitBases and {} rightSplitKmers".format(
                    rightSplitBases,rightSplitKmers))
                self.h5file["/config/"].attrs["numberRightSplitBases"]=rightSplitBases
                self.h5file["/config/"].attrs["numberRightSplitKmers"]=rightSplitKmers

            except KeyboardInterrupt:
                self._logger.debug("caught KeyboardInterrupt")
                #shutdown directly
                shutdown_event.set()
                self._logger.debug("forced exit")
                sys.exit()
            finally:
                #close queus workers
                close_queue(queue_splits)
                #shutdown
                shutdown_event.set()
                #join
                process_list.join() 

            #store stats
            self.frequencyHistogram["kmer"] = dict(histogramKmer)
            self.h5file["/config/"].attrs["numberKmers"]=totalNumberOfKmers.value
            self.h5file["/config/"].attrs["totalKmerFrequencies"]=totalSumOfKmerFrequencies.value
            self.h5file["/config/"].attrs["minimumKmerFrequencies"]=minimumAllKmerFrequencies.value
            self.h5file["/config/"].attrs["maximumKmerFrequencies"]=maximumAllKmerFrequencies.value
            
    def workerList(shutdown_event,filename,k,minimumFrequency,totalNumberOfKmers,totalSumOfKmerFrequencies,
                                        minimumAllKmerFrequencies,maximumAllKmerFrequencies,histogramKmer,
                                        queue_splits):
        logger = logging.getLogger("{}.worker.list".format(__name__))
        try:
            with gzip.open(filename, "rt") as f: 
                logger.debug("parse sorted list to detect right splitting k-mers")
                previousBase = ""
                previousBranch = ""
                previousKmer = ""
                previousNumber = 0                 
                #administration
                rightSplitBases = 0
                rightSplitKmers = 0
                stored={}
                errorNumber = 0
                _totalNumberOfKmers=0
                _totalSumOfKmerFrequencies=0
                _minimumAllKmerFrequencies=2**32
                _maximumAllKmerFrequencies=0
                _histogramKmer = {}
                while (row := f.readline()):
                    line = row.strip().split("\t") 
                    currentBase = line[0][:-1]
                    currentBranch = line[0][-1]
                    currentKmer = line[0]
                    _totalNumberOfKmers+=1
                    if (_totalNumberOfKmers%10000000)==0 and (_totalNumberOfKmers>0):
                        logger.debug("processed {} k-mers".format(_totalNumberOfKmers))
                    currentNumber = int(line[1])
                    _totalSumOfKmerFrequencies += currentNumber
                    _minimumAllKmerFrequencies = min(_minimumAllKmerFrequencies,currentNumber)
                    _maximumAllKmerFrequencies = max(_maximumAllKmerFrequencies,currentNumber)
                    _histogramKmer[currentNumber] = (_histogramKmer.get(currentNumber,0)+1)
                    if currentNumber < minimumFrequency:
                        pass
                    else:
                        if currentBase==previousBase:
                            stored[previousBranch]=previousNumber
                            stored[currentBranch]=currentNumber
                        else:
                            if len(stored)>1:
                                queue_splits.put((previousBase,stored,))
                                stored={}                                        
                        previousBase = currentBase
                        previousBranch = currentBranch
                        previousKmer = currentKmer
                        previousNumber = currentNumber
                if len(stored)>1:
                    queue_splits.put((previousBase,stored,)) 
                #store stats
                totalNumberOfKmers.value=_totalNumberOfKmers
                totalSumOfKmerFrequencies.value=_totalSumOfKmerFrequencies
                minimumAllKmerFrequencies.value=_minimumAllKmerFrequencies
                maximumAllKmerFrequencies.value=_maximumAllKmerFrequencies
                for key in _histogramKmer.keys():
                    histogramKmer[key] = _histogramKmer[key]
                #close queue
                queue_splits.put(None)
                #warning
                if errorNumber>0:
                    logger.warning("skipped {} items in sorted list".format(errorNumber))
                #stats
                minimumAllKmerFrequencies.value=min(minimumAllKmerFrequencies.value,maximumAllKmerFrequencies.value)
                logger.debug("checked {} k-mers with {} total frequency".format(
                    totalNumberOfKmers.value,totalSumOfKmerFrequencies.value))
                logger.debug("frequency k-mers between {} and {}".format(
                    minimumAllKmerFrequencies.value,maximumAllKmerFrequencies.value))
        except Exception as ex:
            logger.error("problem with sorted list: "+str(ex))
    

    def _sort(self, pytablesStorage):                        
        
        def saveToSortedCkmerStorage(ckmer,ckmerType,number,numberOfCkmers):
            ckmerRow = tableSortedKmers.row            
            ckmerRow["ckmer"] = ckmer
            ckmerRow["type"] = ckmerType
            ckmerRow["number"] = number
            ckmerRow.append()
            if ckmerType=="r" or ckmerType=="b":
                saveToDumpBaseStorage(ckmer[:-1],ckmer[-1],number,numberOfCkmers)
            if ckmerType=="l" or ckmerType=="b":
                rkmer=haplotyping.General.reverse_complement(ckmer)
                saveToDumpBaseStorage(rkmer[:-1],rkmer[-1],number,numberOfCkmers)
            return numberOfCkmers+1
        
        def saveToSortedBaseStorage(base,branches,number,numberOfBases):
            baseRow = tableSortedBases.row    
            baseRow["base"] = base
            baseRow["number"] = number
            self.frequencyHistogram["base"][number] = (
                                self.frequencyHistogram["base"].get(number,0)+1)
            for letter in branches.keys():
                baseRow["branches/"+str(letter)+"/number"] = branches[letter][0]
                baseRow["branches/"+str(letter)+"/ckmerLink"] = branches[letter][1]
                ckmerRow = tableSortedKmers[branches[letter][1]]
                ckmer = ckmerRow["ckmer"].decode()
                rkmer=haplotyping.General.reverse_complement(ckmer)
                if ckmer == base+letter:
                    baseReferences[branches[letter][1]][1] = numberOfBases
                if rkmer == base+letter:
                    baseReferences[branches[letter][1]][0] = numberOfBases
            baseRow.append()
            return numberOfBases+1
            
        def saveToDumpBaseStorage(base,branch,number,ckmerLink):   
            baseRow = tableDumpRightSplitBases.row
            baseRow["base"] = base
            baseRow["branch"] = branch
            baseRow["number"] = number
            baseRow["ckmerLink"] = ckmerLink
            baseRow.append()
            
        tableDumpKmers = pytablesStorage.root.dumpCkmer
        
        #create sorted and grouped storage
        tableDumpKmers.flush()                
        numberOfSplittingKmers=tableDumpKmers.shape[0]
        if numberOfSplittingKmers>0:
            
            #create temporary storage
            tableCkmerDef = {
                "ckmer": tables.StringCol(self.k,pos=0),
                "type": tables.StringCol(1,pos=1),
                "rightSplitBaseLink": {
                    "leftSplit": haplotyping.index.Database.getTablesUint(numberOfSplittingKmers,0),
                    "rightSplit": haplotyping.index.Database.getTablesUint(numberOfSplittingKmers,1),
                },
                "number": tables.UInt32Col(pos=3,dflt=0),
            }
            tableSortedKmers = pytablesStorage.create_table(pytablesStorage.root,
                                        "sortedCkmer",tableCkmerDef, "Temporary to store sorted canonical k-mers",
                                        expectedrows=numberOfSplittingKmers)  
            
            #create dump table for bases
            tableDumpBaseDef = {
                "base": tables.StringCol(self.k-1,pos=0),
                "branch": tables.StringCol(1,pos=1),
                "number": tables.UInt32Col(pos=2,dflt=0),
                "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfSplittingKmers,3),
            }
            tableDumpRightSplitBases = pytablesStorage.create_table(pytablesStorage.root,
                                                                   "dumpRightSplitBase",tableDumpBaseDef, 
                                                                   "Temporary to store dump bases",
                                                                       expectedrows=numberOfSplittingKmers)
            
            #store sorted k-mers
            self._logger.debug("sort and group "+str(numberOfSplittingKmers)+" splitting k-mers")
            tableDumpKmers.cols.ckmer.create_csindex()
            tableDumpKmers.flush()
            self._logger.debug("created index to sort k-mers")
            previousCkmer=None
            previousType=None
            previousN=None
            self.maximumNumber = 0
            numberOfCkmers = 0
            for row in tableDumpKmers.itersorted("ckmer",checkCSI=True):
                currentCkmer=row["ckmer"]
                currentType=row["type"]
                currentN=row["number"]
                if previousCkmer and not previousCkmer==currentCkmer:
                    numberOfCkmers = saveToSortedCkmerStorage(previousCkmer.decode(),previousType,
                                                         previousN,numberOfCkmers)
                    self.maximumNumber=max(self.maximumNumber,previousN)
                    previousType=currentType.decode()
                elif previousCkmer and not previousType==currentType:
                    previousType="b"
                else:
                    previousType=currentType.decode()
                previousCkmer=currentCkmer
                previousN=currentN
                if (numberOfCkmers%1000000)==0 and (numberOfCkmers>0):
                    self._logger.debug("processed {} k-mers".format(numberOfCkmers))
            if previousCkmer:        
                numberOfCkmers = saveToSortedCkmerStorage(previousCkmer.decode(),previousType,
                                                     previousN,numberOfCkmers) 
                self.maximumNumber=max(self.maximumNumber,previousN)
            self._logger.debug("in total processed {} k-mers".format(numberOfCkmers))
            
            tableDumpRightSplitBases.flush()                
            numberOfBases=tableDumpRightSplitBases.shape[0]
                
            #store sorted bases
            tableSortedKmers.flush()
            numberOfKmers=tableSortedKmers.shape[0]
            tableBaseDef = {
                "base": tables.StringCol(self.k-1,pos=0),
                "number": tables.UInt32Col(pos=1,dflt=0),    
                "branches": {},                            
            }
            for letter in haplotyping.index.Database.letters:
                tableBaseDef["branches"][letter] = {
                    "number": tables.UInt32Col(pos=0,dflt=0),
                    "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
                }
            tableSortedBases = pytablesStorage.create_table(pytablesStorage.root,
                                                                   "sortedBase",tableBaseDef, 
                                                                   "Temporary to store dump bases",
                                                                       expectedrows=numberOfKmers)
            
            #store sorted bases
            self._logger.debug("sort and group {} bases".format(numberOfBases))
            tableDumpRightSplitBases.cols.base.create_csindex()
            tableDumpRightSplitBases.flush()
            self._logger.debug("created index to sort bases")
            previousBase=None
            previousBranches={}
            previousNumber=0
            numberOfBases = 0
            #reserve memory for links in ckmer table
            ckmerLinkType = np.dtype(haplotyping.index.Database.getUint(numberOfKmers)).type
            baseReferences = np.full((numberOfKmers,2), 0, dtype=ckmerLinkType)
            for row in tableDumpRightSplitBases.itersorted("base",checkCSI=True):
                currentBase=row[0]
                currentBranch=row[1]
                currentNumber=row[2]
                currentCkmerLink=row[3]
                if previousBase and not previousBase==currentBase:
                    numberOfBases = saveToSortedBaseStorage(previousBase.decode(),previousBranches,
                                                            previousNumber,numberOfBases)
                    previousBase=None
                    previousBranches={}
                    previousNumber=0
                    if (numberOfBases%1000000)==0 and (numberOfBases>0):
                        self._logger.debug("processed {} bases".format(numberOfBases))
                previousBase=currentBase
                previousBranches[currentBranch.decode()] = (currentNumber, currentCkmerLink,)
                previousNumber+=currentNumber                    
            if previousBase:
                    numberOfBases = saveToSortedBaseStorage(previousBase.decode(),previousBranches,
                                                            previousNumber,numberOfBases)  
            self._logger.debug("in total processed {} bases".format(numberOfBases))
            #update ckmers
            tableSortedKmers.modify_columns(columns=baseReferences,
                             names=["rightSplitBaseLink/leftSplit","rightSplitBaseLink/rightSplit"])   
                
        else:
            self._logger.warning("no splitting k-mers to sort and group")
            
    def _store(self, pytablesStorage):
        canonicalSplitKmers = 0
        canonicalSplitKmersLeft = 0
        canonicalSplitKmersRight = 0
        canonicalSplitKmersBoth = 0
        tableSortedKmers = pytablesStorage.root.sortedCkmer
        tableSortedBases = pytablesStorage.root.sortedBase
        #create storage    
        tableSortedKmers.flush()
        tableSortedBases.flush()
        numberOfKmers=tableSortedKmers.shape[0]
        numberOfBases=tableSortedBases.shape[0]
        # CKMER STORAGE - don't make the structure unnecessary big
        dtypeCkmerList=[("ckmer","S"+str(self.k)),
                   ("type","S1"),
                   ("number",haplotyping.index.Database.getUint(self.maximumNumber)),
                   ("rightSplitBaseLink",[("leftSplit",haplotyping.index.Database.getUint(numberOfBases)),
                                          ("rightSplit",haplotyping.index.Database.getUint(numberOfBases))]),
                   ("direct",[("link", haplotyping.index.Database.getUint(4*numberOfKmers)),
                              ("left",
                               [("distinct","uint8"),
                                ("number",haplotyping.index.Database.getUint(self.maximumNumber))]),
                              ("right",
                               [("distinct","uint8"),
                                ("number",haplotyping.index.Database.getUint(self.maximumNumber))])
                             ]),
                   ("partition",haplotyping.index.Database.getUint(numberOfKmers)),
                   ("cycle",[("number",haplotyping.index.Database.getUint(self.maximumNumber))]),
                   ("reversal",[("number",haplotyping.index.Database.getUint(self.maximumNumber))]),
                   ("paired",[("link",haplotyping.index.Database.getUint(numberOfKmers)),
                              ("number",haplotyping.index.Database.getUint(self.maximumNumber))])]
        dtCkmer=np.dtype(dtypeCkmerList)
        dsCkmer=self.h5file["/split/"].create_dataset("ckmer",(numberOfKmers,), dtype=dtCkmer, chunks=None, 
                                                      compression="gzip", compression_opts=9)
        self._logger.info("store {} splitting k-mers".format(numberOfKmers))
        #add stored and grouped kmers to the final unchunked storage
        for i in range(0,numberOfKmers,Splits.stepSizeStorage):
            stepData = tableSortedKmers[i:i+Splits.stepSizeStorage]
            dsCkmer[i:i+Splits.stepSizeStorage] = stepData
            for row in stepData:
                canonicalSplitKmers+=1
                if row[1].decode()=="l":
                    canonicalSplitKmersLeft+=1
                elif row[1].decode()=="b":
                    canonicalSplitKmersBoth+=1
                elif row[1].decode()=="r":
                    canonicalSplitKmersRight+=1
        dsCkmer.flush()
        # BASE STORAGE - don't make the structure unnecessary big
        dtypeBaseList=[("base","S"+str(self.k-1)),
                   ("number",haplotyping.index.Database.getUint(self.maximumNumber)),
                   ("branches",[])]
        for i in range(len(haplotyping.index.Database.letters)):
            letter = haplotyping.index.Database.letters[i]
            dtypeBaseList[2][1].append((letter,[("number",haplotyping.index.Database.getUint(self.maximumNumber)),
                                                ("ckmerLink",numberOfKmers)]))
        dtBase=np.dtype(dtypeBaseList)
        dsBase=self.h5file["/split/"].create_dataset("base",(numberOfBases,), dtype=dtBase, chunks=None,
                                                     compression="gzip", compression_opts=9)
        self._logger.info("store {} bases".format(numberOfBases))
        #add stored and grouped bases to the final unchunked storage
        for i in range(0,numberOfBases,Splits.stepSizeStorage):
            stepData = tableSortedBases[i:i+Splits.stepSizeStorage]
            dsBase[i:i+Splits.stepSizeStorage] = stepData
        # HISTOGRAM K-MER STORAGE - don't make the structure unnecessary big
        maximumFrequency = max(self.frequencyHistogram["kmer"].keys())
        maximumNumber = max(self.frequencyHistogram["kmer"].values())
        dtypeFrequencyHistogramKmerList=[
                    ("frequency",haplotyping.index.Database.getUint(maximumFrequency)),
                    ("number",haplotyping.index.Database.getUint(maximumNumber)),]
        dtFrequencyHistogramKmer=np.dtype(dtypeFrequencyHistogramKmerList)
        dsFrequencyHistogramKmer=self.h5file["/histogram/"].create_dataset("kmer",(len(self.frequencyHistogram["kmer"]),), 
                                                                  dtype=dtFrequencyHistogramKmer, chunks=None, 
                                                                  compression="gzip", compression_opts=9)
        self._logger.info("store {} entries k-mer histogram".format(len(self.frequencyHistogram["kmer"])))
        #store histogram
        dsFrequencyHistogramKmer[0:len(self.frequencyHistogram["kmer"])] = list(
            sorted(self.frequencyHistogram["kmer"].items()))
        # HISTOGRAM SPLITTING K-MER STORAGE - don't make the structure unnecessary big
        maximumFrequency = max(self.frequencyHistogram["ckmer"].keys())
        maximumNumber = max(self.frequencyHistogram["ckmer"].values())
        dtypeFrequencyHistogramCkmerList=[
                    ("frequency",haplotyping.index.Database.getUint(maximumFrequency)),
                    ("number",haplotyping.index.Database.getUint(maximumNumber)),]
        dtFrequencyHistogramCkmer=np.dtype(dtypeFrequencyHistogramCkmerList)
        dsFrequencyHistogramCkmer=self.h5file["/histogram/"].create_dataset("ckmer",(len(self.frequencyHistogram["ckmer"]),), 
                                                                  dtype=dtFrequencyHistogramCkmer, chunks=None, 
                                                                  compression="gzip", compression_opts=9)
        self._logger.info("store {} entries splitting k-mer histogram".format(len(self.frequencyHistogram["ckmer"])))
        #store histogram
        dsFrequencyHistogramCkmer[0:len(self.frequencyHistogram["ckmer"])] = list(
            sorted(self.frequencyHistogram["ckmer"].items()))
        # HISTOGRAM SPLITTING K-MER BASES STORAGE - don't make the structure unnecessary big
        maximumFrequency = max(self.frequencyHistogram["base"].keys())
        maximumNumber = max(self.frequencyHistogram["base"].values())
        dtypeFrequencyHistogramBaseList=[
                    ("frequency",haplotyping.index.Database.getUint(maximumFrequency)),
                    ("number",haplotyping.index.Database.getUint(maximumNumber)),]
        dtFrequencyHistogramBase=np.dtype(dtypeFrequencyHistogramBaseList)
        dsFrequencyHistogramBase=self.h5file["/histogram/"].create_dataset("base",(len(self.frequencyHistogram["base"]),), 
                                                                  dtype=dtFrequencyHistogramBase, chunks=None, 
                                                                  compression="gzip", compression_opts=9)
        self._logger.info("store {} entries splitting k-mer bases histogram".format(len(self.frequencyHistogram["base"])))
        #store histogram
        dsFrequencyHistogramBase[0:len(self.frequencyHistogram["base"])] = list(
            sorted(self.frequencyHistogram["base"].items()))
        #store the stats
        self._logger.debug("found {} canonicalSplitKmersLeft".format(canonicalSplitKmersLeft))
        self._logger.debug("found {} canonicalSplitKmersRight".format(canonicalSplitKmersRight))
        self._logger.debug("found {} canonicalSplitKmersBoth".format(canonicalSplitKmersBoth))
        self.h5file["/config/"].attrs["maximumCanonicalSplitFrequency"]=self.maximumNumber
        self.h5file["/config/"].attrs["numberCanonicalSplit"]=canonicalSplitKmers
        self.h5file["/config/"].attrs["numberCanonicalSplitLeft"]=canonicalSplitKmersLeft
        self.h5file["/config/"].attrs["numberCanonicalSplitBoth"]=canonicalSplitKmersBoth
        self.h5file["/config/"].attrs["numberCanonicalSplitRight"]=canonicalSplitKmersRight 
        
#-------------------------------------------------
# COMPUTATION AUTOMATON AND INDEX SPLITTING K-MER
#-------------------------------------------------

    def deleteAutomatonWithIndex(filenameBase, k):
        indexFile = "{}_{}.index.splits".format(filenameBase,k)
        automatonFile = "{}_{}.automaton.splits".format(filenameBase,k)
        automatonStatsFile = "{}_{}.automaton.splits.stats".format(filenameBase,k)
        if os.path.exists(indexFile):
            os.remove(indexFile)
        if os.path.exists(automatonFile):
            os.remove(automatonFile)
        if os.path.exists(automatonStatsFile):
            os.remove(automatonStatsFile)
        
    def createAutomatonWithIndex(h5file, filenameBase, k):
        k = min(h5file["/config"].attrs["k"],k)
        logger = logging.getLogger(__name__)
        logger.info("create automaton with k' = {}".format(k))
        indexFile = "{}_{}.index.splits".format(filenameBase,k)
        automatonFile = "{}_{}.automaton.splits".format(filenameBase,k)
        automatonStatsFile = "{}_{}.automaton.splits.stats".format(filenameBase,k)
        if os.path.exists(indexFile) and os.path.exists(automatonFile)and os.path.exists(automatonStatsFile):
            logger.debug("detected previously generated automaton and index")
            with open(automatonStatsFile, "rb") as f:
                stats = pickle.load(f)
            logger.debug("automaton with {} words, size {} MB".format(stats["words_count"],
                                                                          round(stats["real_size"]/1048576)))
            logger.debug("associated index filesize {} MB".format(round(os.stat(indexFile).st_size/1048576)))
            return (max(stats["real_size"],stats["total_size"]),indexFile, automatonFile)
        else:
            logger.debug("generate automaton and index")
            tmpAutomatonFile = "{}.tmp".format(automatonFile)
            tmpIndexFile = "{}.tmp".format(indexFile)
            try:
                process = psutil.Process(os.getpid())
                memoryBefore = process.memory_info().rss
                automatonSplits = ahocorasick.Automaton()
                numberOfKmers = h5file["/split/ckmer"].shape[0]
                kmers = h5file["/split/ckmer"]
                with open(tmpIndexFile, "w") as f:
                    for i in range(0,numberOfKmers,Splits.stepSizeStorage):
                        kmerSubset = kmers[i:i+Splits.stepSizeStorage]
                        for j in range(len(kmerSubset)):
                            id = i + j
                            row = kmerSubset[j]                
                            kmer = row[0].decode()
                            #store in index
                            f.write(kmer)
                            #build automaton for k'-mers with k'<=k describing:
                            #- number of canonical k-mers to check starting from the provided index-location
                            #- index-location of the first matching canonical k-mer
                            rkmer = haplotyping.General.reverse_complement(kmer[-k:])
                            kmer = kmer[:k]
                            if (not automatonSplits.exists(kmer)):
                                automatonSplits.add_word(kmer,(1,id))
                            else:
                                value=list(automatonSplits.get(kmer))
                                if value[0]==0:
                                    value[1]=id
                                    value[0]=1
                                else:
                                    value[0]=id-value[1]+1
                                automatonSplits.add_word(kmer,tuple(value))
                            if not kmer==rkmer:
                                if (not automatonSplits.exists(rkmer)):
                                    automatonSplits.add_word(rkmer,(0,id))
                #create and store automaton
                automatonSplits.make_automaton()
                memoryAfter = process.memory_info().rss
                stats = automatonSplits.get_stats()
                stats["real_size"] = memoryAfter - memoryBefore
                logger.debug("automaton with {} words, size {} MB".format(stats["words_count"],
                                                                          round(stats["real_size"]/1048576)))
                logger.debug("associated index filesize {} MB".format(round(os.stat(tmpIndexFile).st_size/1048576)))
                automatonSplits.save(tmpAutomatonFile, pickle.dumps)
                #clear
                automatonSplits.clear()
                del automatonSplits
                #copy all
                shutil.copyfile(tmpAutomatonFile, automatonFile)
                shutil.copyfile(tmpIndexFile, indexFile)
                #also store stats
                with open(automatonStatsFile, "wb") as f:
                    pickle.dump(stats, f)
            except Exception as ex:
                logger.debug("problem while creating automaton: {}".format(ex))
            finally:
                if os.path.exists(tmpIndexFile):
                    os.remove(tmpIndexFile)
                if os.path.exists(tmpAutomatonFile):
                    os.remove(tmpAutomatonFile)
                    
                
        return (max(stats["real_size"],stats["total_size"]),indexFile, automatonFile)
        
   

    