import logging, h5py, tables, gzip, time
import os, sys, math, signal, psutil
import re, haplotyping
import numpy as np
import haplotyping.index.storage
import haplotyping.index.splits
import haplotyping.index.database
import multiprocessing as mp
from threading import Event
from queue import Empty

class Connections:
    
    """
    Internal use, parse read files and store results in database
    """
    
    stepSizeStorage = 1000000
    
    def __init__(self, unpairedReadFiles, pairedReadFiles, h5file, filenameBase, 
                 indexType=None, debug=False, keepTemporaryFiles=False):
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        if len(unpairedReadFiles)>0:
            self._logger.info("parse "+str(len(unpairedReadFiles))+" unpaired readfile(s)")
        if len(pairedReadFiles)>0:
            self._logger.info("parse "+str(2*len(pairedReadFiles))+" paired readfile(s)")                
                    
        self.debug = debug
        self.keepTemporaryFiles = keepTemporaryFiles
        self.indexType = indexType
        self.filenameBase = filenameBase
        
        #set variables
        self.k = h5file["/config"].attrs["k"]
        self.automatonKmerSize = h5file["/config"].attrs["automatonKmerSize"]
        self.maximumFrequency = h5file["/config"].attrs["maximumCanonicalSplitFrequency"]
        self.minimumFrequency = h5file["/config"].attrs["minimumCanonicalSplitFrequency"]
        self.maximumMemory = h5file["/config"].attrs["maximumMemory"]
        self.maximumProcesses = h5file["/config"].attrs["maximumProcesses"]
        self.numberOfKmers = h5file["/split/ckmer"].shape[0]
        self.totalNumberOfKmers = h5file["/config"].attrs["numberKmers"]
        self.h5file = h5file
        self.unpairedReadFiles = unpairedReadFiles
        self.pairedReadFiles = pairedReadFiles
        self.temporaryMergedRelationsFile = None
        
        #statistics
        self.readLengthMinimum=None
        self.readLengthMaximum=None
        self.readPairedTotal=0
        self.readUnpairedTotal=0
        self.readTotal=0
        self.totalReadLength=0
        self.processReadsTime=0
        
        #estimated size
        self.estimatedMaximumReadLength = 0
        
        #check existence group
        if not "/relations" in h5file:
            h5file.create_group("/relations")
        
        #create relations datasets
        if "/relations/direct" in h5file:
            self._logger.warning("direct relation dataset already exists in hdf5 storage")
        elif "/relations/cycle" in h5file:
            self._logger.warning("cycle relation dataset already exists in hdf5 storage")
        elif "/relations/reversal" in h5file:
            self._logger.warning("reversal relation dataset already exists in hdf5 storage")
        else:
            #theoretical maximum is 2 * #letters, however read-errors will introduce additional connections
            #partly they will be pre-filtered, but a buffer seems sensible (?)
            self.numberDirectArray = 2*len(haplotyping.index.Database.letters)
            h5file["/config/"].attrs["numberDirectArray"] = self.numberDirectArray
                                    
            #process
            try:                                                
                #create automaton and index
                automatonKmerSize = (math.ceil((self.k+1)/2) 
                                     if self.automatonKmerSize==0 else min(self.k,self.automatonKmerSize))  
                self.automatonKmerSize = automatonKmerSize
                (automatonMemory,indexFile, automatonFile) = haplotyping.index.splits.Splits.createAutomatonWithIndex(
                    self.h5file, filenameBase, automatonKmerSize)
                #process
                pytablesFile = filenameBase+"_tmp_connections_merge.h5"
                if os.path.exists(pytablesFile):
                    os.remove(pytablesFile)
                self._logger.debug("store temporary in "+pytablesFile)    
                #create datasets
                with tables.open_file(pytablesFile, mode="w", title="Temporary storage") as pytablesStorage:
                    self._processReadFiles(indexFile, automatonFile, automatonMemory, pytablesStorage)
                    self._storeDirect(pytablesStorage)
                    self.h5file.flush()
                    if not self.indexType==haplotyping.index.database.Database.ONLYDIRECTCONNECTIONS:
                        self._processReads(pytablesStorage)
                        self._storeReads(pytablesStorage)
                        self.h5file.flush()     
            except Exception as e:
               self._logger.error("problem occurred while processing reads: "+str(e))
            finally:
                try:
                    if not self.keepTemporaryFiles:
                        os.remove(pytablesFile)
                        haplotyping.index.splits.Splits.deleteAutomatonWithIndex(filenameBase, automatonKmerSize)
                except:
                    self._logger.error("problem removing files")

#--------------
# Handle Queue
#--------------    
                    
    def _collect_queue(queue_entry):
        time.sleep(0.1)
        while True:
            try:
                item = queue_entry.get(block=True, timeout=1)
            except Empty:
                break

    def _close_queue(queue_entry):
        time.sleep(0.1)
        queue_entry.close()
        queue_entry.join_thread()

    def _collect_and_close_queue(queue_entry):
        entries = []
        queue_entry.put(None)
        time.sleep(0.1)
        while True:
            try:
                item = queue_entry.get(block=True, timeout=1)
                if item==None:
                    break
                else:
                    entries.append(item)
            except Empty:
                break 
        #now also close this queue
        Connections._close_queue(queue_entry)
        return entries
    
#---------------
# Handle Memory
#---------------

    def _processMemory():
        process = psutil.Process(os.getpid())
        memory = process.memory_info().rss
        for child in process.children(recursive=True):
            memory += child.memory_info().rss
        return memory
    
#-----------------------------------
# Main functions Direct Connections
#-----------------------------------

    def _processReadFiles(self, indexFile, automatonFile, automatonMemory, pytablesStorage):
          
        def estimateIndexMemory(nWorkersAutomaton,workerAutomatonMemory,nWorkersMatches,
                           workerMatchesMemory,nWorkersIndex,workerIndexMemory):
            return ((nWorkersAutomaton*workerAutomatonMemory) + 
                    (nWorkersMatches*workerMatchesMemory) + 
                    (nWorkersIndex*workerIndexMemory))
                
        #get method
        self._logger.debug("using method '{}' for multiprocessing".format(mp.get_start_method()))
        process = psutil.Process(os.getpid())
        self._logger.debug("initially used memory {} MB".format(math.ceil(Connections._processMemory()/1048576)))
        
        self._logger.debug("memory info: {}".format(psutil.Process(os.getpid()).memory_info()))
        
        #compute size shared memory k-mer index (to confirm results from automaton)
        shm_index_size = os.path.getsize(indexFile)
        self._logger.debug("size shared memory {} MB k-mer index".format(math.ceil(shm_index_size/1048576)))
                    
        #compute size shared memory k-mer type, number and bases
        shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(self.numberOfKmers)).type
        shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(self.maximumFrequency)).type
        shm_kmer_size = self.numberOfKmers*(1+shm_kmer_number(0).nbytes+(2*shm_kmer_link(0).nbytes))
        self._logger.debug("size shared memory {} MB k-mer properties".format(math.ceil(shm_kmer_size/1048576)))
        
        #shutdown
        shutdown_event = mp.Event()
        
        qsize = 10000
        queue_start = mp.Queue(qsize)
        queue_automaton = mp.Queue(qsize)
        queue_index = mp.Queue(qsize)
        queue_matches = mp.Queue(qsize)
        queue_finished = mp.Queue()
        queue_storageDirect = mp.Queue()
        queue_storageReads = mp.Queue()
        
        #estimate worker automaton memory
        workerAutomatonMemory = automatonMemory
        #estimate worker matches memory
        workerMatchesDtypeEntry = haplotyping.index.storage.Storage.worker_matches_dtype(
            self.numberOfKmers,self.maximumFrequency,self.estimatedMaximumReadLength,self.numberDirectArray)
        workerMatchesMemory = self.numberOfKmers * np.dtype(workerMatchesDtypeEntry).itemsize
        #estimate worker index memory (shared)
        workerIndexMemory = 0
        workerSharedMemory = (shm_kmer_size+shm_index_size)
        
        self._logger.debug("estimated shared memory: {} MB".format(math.ceil(workerSharedMemory/1048576)))
        self._logger.debug("estimated memory automaton worker: {} MB".format(math.ceil(workerAutomatonMemory/1048576)))
        self._logger.debug("estimated memory index worker: {} MB".format(math.ceil(workerIndexMemory/1048576)))
        self._logger.debug("estimated memory matches worker: {} MB".format(math.ceil(workerMatchesMemory/1048576)))
        
        #worker requirements
        nWorkers = max(3,mp.cpu_count()-1) if self.maximumProcesses==0 else self.maximumProcesses - 1
        
        #memory requirements
        if self.maximumMemory>0:
            maximumMemory = min(psutil.virtual_memory().available + process.memory_info().rss, self.maximumMemory)
        else:
            maximumMemory = round(0.95*psutil.virtual_memory().available) + process.memory_info().rss    
            
            
        #compute maximum number of automaton workers (high memory usage)
        #assume that ideal ratio workers is 1:2:4 (to be verified/computed?)
        nWorkersAutomaton = max(1,min(math.floor(nWorkers/3),
                                      math.floor((maximumMemory-workerSharedMemory)/
                                             estimateIndexMemory(1,workerAutomatonMemory,
                                                                 2,workerMatchesMemory,
                                                                 4,workerIndexMemory))))
        #auto distribute other workers within limits
        nWorkersMatches = math.floor((nWorkers - nWorkersAutomaton)/2)
        nWorkersIndex = nWorkers - nWorkersAutomaton - nWorkersMatches        
        #increment if processes available
        nWorkersLeft = nWorkers - nWorkersAutomaton - nWorkersMatches - nWorkersIndex
        while(nWorkersLeft>0):
            if(nWorkersLeft>0):
                nWorkersMatches+=1
                nWorkersLeft-=1
            if(nWorkersLeft>0):
                nWorkersIndex+=1
                nWorkersLeft-=1                                         
        #reduce to fit memory requirements   
        while (estimatedMemory := estimateIndexMemory(nWorkersAutomaton,workerAutomatonMemory,
                             nWorkersMatches,workerMatchesMemory,
                             nWorkersIndex,workerIndexMemory)) + workerSharedMemory > maximumMemory:
            if (nWorkersAutomaton==1) and (nWorkersMatches==1) and (nWorkersIndex==1):
                raise Exception("not enough memory available, required: {} MB".format(round(estimatedMemory/1048576)))
            elif (nWorkersAutomaton==1) and (nWorkersMatches==1):
                nWorkersIndex -= 1
            elif (nWorkersMatches>nWorkersAutomaton):
                nWorkersMatches -= 1
                nWorkersIndex = nWorkersIndex+1
            else:
                nWorkersAutomaton-=1
                nWorkersIndex = nWorkersIndex+1
        #don't oversize the index workers, maximum two times matches workers
        nWorkersIndex = min(nWorkersIndex,2*nWorkersMatches)
        #final calculation memory estimation
        estimatedMemory = estimateIndexMemory(nWorkersAutomaton,workerAutomatonMemory,
                             nWorkersMatches,workerMatchesMemory,nWorkersIndex,workerIndexMemory) + workerSharedMemory
        
        self._logger.debug("start {} processes to parse reads with reduced automaton".format(nWorkersAutomaton))
        self._logger.debug("start {} processes to check index".format(nWorkersIndex))
        self._logger.debug("start {} processes to process matches".format(nWorkersMatches))
        self._logger.debug("estimated total memory usage: {} MB".format(math.ceil(estimatedMemory/1048576)))
                
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool_automaton = mp.Pool(nWorkersAutomaton, haplotyping.index.storage.Storage.workerAutomaton, 
                             (shutdown_event,queue_start,queue_automaton,queue_index,queue_finished,
                              self.k,self.automatonKmerSize,automatonFile,))
        #first start automatons, one at a time, because of memory peak
        startedAutomaton = 0
        queue_start.put("automaton")
        while not (startedAutomaton==nWorkersAutomaton):
            try:
                item = queue_finished.get(block=True, timeout=1)
                if item=="automaton:started":
                    startedAutomaton+=1
                    self._logger.debug("{} of {} automatons started".format(startedAutomaton,nWorkersAutomaton))
                    if startedAutomaton<nWorkersAutomaton:
                        queue_start.put("automaton")
                else:
                    self._logger.error("unexpected start value in finished queue: {}".format(item))
            except:
                pass
            time.sleep(1)
            
        #create shared memory k-mer index (to confirm results from automaton)
        shm_index = mp.shared_memory.SharedMemory(create=True, size=shm_index_size)
        with open(indexFile, "r") as f:
            shm_index.buf[0:shm_index_size] = bytes(f.read(),"utf-8")
        self._logger.debug("created shared memory {} MB k-mer index".format(math.ceil(shm_index_size/1048576)))
        
        #create shared memory k-mer type, number and bases
        shm_kmer = mp.shared_memory.SharedMemory(create=True, size=shm_kmer_size)
        kmer_properties = np.ndarray((self.numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                           ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm_kmer.buf)
        ckmerLink = 0
        for i in range(0,self.numberOfKmers,Connections.stepSizeStorage):
            ckmers = self.h5file["/split/ckmer"][i:min(self.numberOfKmers,i+Connections.stepSizeStorage)]
            stepData = [(row[1],row[2],row[3][0],row[3][1],) for row in ckmers]
            kmer_properties[ckmerLink:ckmerLink+len(stepData)] = stepData
            ckmerLink+=len(stepData)
            del stepData
        self._logger.debug("created shared memory {} MB k-mer properties".format(math.ceil(shm_kmer_size/1048576)))

        #now start other workers
        pool_index = mp.Pool(nWorkersIndex, haplotyping.index.storage.Storage.workerIndex, 
                             (shutdown_event,queue_index,queue_matches,queue_storageReads,queue_finished,
                              self.filenameBase,self.numberOfKmers,self.k,
                              self.indexType,shm_index.name))
        pool_matches = mp.Pool(nWorkersMatches, haplotyping.index.storage.Storage.workerMatches, 
                               (shutdown_event,queue_matches,queue_storageDirect,queue_finished,
                                self.filenameBase,self.numberOfKmers,self.maximumFrequency,
                                self.estimatedMaximumReadLength,self.numberDirectArray,
                                self.indexType,shm_kmer.name))
        signal.signal(signal.SIGINT, original_sigint_handler)

        try:
            #process and register unpaired read files
            if not "unpairedReads" in self.h5file["/config/"].keys():
                dtypeList = [("file","S255"),("readLength","uint64"),
                     ("readNumber","uint64"),("totalReadLength","uint64"),("processTime","uint32")]
            
                ds = self.h5file["/config/"].create_dataset("unpairedReads",(len(self.unpairedReadFiles),),
                                                  dtype=np.dtype(dtypeList),chunks=None, 
                                                  compression="gzip", compression_opts=9)
                for i in range(len(self.unpairedReadFiles)):
                    self._logger.debug("process {}".format(os.path.basename(self.unpairedReadFiles[i])))
                    (readLength,readNumber,totalReadLength,processTime) = self._processReadFile(
                                                                   self.unpairedReadFiles[i], 
                                                                   queue_automaton, queue_index, queue_matches)
                    ds[i] = (self.unpairedReadFiles[i],
                             readLength,readNumber,totalReadLength,int(processTime))
            else:
                self._logger.error("unpairedReads already (partly) processed")

            #process and register paired read files
            if not "pairedReads" in self.h5file["/config/"].keys():
                dtypeList = [("file0","S255"),("file1","S255"),("readLength","uint64"),
                         ("readNumber","uint64"),("totalReadLength","uint64"),("processTime","uint32")]
                ds = self.h5file["/config/"].create_dataset("pairedReads",(len(self.pairedReadFiles),),
                                                  dtype=np.dtype(dtypeList),chunks=None, 
                                                  compression="gzip", compression_opts=9)
                for i in range(len(self.pairedReadFiles)):
                    self._logger.debug("process {} and {}".format(os.path.basename(self.pairedReadFiles[i][0]),
                                                           os.path.basename(self.pairedReadFiles[i][1])))
                    (readLength,readNumber,totalReadLength,processTime) = self._processPairedReadFiles(
                                                                     self.pairedReadFiles[i][0],
                                                                     self.pairedReadFiles[i][1], 
                                                                     queue_automaton, queue_index, queue_matches)
                    ds[i] = (self.pairedReadFiles[i][0],self.pairedReadFiles[i][1],
                             readLength,readNumber,totalReadLength,int(processTime))
            else:
                self._logger.error("pairedReads already (partly) processed")                    
            
            #now wait until queues are empty
            while not (queue_automaton.empty() and queue_index.empty() and queue_matches.empty()):
                time.sleep(1)
                
            #then trigger stopping by sending enough Nones
            for i in range(pool_automaton._processes):
                queue_automaton.put(None)
            for i in range(pool_index._processes):
                queue_index.put(None)
            for i in range(pool_matches._processes):
                queue_matches.put(None)
                
            #now wait until everyone is finished
            finishedAutomaton=0
            finishedIndex=0
            finishedMatches=0
            totalCanonicalSplitFrequencies=0
            while not (finishedAutomaton==nWorkersAutomaton and finishedIndex==nWorkersIndex
                       and finishedMatches==nWorkersMatches):
                try:
                    item = queue_finished.get(block=True, timeout=1)
                    if item.startswith("automaton:ended"):
                        finishedAutomaton+=1
                    elif item.startswith("index:ended"):
                        finishedIndex+=1
                        totalCanonicalSplitFrequencies+=int(item.split(":")[3])
                    elif item.startswith("matches:ended"):
                        finishedMatches+=1
                    else:
                        self._logger.error("unexpected value in finished queue: {}".format(item))
                except:
                    pass
                time.sleep(1)
            #store total found canonical k-mers (doesn't fully match results from kmc because of overlapping sequences)
            self.h5file["/config/"].attrs["totalCanonicalSplitFrequencies"] = totalCanonicalSplitFrequencies
                
        except KeyboardInterrupt:
            self._logger.debug("caught KeyboardInterrupt")
            #collect from queue
            Connections._collect_queue(queue_automaton)
            Connections._collect_queue(queue_index)
            Connections._collect_queue(queue_matches)
            #shutdown directly
            shutdown_event.set()
            #release memory
            shm_index.close()
            shm_kmer.close()
            shm_index.unlink()
            shm_kmer.unlink()
            sys.exit()
        finally:
            #terminate pools
            pool_automaton.terminate()
            pool_index.terminate()
            pool_matches.terminate()
            #close queus workers
            Connections._close_queue(queue_automaton)
            Connections._close_queue(queue_index)
            Connections._close_queue(queue_matches)
            #shutdown
            shutdown_event.set()
            #join pools
            pool_automaton.join()
            pool_index.join()
            pool_matches.join()
            #release memory
            shm_index.close()
            shm_index.unlink()
            #collect created files 
            storageDirectFiles = Connections._collect_and_close_queue(queue_storageDirect)
            self.storageReadFiles = Connections._collect_and_close_queue(queue_storageReads)            
            
        #now all files are created, so merging can start
        self._logger.debug("merge {} files with direct connections".format(len(storageDirectFiles)))
            
        #reset shutdown event
        shutdown_event.clear()
        
        #assume memory is no problem
        nWorkersMerges = max(1,mp.cpu_count()-1) if self.maximumProcesses==0 else self.maximumProcesses - 1
        self._logger.debug("start {} processes to merge stored data".format(nWorkersMerges))
        
        queue_ranges = mp.Queue(nWorkersMerges)
        queue_merges = mp.Queue(nWorkersMerges)
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool_merges = mp.Pool(nWorkersMerges, haplotyping.index.storage.Storage.workerMergeDirect, 
                               (shutdown_event,queue_ranges,queue_merges,
                                storageDirectFiles,
                                self.filenameBase,self.numberOfKmers,
                                self.maximumFrequency,self.minimumFrequency,
                                shm_kmer.name))
        signal.signal(signal.SIGINT, original_sigint_handler)
        
        try:
            #fill queue
            mergeNumber = math.ceil(self.numberOfKmers/nWorkersMerges)
            for mergeStart in range(0,self.numberOfKmers,mergeNumber):
                queue_ranges.put((mergeStart,mergeNumber,))

            #then trigger stopping by sending enough Nones
            for i in range(pool_merges._processes):
                queue_ranges.put(None)
                
            #now wait    
            while not (queue_ranges.empty()):
                time.sleep(1)
                
            #clean
            if not self.keepTemporaryFiles:
                for item in storageDirectFiles:
                    os.remove(item)
                                
            #collect created merges
            mergeFiles = []
            while True:
                try:
                    item = queue_merges.get(block=True, timeout=1)
                    if item==None:
                        break
                    elif isinstance(item,str):
                        mergeFiles.append(item)
                except Empty:
                    break 
            #now also close this queue
            Connections._close_queue(queue_merges)
            
            self._logger.debug("combine {} merged files".format(len(mergeFiles)))
            
            #combine
            haplotyping.index.storage.Storage.combineDirectMerges(
                mergeFiles, pytablesStorage, self.numberOfKmers, self.maximumFrequency)                        
            
            #clean
            if not self.keepTemporaryFiles:
                for item in mergeFiles:
                    os.remove(item)
                
        except KeyboardInterrupt:
            self._logger.debug("caught KeyboardInterrupt")
            shutdown_event.set()
            #release memory
            shm_kmer.close()
            shm_kmer.unlink()
            sys.exit()
        finally:
            #terminate pool
            pool_merges.terminate()
            #close queus workers
            Connections._close_queue(queue_ranges)
            Connections._close_queue(queue_merges)
            #shutdown
            shutdown_event.set()
            #join pool
            pool_merges.join()
            #release memory
            shm_kmer.close()
            shm_kmer.unlink()
            
        self._logger.debug("created merged file")
        
        self._logger.debug("process {} files with read information".format(len(self.storageReadFiles)))        
            
    def _processReadFile(self, filename: str, queue_automaton, queue_index, queue_matches):
        startTime = time.time()
        open_fn = gzip.open if filename.endswith(".gz") else open
        with open_fn(filename, "rt") as f:
            readLengthMinimum=None
            readLengthMaximum=None
            readNumber=0
            totalReadLength=0
            while True:               
                identifier = f.readline().rstrip()
                sequence = f.readline().rstrip()  
                plusline = f.readline().rstrip()
                quality = f.readline().rstrip()
                if sequence:
                    if not (identifier.startswith("@") and plusline.startswith("+") and len(sequence)==len(quality)):
                        self._logger.error("invalid fastq-file {}, line {}".format(filename,4*readNumber))
                        break
                    else:
                        readLengthMinimum=(len(sequence) if readLengthMinimum==None 
                                           else min(readLengthMinimum,len(sequence)))
                        readLengthMaximum=(len(sequence) if readLengthMaximum==None 
                                           else max(readLengthMaximum,len(sequence)))
                        readNumber+=1
                        totalReadLength+=len(sequence)
                        queue_automaton.put(sequence)
                        if readNumber%1000000==0:
                            self._logger.debug("- processed {} reads".format(readNumber)) 
                else:
                    break
            endTime = time.time()
            if readNumber>0:
                if readLengthMinimum>0:
                    self.readLengthMinimum=(readLengthMinimum if self.readLengthMinimum==None 
                                            else min(self.readLengthMinimum,readLengthMinimum))
                self.readLengthMaximum=(readLengthMaximum if self.readLengthMaximum==None 
                                        else max(self.readLengthMaximum,readLengthMaximum))
                self.readUnpairedTotal+=readNumber
                self.readTotal+=readNumber
                self.totalReadLength+=totalReadLength
            self.processReadsTime+=endTime-startTime
        self._logger.info("processed {} reads".format(readNumber)) 
        return (readLengthMaximum,readNumber,totalReadLength,endTime-startTime)

    def _processPairedReadFiles(self, filename0: str, filename1: str, queue_automaton, queue_index, queue_matches):
        startTime = time.time()
        open_fn0 = gzip.open if filename0.endswith(".gz") else open
        open_fn1 = gzip.open if filename1.endswith(".gz") else open
        with open_fn0(filename0, "rt") as f0, open_fn1(filename1, "rt") as f1:
            readLengthMinimum=None
            readLengthMaximum=None
            readNumber=0
            totalReadLength=0
            while True:
                #first of pair
                identifier0 = f0.readline().rstrip() 
                sequence0 = f0.readline().rstrip()  
                plusline0 = f0.readline().rstrip() 
                quality0 = f0.readline().rstrip() 
                #second of pair
                identifier1 = f1.readline().rstrip() 
                sequence1 = haplotyping.General.reverse_complement(f1.readline().rstrip())
                plusline1 = f1.readline().rstrip() 
                quality1 = f1.readline().rstrip() 
                #process
                if sequence0 and sequence1:
                    if not (identifier0.startswith("@") and plusline0.startswith("+") and len(sequence0)==len(quality0)):
                        self._logger.error("invalid fastq-file {}, line {}".format(filename0,2*readNumber))
                        break
                    elif not (identifier1.startswith("@") and plusline1.startswith("+") and len(sequence1)==len(quality1)):
                        self._logger.error("invalid fastq-file {}, line {}".format(filename1,2*readNumber))
                        break
                    else:
                        readLengthMinimum=(min(len(sequence0),len(sequence1)) if readLengthMinimum==None 
                                               else min(readLengthMinimum,len(sequence0),len(sequence1)))
                        readLengthMaximum=(max(len(sequence0),len(sequence1)) if readLengthMaximum==None 
                                           else max(readLengthMaximum,len(sequence0),len(sequence1)))
                        readNumber+=2  
                        totalReadLength+=(len(sequence0)+len(sequence1))
                        if sequence1[0:self.k] in sequence0:
                            pos = sequence0.find(sequence1[0:self.k])
                            rpos = sequence0.rfind(sequence1[0:self.k])
                            if pos==rpos:
                                match = sequence0[pos:]
                                if sequence1[0:len(match)]==match:
                                    #process as single read because of minimal glue match of size k
                                    queue_automaton.put(sequence0[0:pos]+sequence1) 
                                else:
                                    queue_automaton.put((sequence0,sequence1,))                                
                            else:
                                queue_automaton.put((sequence0,sequence1,))                            
                        else:
                            queue_automaton.put((sequence0,sequence1,))                        
                        if readNumber%1000000==0:
                            self._logger.debug("- processed {} paired reads".format(readNumber))                     
                else:
                    break
            endTime = time.time()
            if readNumber>0:
                if readLengthMinimum>0:
                    self.readLengthMinimum=(readLengthMinimum if self.readLengthMinimum==None 
                                            else min(self.readLengthMinimum,readLengthMinimum))
                self.readLengthMaximum=(readLengthMaximum if self.readLengthMaximum==None 
                                        else max(self.readLengthMaximum,readLengthMaximum))
                self.readPairedTotal+=readNumber
                self.readTotal+=readNumber
                self.totalReadLength+=totalReadLength
            self._logger.info("processed {} paired reads".format(readNumber)) 
            self.processReadsTime+=endTime-startTime
        return (readLengthMaximum,readNumber,totalReadLength,endTime-startTime)
        
        
    def _storeDirect(self, pytablesStorage):
        self.h5file["/config/"].attrs["minimumReadLength"]=self.readLengthMinimum
        self.h5file["/config/"].attrs["maximumReadLength"]=self.readLengthMaximum
        self.h5file["/config/"].attrs["numberReadsPaired"]=self.readPairedTotal
        self.h5file["/config/"].attrs["numberReadsUnpaired"]=self.readUnpairedTotal
        self.h5file["/config/"].attrs["numberReads"]=self.readTotal
        self.h5file["/config/"].attrs["totalReadLength"]=self.totalReadLength
        self.h5file["/config/"].attrs["timeProcessReads"]=int(np.ceil(self.processReadsTime))        
        
        #store merged data
        haplotyping.index.storage.Storage.storeMergedDirect(
            self.h5file, pytablesStorage, 
            self.numberOfKmers, self.minimumFrequency)

#-------------------------------------
# Main functions Indirect Connections
#-------------------------------------
        
    def _processReads(self, pytablesStorage):   
        
        if len(self.storageReadFiles)>0:
            
            #shutdown
            shutdown_event = mp.Event()
        
            #partition k-mers based on direct connections
            maxNumberOfPartitions = int(self.numberOfKmers ** (2/3))
            self.numberOfPartitions = haplotyping.index.storage.Storage.partitionKmers(self.h5file, pytablesStorage,
                                                              maxNumberOfPartitions)
            #prepare shared memory with k-mer properties and direct connections
            numberOfDirect = self.h5file["/relations/direct"].shape[0]
            shm_kmer_partition = np.dtype(haplotyping.index.Database.getUint(self.numberOfPartitions)).type
            shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(self.maximumFrequency)).type
            shm_kmer_reference = np.dtype(haplotyping.index.Database.getUint(numberOfDirect)).type
            shm_kmer_numberLeft = np.dtype("uint8").type
            shm_kmer_numberRight = np.dtype("uint8").type
            shm_kmer_size = self.numberOfKmers*(shm_kmer_partition(0).nbytes+
                                                shm_kmer_number(0).nbytes+
                                                shm_kmer_reference(0).nbytes+
                                                shm_kmer_numberLeft(0).nbytes+
                                                shm_kmer_numberRight(0).nbytes)
            shm_direct_kmer = np.dtype(haplotyping.index.Database.getUint(self.numberOfKmers)).type
            shm_direct_size = numberOfDirect*shm_direct_kmer(0).nbytes
            self._logger.debug("size shared memory {} MB k-mer properties".format(math.ceil(shm_kmer_size/1048576)))
            self._logger.debug("size shared memory {} MB direct connections".format(math.ceil(shm_direct_size/1048576)))

            shm_kmer = mp.shared_memory.SharedMemory(create=True, size=shm_kmer_size)
            shm_direct = mp.shared_memory.SharedMemory(create=True, size=shm_direct_size)
            kmer_properties = np.ndarray((self.numberOfKmers,), 
                                         dtype=[("partition",shm_kmer_partition),
                                                ("number",shm_kmer_number),
                                                ("reference",shm_kmer_reference),
                                                ("numberLeft",shm_kmer_numberLeft),
                                                ("numberRight",shm_kmer_numberRight)], 
                                         buffer=shm_kmer.buf)
            direct_properties = np.ndarray((numberOfDirect,), 
                                         dtype=shm_direct_kmer, 
                                         buffer=shm_direct.buf)
            ckmerLink = 0
            for i in range(0,self.numberOfKmers,Connections.stepSizeStorage):
                #get k-mer data
                ckmers = self.h5file["/split/ckmer"][i:min(self.numberOfKmers,i+Connections.stepSizeStorage)]
                stepData = [(row[5],row[2],row[4][0],row[4][1][0],row[4][2][0]) for row in ckmers]
                kmer_properties[ckmerLink:ckmerLink+len(stepData)] = stepData
                ckmerLink+=len(stepData)
                del stepData
                #get direct data
                directStart = numberOfDirect
                directEnd = 0
                for row in ckmers:
                    if row[4][1][1]==0 and row[4][2][1]==0:
                        continue
                    else:
                        directStart=min(directStart,row[4][0])
                        directEnd=max(directEnd,row[4][0]+row[4][1][1]+row[4][2][1])
                if directStart<directEnd:
                    direct = self.h5file["/relations/direct"][directStart:directEnd]
                    stepData = [row[1][0] for row in direct]
                    direct_properties[directStart:directEnd] = stepData
                    del stepData
                                
            self._logger.debug("created shared memory {} MB k-mer properties".format(math.ceil(shm_kmer_size/1048576)))
            self._logger.debug("created shared memory {} MB direct connections".format(math.ceil(shm_direct_size/1048576)))
            
            #worker requirements
            nWorkers = max(3,mp.cpu_count()-1) if self.maximumProcesses==0 else self.maximumProcesses - 1
            nWorkersReads = min(len(self.storageReadFiles),nWorkers)            
            self._logger.debug("start {} processes to process reads".format(nWorkersReads))
            
            #queues
            queue_rawReads = mp.Queue(len(self.storageReadFiles))
            queue_filteredReads = mp.Queue()
            queue_finished = mp.Queue()            
            
            #maximum number of splitting k-mers in read
            maximumReadLength = self.h5file["/config/"].attrs["maximumReadLength"] - self.k + 1
            
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            pool_reads = mp.Pool(nWorkersReads, haplotyping.index.storage.Storage.workerProcessReads, 
                                 (shutdown_event,queue_rawReads,queue_filteredReads,queue_finished,
                                  self.filenameBase,self.numberOfKmers,
                                  self.numberOfPartitions,numberOfDirect,self.maximumFrequency,maximumReadLength,
                                  shm_kmer.name,shm_direct.name))
            signal.signal(signal.SIGINT, original_sigint_handler)
            
            try:
                #fill queue
                for item in self.storageReadFiles:
                    queue_rawReads.put(item)
                #trigger stopping by sending enough Nones
                for i in range(pool_reads._processes):
                    queue_rawReads.put(None)

                #now wait    
                while not (queue_rawReads.empty()):
                    time.sleep(1)

                #clean
                if not self.keepTemporaryFiles:
                   for item in self.storageReadFiles:
                       os.remove(item)
                del self.storageReadFiles

            except KeyboardInterrupt:
                self._logger.debug("caught KeyboardInterrupt")
                shutdown_event.set()
                #release memory
                shm_kmer.close()
                shm_direct.close()
                shm_kmer.unlink()
                shm_direct.unlink()
                sys.exit()
            finally:
                #terminate pool
                pool_reads.terminate()
                #close queus workers
                Connections._close_queue(queue_rawReads)
                Connections._close_queue(queue_finished)
                #shutdown
                shutdown_event.set()
                #join pool
                pool_reads.join()
                #release memory
                shm_kmer.close()
                shm_direct.close()
                shm_kmer.unlink()
                shm_direct.unlink()
                #get filtered readfiles
                storageFilteredReadFiles = Connections._collect_and_close_queue(queue_filteredReads)   
                
            #compute merged paired from storageFilteredReadFiles
            haplotyping.index.storage.Storage.combineFilteredPairs(
                storageFilteredReadFiles, pytablesStorage, self.numberOfKmers)

            #get sizes for partitions
            partitionSizes = [[0,0] for i in range(self.numberOfPartitions)]
            for item in storageFilteredReadFiles:
                with tables.open_file(item, mode="r") as pytablesStorageFiltered:
                    nReads = pytablesStorageFiltered.root.readPartitionInfo.shape[0]
                    for i in range(0,nReads,Connections.stepSizeStorage):
                        stepData = pytablesStorageFiltered.root.readPartitionInfo[i:i+Connections.stepSizeStorage]
                        for row in stepData:
                            partitionSizes[row[1]][0]+=1
                            partitionSizes[row[1]][1]+=row[0]   
            
            #shutdown
            shutdown_event = mp.Event()
        
            #worker requirements
            nWorkers = max(3,mp.cpu_count()-1) if self.maximumProcesses==0 else self.maximumProcesses - 1
            nWorkersMerges = min(len(storageFilteredReadFiles),nWorkers)     
            self._logger.debug("start {} processes to merge reads".format(nWorkersMerges))
            
            queue_ranges = mp.Queue(nWorkersMerges)
            queue_merges = mp.Queue(nWorkersMerges)
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            pool_merges = mp.Pool(nWorkersMerges, haplotyping.index.storage.Storage.workerMergeReads, 
                                   (shutdown_event,queue_ranges,queue_merges,
                                    storageFilteredReadFiles,partitionSizes,self.filenameBase,
                                    self.numberOfKmers,self.numberOfPartitions,maximumReadLength))
            signal.signal(signal.SIGINT, original_sigint_handler)
        
            try:
                #fill queue
                mergeNumber = math.ceil(self.numberOfPartitions/nWorkersMerges)
                for mergeStart in range(0,self.numberOfPartitions,mergeNumber):
                    queue_ranges.put((mergeStart,mergeNumber,))
                    
                #then trigger stopping by sending enough Nones
                for i in range(pool_merges._processes):
                    queue_ranges.put(None)

                #now wait    
                while not (queue_ranges.empty()):
                    time.sleep(1)

                #collect created merges
                mergeFiles = []
                while True:
                    try:
                        item = queue_merges.get(block=True, timeout=1)
                        if item==None:
                            break
                        elif isinstance(item,str):
                            mergeFiles.append(item)
                    except Empty:
                        if len(mergeFiles)==nWorkersMerges:
                            break
                        else:
                            time.sleep(1)
                #now also close this queue
                Connections._close_queue(queue_merges)

                self._logger.debug("combine {} merged files".format(len(mergeFiles)))

                #combine
                haplotyping.index.storage.Storage.combineReadMerges(
                    sorted(mergeFiles), pytablesStorage, 
                    self.numberOfKmers,self.numberOfPartitions,maximumReadLength)

                #clean mergeFiles
                if not self.keepTemporaryFiles:
                    for item in mergeFiles:
                        os.remove(item)
                        
            except KeyboardInterrupt:
                self._logger.debug("caught KeyboardInterrupt")
                shutdown_event.set()
                sys.exit()
            finally:
                #terminate pool
                pool_merges.terminate()
                #close queus workers
                Connections._close_queue(queue_ranges)
                #shutdown
                shutdown_event.set()
                #join pool
                pool_merges.join()
                
            #clean storageFilteredReadFiles
            if not self.keepTemporaryFiles:
                for item in storageFilteredReadFiles:
                    os.remove(item)
                    
    def _storeReads(self, pytablesStorage):
        #store merged data
        haplotyping.index.storage.Storage.storeMergedReads(
            self.h5file, pytablesStorage, 
            self.numberOfKmers,self.numberOfPartitions)
        
    
                
            
            
            
            
