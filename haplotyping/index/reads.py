import logging, h5py, tables, gzip, time
import os, sys, math, signal
import re, haplotyping
import numpy as np
import haplotyping.index.storage
import multiprocessing as mp
from queue import Empty

class Reads:
    
    """
    Internal use, parse read files and store results in database
    """
    
    def __init__(self, unpairedReadFiles, pairedReadFiles, h5file, filenameBase, debug):
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        if len(unpairedReadFiles)>0:
            self._logger.info("parse "+str(len(unpairedReadFiles))+" unpaired readfile(s)")
        if len(pairedReadFiles)>0:
            self._logger.info("parse "+str(2*len(pairedReadFiles))+" paired readfile(s)")
        
        #self.problemPattern = re.compile(r"["+"".join(haplotyping.index.Database.letters)+
        #                                 "][^"+"".join(haplotyping.index.Database.letters)+"]+")
        
        self.debug = debug
        self.filenameBase = filenameBase
        self.automatonFile = filenameBase+".automaton.splits"
        self.indexFile = filenameBase+".index.splits"
        self.connectedFile = filenameBase+".connected.splits"
        
        #automaton and index should be available
        assert os.path.exists(self.automatonFile)
        assert os.path.exists(self.indexFile)
        
        #set variables
        self.k = h5file["/config"].attrs["k"]
        self.automatonKmerSize = h5file["/config"].attrs["automatonKmerSize"]
        self.maxFrequency = h5file["/config"].attrs["maxFrequency"]
        self.maxProcesses = h5file["/config"].attrs["maxProcesses"]
        self.maxProcessesAutomaton = h5file["/config"].attrs["maxProcessesAutomaton"]
        self.maxProcessesIndex = h5file["/config"].attrs["maxProcessesIndex"]
        self.maxProcessesMatches = h5file["/config"].attrs["maxProcessesMatches"]
        self.maxProcessesConnections = h5file["/config"].attrs["maxProcessesConnections"]
        self.maxProcessesMerges = h5file["/config"].attrs["maxProcessesMerges"]
        self.h5file = h5file
        self.unpairedReadFiles = unpairedReadFiles
        self.pairedReadFiles = pairedReadFiles
        self.numberOfKmers = self.h5file["/split/ckmer"].shape[0]
        
        self.readLengthMinimum=None
        self.readLengthMaximum=None
        self.readPairedTotal=0
        self.readUnpairedTotal=0
        self.readTotal=0
        self.processReadsTime=0
        
        self.maxConnectedIndexSize = 0
        self.maxConnectedIndexLength = 0
        self.maxConnectedIndexNumber = 0
        self.maxConnectedPairNumber = 0        
        
        self.dumpDirect = {}
        
        #theoretical maximum is 2 * #letters, however read-errors will introduce additional connections
        #partly they will be pre-filtered, but a buffer seems sensible (?)
        self.arrayNumberDirect = 2*len(haplotyping.index.Database.letters)
        self.arrayNumberConnection = 10
        
        #get config
        self.minimumFrequency = h5file["/config"].attrs["minimumFrequency"]
        
        #check existence group
        if not "/relations" in h5file:
            h5file.create_group("/relations")
        if not "/connections" in h5file:
            h5file.create_group("/connections")    
        
        #create relations and connections datasets
        if "/relations/direct" in h5file:
            self._logger.warning("direct relation dataset already exists in hdf5 storage")
        elif "/relations/cycle" in h5file:
            self._logger.warning("cycle relation dataset already exists in hdf5 storage")
        elif "/relations/reversal" in h5file:
            self._logger.warning("reversal relation dataset already exists in hdf5 storage")
        elif "/connections/list" in h5file:
            self._logger.warning("connections list dataset already exists in hdf5 storage")
        else:
            #process
            self._processReadFiles()
            self._store()
            #flush
            self.h5file.flush()                                       
                
              
    def _processReadFiles(self):
        
        def close_queue(queue_entry):
            try:
                item = queue_entry.get(block=False)
                while item:
                    try:
                        queue_entry.get(block=False)
                    except:
                        break
            except:
                pass
            queue_entry.close()
            queue_entry.join_thread()
            
        def collect_and_close_queue(queue_entry):
            entries = []
            queue_entry.put(None)
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
            close_queue(queue_entry)
            return entries
            
        mp.set_start_method("fork")
        
        #create shared memory k-mer index (to confirm results from automaton)
        shm_index_size = os.path.getsize(self.indexFile)
        shm_index = mp.shared_memory.SharedMemory(create=True, size=shm_index_size)
        with open(self.indexFile, "r") as f:
            shm_index.buf[0:shm_index_size] = bytes(f.read(),"utf-8")
        self._logger.debug("created shared memory {} bytes k-mer index".format(shm_index_size))
            
        #create shared memory k-mer type, number and bases
        shm_kmer_link = np.dtype(haplotyping.index.Database.getUint(self.numberOfKmers)).type
        shm_kmer_number = np.dtype(haplotyping.index.Database.getUint(self.maxFrequency)).type
        shm_kmer_size = self.numberOfKmers*(1+shm_kmer_number(0).nbytes+(2*shm_kmer_link(0).nbytes))
        shm_kmer = mp.shared_memory.SharedMemory(create=True, size=shm_kmer_size)
        kmer_properties = np.ndarray((self.numberOfKmers,), dtype=[("type","S1"),("number",shm_kmer_number),
                           ("left",shm_kmer_link),("right",shm_kmer_link)], buffer=shm_kmer.buf)
        ckmerLink = 0
        stepSizeStorage=1000000
        for i in range(0,self.numberOfKmers,stepSizeStorage):
            ckmers = self.h5file["/split/ckmer"][i:min(self.numberOfKmers,i+stepSizeStorage)]
            for row in ckmers:
                kmer_properties[ckmerLink] = (row[1],row[2],row[3][0],row[3][1],)
                ckmerLink+=1
        self._logger.debug("created shared memory {} bytes k-mer properties".format(shm_kmer_size))

        shutdown_event = mp.Event()
        
        qsize = 10000
        queue_automaton = mp.Queue(qsize)
        queue_index = mp.Queue(qsize)
        queue_matches = mp.Queue(qsize)
        queue_connections = mp.Queue(qsize)
        queue_finished = mp.Queue()
        queue_storageDirect = mp.Queue()
        queue_storageConnections = mp.Queue()
        
        #auto distribute within limits
        nWorkers = self.maxProcesses - 1
        nWorkersAutomaton = min(self.maxProcessesAutomaton,math.floor(nWorkers/4))
        nWorkersMatches = min(self.maxProcessesMatches,math.floor((nWorkers - nWorkersAutomaton)/3))
        nWorkersIndex = min(self.maxProcessesIndex,math.floor((nWorkers - nWorkersAutomaton - nWorkersMatches)/2))
        nWorkersConnections = min(self.maxProcessesConnections,
                                  (nWorkers - nWorkersAutomaton - nWorkersMatches - nWorkersIndex))
        
        nWorkersLeft = nWorkers - nWorkersAutomaton - nWorkersMatches - nWorkersIndex - nWorkersConnections
        while(nWorkersLeft>0):
            print("change")
            break
        
        self._logger.debug("start {} processes to parse reads with reduced automaton".format(nWorkersAutomaton))
        self._logger.debug("start {} processes to check matches with index".format(nWorkersIndex))
        self._logger.debug("start {} processes to process matches".format(nWorkersMatches))
        self._logger.debug("start {} processes to process connections".format(nWorkersConnections))
        
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool_automaton = mp.Pool(nWorkersAutomaton, haplotyping.index.storage.Storage.worker_automaton, 
                             (shutdown_event,queue_automaton,queue_index,queue_finished,
                              self.k,self.automatonKmerSize,self.automatonFile,))
        pool_index = mp.Pool(nWorkersIndex, haplotyping.index.storage.Storage.worker_index, 
                             (shutdown_event,queue_index,queue_matches,queue_finished,
                              self.k,shm_index.name))
        pool_matches = mp.Pool(nWorkersMatches, haplotyping.index.storage.Storage.worker_matches, 
                               (shutdown_event,queue_matches,queue_connections,queue_storageDirect,queue_finished,
                                self.filenameBase,self.numberOfKmers,self.arrayNumberDirect))
        pool_connections = mp.Pool(nWorkersConnections, haplotyping.index.storage.Storage.worker_connections, 
                               (shutdown_event,queue_connections,queue_storageConnections,queue_finished,
                                self.filenameBase,self.numberOfKmers,
                                self.arrayNumberConnection,self.maxFrequency,shm_kmer.name))
        signal.signal(signal.SIGINT, original_sigint_handler)
        
        try:
            #process and register unpaired read files
            if not "unpairedReads" in self.h5file["/config/"].keys():
                dtypeList = [("file","S255"),("readLength","uint64"),
                     ("readNumber","uint64"),("processTime","uint32")]
            
                ds = self.h5file["/config/"].create_dataset("unpairedReads",(len(self.unpairedReadFiles),),
                                                  dtype=np.dtype(dtypeList),chunks=None)
                for i in range(len(self.unpairedReadFiles)):
                    (readLength,readNumber,processTime) = self._processReadFile(self.unpairedReadFiles[i], 
                                                                   queue_automaton, queue_index, queue_matches)
                    ds[i] = (self.unpairedReadFiles[i],
                             readLength,readNumber,int(processTime))
            else:
                self._logger.error("unpairedReads already (partly) processed")

            #process and register paired read files
            if not "pairedReads" in self.h5file["/config/"].keys():
                dtypeList = [("file0","S255"),("file1","S255"),("readLength","uint64"),
                         ("readNumber","uint64"),("processTime","uint32")]
                ds = self.h5file["/config/"].create_dataset("pairedReads",(len(self.pairedReadFiles),),
                                                  dtype=np.dtype(dtypeList),chunks=None)
                for i in range(len(self.pairedReadFiles)):
                    (readLength,readNumber,processTime) = self._processPairedReadFiles(self.pairedReadFiles[i][0],
                                                                     self.pairedReadFiles[i][1], 
                                                                     queue_automaton, queue_index, queue_matches)
                    ds[i] = (self.pairedReadFiles[i][0],self.pairedReadFiles[i][1],
                             readLength,readNumber,int(processTime))
            else:
                self._logger.error("pairedReads already (partly) processed")        

            #now wait until queues are empty
            while not (queue_automaton.empty() and queue_index.empty() and queue_matches.empty() 
                       and queue_connections.empty()):
                time.sleep(1)
                
            #then trigger stopping by sending enough Nones
            for i in range(pool_automaton._processes):
                queue_automaton.put(None)
            for i in range(pool_index._processes):
                queue_index.put(None)
            for i in range(pool_matches._processes):
                queue_matches.put(None)
            for i in range(pool_connections._processes):
                queue_connections.put(None)
                
            #now wait until everyone is finished
            finishedAutomaton=0
            finishedIndex=0
            finishedMatches=0
            finishedConnections=0
            while not (finishedAutomaton==nWorkersAutomaton and finishedIndex==nWorkersIndex
                       and finishedMatches==nWorkersMatches and finishedConnections==nWorkersConnections):
                item = queue_finished.get(block=True, timeout=1)
                if item=="automaton":
                    finishedAutomaton+=1
                elif item=="index":
                    finishedIndex+=1
                elif item=="matches":
                    finishedMatches+=1
                elif item=="connections":
                    finishedConnections+=1
                else:
                    self._logger.error("unexpected value in finished queue: {}".format(item))
                time.sleep(1)
                
        except KeyboardInterrupt:
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            self._logger.debug("caught KeyboardInterrupt")
            #shutdown directly
            shutdown_event.set()
            #release memory
            shm_index.close()
            shm_index.unlink()
            signal.signal(signal.SIGINT, original_sigint_handler)
        finally:
            #terminate pools
            pool_automaton.terminate()
            pool_index.terminate()
            pool_matches.terminate()
            pool_connections.terminate()
            #close queus workers
            close_queue(queue_automaton)
            close_queue(queue_index)
            close_queue(queue_matches)
            close_queue(queue_connections)
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
            storageDirectFiles = collect_and_close_queue(queue_storageDirect)
            storageConnectionsFiles = collect_and_close_queue(queue_storageConnections)

        #now all files are created, so merging can start
        self._logger.debug("merge {} files with direct connections".format(len(storageDirectFiles)))
        self._logger.debug("merge {} files with additional connections".format(len(storageConnectionsFiles)))
            
        #reset shutdown event
        shutdown_event.clear()
        
        nWorkersMerges = min(self.maxProcessesMerges,self.maxProcesses - 1)
        self._logger.debug("start {} processes to merge stored data".format(nWorkersMerges))
        
        queue_ranges = mp.Queue(nWorkersMerges)
        queue_merges = mp.Queue(nWorkersMerges)
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool_merges = mp.Pool(nWorkersMerges, haplotyping.index.storage.Storage.worker_merges, 
                               (shutdown_event,queue_ranges,queue_merges,storageDirectFiles,
                                self.filenameBase,self.numberOfKmers,self.arrayNumberDirect,
                                self.maxFrequency,self.minimumFrequency,shm_kmer.name))
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
            for item in storageDirectFiles:
                os.remove(item)
            #for item in storageConnectionsFiles:
            #    os.remove(item)
                
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
            close_queue(queue_merges)
            
            self._logger.debug("combine {} merged files".format(len(mergeFiles)))
            
            #combine
            haplotyping.index.storage.Storage.combine_merges(mergeFiles,self.filenameBase,self.numberOfKmers)
            
            #clean
            for item in mergeFiles:
                os.remove(item)
                
        except KeyboardInterrupt:
            original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
            self._logger.debug("caught KeyboardInterrupt")
            shutdown_event.set()
            #release memory
            shm_kmer.close()
            shm_kmer.unlink()
            signal.signal(signal.SIGINT, original_sigint_handler)
        finally:
            #terminate pool
            pool_merges.terminate()
            #close queus workers
            close_queue(queue_ranges)
            close_queue(queue_merges)
            #shutdown
            shutdown_event.set()
            #join pool
            pool_merges.join()
            #release memory
            shm_kmer.close()
            shm_kmer.unlink()
            
        self._logger.debug("created merged file")
            
    def _processReadFile(self, filename: str, queue_automaton, queue_index, queue_matches):
        startTime = time.time()
        with gzip.open(filename, "rt") as f, open(self.connectedFile, "a") as fc:
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
                    queue_automaton.put(sequence)
                    #for matches in matchesList:
                    #    self._processMatches(matches, fc)
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
        self._logger.info("processed {} reads".format(readNumber)) 
        return (readLengthMaximum,readNumber,endTime-startTime)

    def _processPairedReadFiles(self, filename0: str, filename1: str, queue_automaton, queue_index, queue_matches):
        startTime = time.time()
        with gzip.open(filename0, "rt") as f0, gzip.open(filename1, "rt") as f1, open(self.connectedFile, "a") as fc:
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
                    if sequence1[0:31] in sequence0:
                        pos = sequence0.find(sequence1[0:self.k])
                        rpos = sequence0.rfind(sequence1[0:self.k])
                        if pos==rpos:
                            match = sequence0[pos:]
                            if sequence1[0:len(match)]==match:
                                #process as single read because of minimal glue match of size k
                                queue_automaton.put((sequence0[0:pos]+sequence1,))
                            else:
                                queue_automaton.put((sequence0,sequence1,))
                        else:
                            queue_automaton.put((sequence0,sequence1,))
                    else:
                        queue_automaton.put((sequence0,sequence1,))
                    if readNumber%10000==0:
                        self._logger.debug("- processed {} paired reads".format(readNumber)) 
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
            self._logger.info("processed {} paired reads".format(readNumber)) 
            self.processReadsTime+=endTime-startTime
        return (readLengthMaximum,readNumber,endTime-startTime)
        
    def _store(self):
        self.h5file["/config/"].attrs["readLengthMinimum"]=self.readLengthMinimum
        self.h5file["/config/"].attrs["readLengthMaximum"]=self.readLengthMaximum
        self.h5file["/config/"].attrs["readPairedTotal"]=self.readPairedTotal
        self.h5file["/config/"].attrs["readUnpairedTotal"]=self.readUnpairedTotal
        self.h5file["/config/"].attrs["readTotal"]=self.readTotal
        self.h5file["/config/"].attrs["processReadsTime"]=int(np.ceil(self.processReadsTime))        
        
        #store merged direct data
        haplotyping.index.storage.Storage.store_merged_direct(self.h5file, self.filenameBase, self.numberOfKmers, 
                                                       self.minimumFrequency)

    
    