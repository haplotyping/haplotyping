import logging, h5py, os, sys, glob, re, math, shutil
import numpy as np, tables, gzip, csv
import multiprocessing as mp

import haplotyping
import haplotyping.index.splits
import haplotyping.index.reads

class Database:
    
    """
    Constructing splitting k-mer database
    
    Parameters
    ----------------------
    k : int
        The used k-mer size, should correspond with the k-mer size of the sorted k-mer list

    name : str
        Name of the variety

    filenameBase : str
        Base name and location for the splitting k-mer database and the temporary databases

    sortedIndexFile: str
        Location of the sorted k-mer list
        The k-mer size for k-mers in this list should correspond with the provided k parameter

    readFiles: optional, default is empty list
        A list with locations from readFiles
        Optional, but at least readFiles or pairedReadFiles should be defined and non-empty

    pairedReadFiles: optional, default is empty list
        A list with pairs of locations from paired readFiles
        Optional, but at least readFiles or pairedReadFiles should be defined and non-empty

    minimumFrequency: int, optional, default is 2
        The minimum frequency for a splitting k-mer.
        Adjusting this number will change the amount of read error based results but 
        can also introduce gaps; With higher read depth or lower error rate, this number can possibly be increased
        
    maximumProcesses: int, optional, default is 5
        The maximum number of processes including main process, should be at least 5
        
    maximumProcessesAutomaton: int, optional, default is 1
        The maximum number of processes used for direct parsing of reads using the automaton
        Because shared memory can't be used, multiple copies of the automaton have to be in memory
        Be carefull if available memory is limited
    
    maximumProcessesIndex: int, optional
        The maximum number of processes used for creating a (partial) index from the matches
        Making use of shared memory to load index
        
    maximumProcessesMatches: int, optional
        The maximum number of processes used for creating a (partial) index from the matches
        Collecting and organising data is done mostly in memory
        Be carefull if available memory is limited
    
    maximumProcessesConnections: int, optional
        The maximum number of processes used for creating a (partial) index from the matches
        Collecting and organising data is done mostly in memory
        Be carefull if available memory is limited
    
    maximumProcessesMerges: int, optional
        The maximum number of processes used for merging operations
    
    automatonKmerSize: int optional, default is None
        The used reduced k-mer size, if undefined automatically set to just above half the k-mer size
        
    onlySplittingKmers: bool, optional, default is False
        Only the first step will be executed, using a single process.
    
    debug: bool, optional, default is False
        Only use this when debugging or extending the code.      
        
    keepTemporaryFiles: bool, optional, default is False
        Only use this when debugging or extending the code.      
        
    """

    version = "20220916"
    
    def __init__(self,
                 k: int, 
                 name: str, 
                 filenameBase: str,
                 sortedIndexFile: str,
                 readFiles=[], pairedReadFiles=[],
                 minimumFrequency: int = 2,
                 maximumProcesses: int = 5,
                 maximumProcessesAutomaton: int = 1,
                 maximumProcessesIndex: int = None,
                 maximumProcessesMatches: int = None,
                 maximumProcessesConnections: int = None,
                 maximumProcessesMerges: int = None,
                 automatonKmerSize: int = None,
                 onlySplittingKmers: bool = False,
                 debug: bool = False,
                 keepTemporaryFiles: bool=False):  
        
        
        
        """
        Internal use only: initialize
        """
        
        if mp.get_start_method()=="spawn":
            frame = sys._getframe()
            while frame:
                if "__name__" in frame.f_locals.keys():
                    if not frame.f_locals["__name__"]=="__main__":
                        return
                    break                    
                frame = frame.f_back
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("create data storage for "+str(name)+" in "+filenameBase)
        
        #store variables
        self.k=k
        self.automatonKmerSize=math.ceil((self.k+1)/2) if automatonKmerSize==None else math.min(self.k,automatonKmerSize)
        self.name=name
        self.minimumFrequency = minimumFrequency
        self.debug = debug
        self.keepTemporaryFiles = keepTemporaryFiles
        self.filenameBase = filenameBase
        self.maximumProcesses = maximumProcesses
        self.maximumProcessesAutomaton = (maximumProcesses if maximumProcessesAutomaton==None 
                                          else maximumProcessesAutomaton)
        self.maximumProcessesIndex = (maximumProcesses if maximumProcessesIndex==None 
                                      else maximumProcessesIndex)
        self.maximumProcessesMatches = (maximumProcesses if maximumProcessesMatches==None 
                                        else maximumProcessesMatches)
        self.maximumProcessesConnections = (maximumProcesses if maximumProcessesConnections==None 
                                            else maximumProcessesConnections)
        self.maximumProcessesMerges = (maximumProcesses if maximumProcessesMerges==None 
                                       else maximumProcessesMerges)
        
        #check boundaries number of processes
        assert self.maximumProcesses>=5
        assert self.maximumProcessesAutomaton>=1
        assert self.maximumProcessesMatches>=1
        assert self.maximumProcessesIndex>=1
        assert self.maximumProcessesConnections>=1
        assert self.maximumProcessesMerges>=1
        
        if (not onlySplittingKmers) and (len(readFiles)==0) and (len(pairedReadFiles)==0):
            self._logger.error("no read files provided")
        else:                
            #define filenames
            filename = filenameBase+".h5"            
            
            #skip finished steps
            if self.debug:
                filename_splits = filenameBase+".splits.h5"
                filename_distances = filenameBase+".distances.h5"
                if os.path.exists(filename_distances):
                    shutil.copyfile(filename_distances, filename)
                elif os.path.exists(filename_splits):
                    shutil.copyfile(filename_splits, filename)

            #use tables for temporary file, and h5py for final
            with h5py.File(filename,"a") as h5file:

                #set or check config
                if not "/config" in h5file:
                    h5file.create_group("/config")
                    h5file["/config"].attrs["k"] = self.k
                    h5file["/config"].attrs["version"] = self.version
                    h5file["/config"].attrs["automatonKmerSize"] = self.automatonKmerSize
                    h5file["/config"].attrs["name"] = self.name
                    h5file["/config"].attrs["debug"] = self.debug
                    h5file["/config"].attrs["minimumFrequency"] = self.minimumFrequency
                    h5file.flush()
                else:
                    assert h5file["/config"].attrs["k"] == self.k
                    assert h5file["/config"].attrs["version"] == self.version
                    assert h5file["/config"].attrs["automatonKmerSize"] == self.automatonKmerSize
                    assert h5file["/config"].attrs["name"] == self.name
                    assert h5file["/config"].attrs["debug"] == self.debug
                    assert h5file["/config"].attrs["minimumFrequency"] == self.minimumFrequency
                    
                #these settings are allowed to change in secondary runs
                h5file["/config"].attrs["maximumProcesses"] = self.maximumProcesses
                h5file["/config"].attrs["maximumProcessesAutomaton"] = self.maximumProcessesAutomaton
                h5file["/config"].attrs["maximumProcessesIndex"] = self.maximumProcessesIndex
                h5file["/config"].attrs["maximumProcessesMatches"] = self.maximumProcessesMatches
                h5file["/config"].attrs["maximumProcessesConnections"] = self.maximumProcessesConnections
                h5file["/config"].attrs["maximumProcessesMerges"] = self.maximumProcessesMerges

                #get splitting k-mers from index   
                if not ("/split" in h5file and "/histogram" in h5file):
                    if not os.path.exists(sortedIndexFile):
                        self._logger.error("no sorted k-mer list provided")
                    else:
                        self._logger.debug("get splitting k-mers from the provided index")
                        haplotyping.index.splits.Splits(sortedIndexFile, h5file, 
                                                        self.filenameBase, self.debug, self.keepTemporaryFiles)    
                        h5file.flush()
                        #backup
                        if self.debug:
                            shutil.copyfile(filename, filename_splits)
                else:
                    self._logger.debug("detected splitting k-mers from previous run")
                    
                #parse read files and store distances
                if (not onlySplittingKmers) and ("/split" in h5file) and ("/histogram" in h5file):
                    if not ("/relations" in h5file and "/connections" in h5file):
                        self._logger.debug("parse read files and store distances in database")
                        haplotyping.index.reads.Reads(readFiles,pairedReadFiles, h5file, 
                                                      self.filenameBase, self.debug, self.keepTemporaryFiles)
                        h5file.flush()
                        #backup
                        if self.debug:
                            shutil.copyfile(filename, filename_distances)
                    else:
                        self._logger.debug("detected distances from previous run")

                #finished
                h5file.flush()
                
                if self.debug and not self.keepTemporaryFiles:
                    if os.path.exists(filename_distances):
                        os.remove(filename_distances)
                    if os.path.exists(filename_splits):
                        os.remove(filename_splits)

                
    letters = ["A","C","G","T"]
        
    def detectReadFiles(location: str, recursive=True):
        unpairedReadFiles = []
        pairedReadFiles = []
        #get files
        patterns = ["*.fq.gz","*.fastq.gz"]
        if recursive:
            allReadFiles = [y for x in os.walk(location) 
                for p in patterns
                for y in glob.glob(os.path.join(x[0], p))]
        else:
            allReadFiles = [y for p in patterns
                for y in glob.glob(os.path.join(location, p))]
        allReadFiles.sort()
        #split between regular and paired
        processed = set()
        for filename in allReadFiles:
            if filename in processed:
                pass
            else:
                basename = os.path.basename(filename)
                dirname = filename[:len(filename)-len(basename)]  
                if re.search(r'_R1_001\.', basename):
                    filename0 = dirname + basename
                    filename1 = dirname + re.sub(r'_R1_001\.', 
                                                 '_R2_001.', basename)
                    if filename0==filename and filename1 in allReadFiles and not filename in processed:
                        processed.add(filename0)
                        processed.add(filename1)
                        pairedReadFiles.append((filename0,filename1,))
                    else:
                        processed.add(filename)
                        unpairedReadFiles.append(filename)
                elif re.search(r'_1P\.', basename):
                    filename0 = dirname + basename
                    filename1 = dirname + re.sub(r'_1P\.', 
                                                 '_2P.', basename)
                    if filename0==filename and filename1 in allReadFiles and not filename in processed:
                        processed.add(filename0)
                        processed.add(filename1)
                        pairedReadFiles.append((filename0,filename1,))
                    else:
                        processed.add(filename)
                        unpairedReadFiles.append(filename)
                elif re.search(r'_R1\.', basename):
                    filename0 = dirname + basename
                    filename1 = dirname + re.sub(r'_R1\.', 
                                                 '_R2.', basename)
                    if filename0==filename and filename1 in allReadFiles and not filename in processed:
                        processed.add(filename0)
                        processed.add(filename1)
                        pairedReadFiles.append((filename0,filename1,))
                    else:
                        processed.add(filename)
                        unpairedReadFiles.append(filename)
                else:
                    processed.add(filename)
                    unpairedReadFiles.append(filename)

        return (unpairedReadFiles, pairedReadFiles, allReadFiles)
    
    def detectKmerSize(location: str):
        try:
            with gzip.open(location, "rt") as f: 
                reader = csv.reader(f, delimiter="\t")
                for line in reader:
                    return len(str(line[0]))
        except OSError as ex:
            raise Exception("problem with sorted list: {}".format(ex))
    
    def getTablesUint(maximumValue, position):
        if maximumValue<=np.iinfo(np.uint8).max:
            return tables.UInt8Col(pos=position)
        elif maximumValue<=np.iinfo(np.uint16).max:
            return tables.UInt16Col(pos=position)
        elif maximumValue<=np.iinfo(np.uint32).max:
            return tables.UInt32Col(pos=position)
        else:
            return tables.UInt64Col(pos=position)
        
    def getTablesUintAtom(maximumValue):
        if maximumValue<=np.iinfo(np.uint8).max:
            return tables.UInt8Atom()
        elif maximumValue<=np.iinfo(np.uint16).max:
            return tables.UInt16Atom()
        elif maximumValue<=np.iinfo(np.uint32).max:
            return tables.UInt32Atom()
        else:
            return tables.UInt64Atom()
        
    def getUint(maximumValue):
        if maximumValue<=np.iinfo(np.uint8).max:
            return "uint8"
        elif maximumValue<=np.iinfo(np.uint16).max:
            return "uint16"
        elif maximumValue<=np.iinfo(np.uint32).max:
            return "uint32"
        else:
            return "uint64"

