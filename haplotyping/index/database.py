import logging, h5py, tables, gzip, csv, os, time
import ahocorasick
import numpy as np

import haplotyping
import haplotyping.index.splits
import haplotyping.index.relations

class Database:
    
    """constructing splitting k-mer database"""
    
    def __init__(self, k: int, name: str, filenameBase: str,
                 sortedIndexFile: str,
                 readFiles,pairedReadFiles,
                 minimumFrequency=2):
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("create data storage for "+str(name)+" in "+filenameBase)
        
        #store variables
        self.k=k
        self.name=name
        self.minimumFrequency = minimumFrequency
        self.filenameBase = filenameBase
        self.processReadsTime = 0
        self.readLengthMinimum=None
        self.readLengthMaximum=None
        self.readUnpairedTotal=0
        self.readPairedTotal=0
        self.readTotal=0        
                
        #define filenames
        filename = filenameBase+".h5"
        
        #if os.path.exists(filename):
        #    os.remove(filename)

        #use tables for temporary file, and h5py for final
        with h5py.File(filename,"a") as h5file:
            
            #config
            if not "/config" in h5file:
                h5file.create_group("/config")
            h5file["/config"].attrs["k"] = self.k
            h5file["/config"].attrs["name"] = self.name
            h5file["/config"].attrs["minimumFrequency"] = self.minimumFrequency
            h5file.flush()
            
            #get splitting k-mers from index
            haplotyping.index.splits.Splits(sortedIndexFile, h5file)   
            
            #parse read files
            haplotyping.index.relations.Relations(readFiles,pairedReadFiles, h5file)   
          
