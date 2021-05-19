import logging, h5py, tables, gzip, csv, os, time
import ahocorasick
import numpy as np

import haplotyping
import haplotyping.index.splits
import haplotyping.index.relations

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
        Base name and location for the splitting k-mer database (and the temporary databases if in debug mode)

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
        
    debug: bool, optional, default is False
        Only use this when debugging or extending the code.
        Running in debug mode will create temporary files using the location defined by filenameBase.
        These temporary files can be used for further inspection and will not be deleted afterwards.         
        
    """
        
    
    def __init__(self, k: int, name: str, 
                 filenameBase: str,
                 sortedIndexFile: str,
                 readFiles=[],pairedReadFiles=[],
                 minimumFrequency: int = 2,
                 debug: bool = False):     
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("create data storage for "+str(name)+" in "+filenameBase)
        
        #store variables
        self.k=k
        self.name=name
        self.minimumFrequency = minimumFrequency
        self.debug = debug
        self.filenameBase = filenameBase
        
        if len(readFiles)==0 and len(pairedReadFiles)==0:
            self._logger.error("no read files provided")
        elif not os.path.exists(sortedIndexFile):
            self._logger.error("no sorted k-mer list provided")
        else:                
            #define filename
            filename = filenameBase+".h5"

            if os.path.exists(filename):
                os.remove(filename)

            #use tables for temporary file, and h5py for final
            with h5py.File(filename,"a") as h5file:

                #config
                if not "/config" in h5file:
                    h5file.create_group("/config")
                h5file["/config"].attrs["k"] = self.k
                h5file["/config"].attrs["name"] = self.name
                h5file["/config"].attrs["debug"] = self.debug
                h5file["/config"].attrs["minimumFrequency"] = self.minimumFrequency
                h5file.flush()

                #get splitting k-mers from index
                self._logger.debug("get splitting k-mers from the provided index")
                haplotyping.index.splits.Splits(sortedIndexFile, h5file, self.filenameBase, self.debug)   

                #parse read files
                self._logger.debug("parse read files and store results in database")
                haplotyping.index.relations.Relations(readFiles,pairedReadFiles, h5file, self.filenameBase, self.debug)   
                
                #finished
                h5file.flush()

