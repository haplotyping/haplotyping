import logging, h5py, os, glob, re, math

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
                 automatonKmerSize: int = None,
                 debug: bool = False):     
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("create data storage for "+str(name)+" in "+filenameBase)
        
        #store variables
        self.k=k
        self.automatonKmerSize=math.ceil((self.k+1)/2) if automatonKmerSize==None else math.min(self.k,automatonKmerSize)
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

                #set or check config
                if not "/config" in h5file:
                    h5file.create_group("/config")
                    h5file["/config"].attrs["k"] = self.k
                    h5file["/config"].attrs["automatonKmerSize"] = self.automatonKmerSize
                    h5file["/config"].attrs["name"] = self.name
                    h5file["/config"].attrs["debug"] = self.debug
                    h5file["/config"].attrs["minimumFrequency"] = self.minimumFrequency
                    h5file.flush()
                else:
                    assert h5file["/config"].attrs["k"] == self.k
                    assert h5file["/config"].attrs["automatonKmerSize"] == self.automatonKmerSize
                    assert h5file["/config"].attrs["name"] == self.name
                    assert h5file["/config"].attrs["debug"] == self.debug
                    assert h5file["/config"].attrs["minimumFrequency"] == self.minimumFrequency         
                    
                #get splitting k-mers from index                
                self._logger.debug("get splitting k-mers from the provided index")
                haplotyping.index.splits.Splits(sortedIndexFile, h5file, self.filenameBase, self.debug)                         
                                                
                #parse read files if splitting k-mers were found
                if h5file["/config/"].attrs["canonicalSplitKmersBoth"]>0:
                    self._logger.debug("parse read files and store results in database")
                    haplotyping.index.relations.Relations(readFiles,pairedReadFiles, h5file, self.filenameBase, self.debug)  
                else:
                    self._logger.error("no splitting k-mers were found")
                
                #finished
                h5file.flush()
                
                
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
                if re.search(r'([^a-zA-Z0-9])R1([^a-zA-Z0-9])', basename):
                    filename0 = dirname + basename
                    filename1 = dirname + re.sub(r'([^a-zA-Z0-9])R1([^a-zA-Z0-9])', 
                                                 '\\1R2\\2', basename)
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

