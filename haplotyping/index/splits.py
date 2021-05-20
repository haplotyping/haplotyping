import logging, h5py, tables, gzip, csv
import os, sys, tempfile, numpy as np
import haplotyping

class Splits:
    
    """
    Internal use, get splitting k-mers from index
    """
    
    def __init__(self, sortedIndexFile: str, h5file, filenameBase, debug=False):
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("parse k-mers to find right splitting k-mers")
        
        self.stepSizeStorage=10000
        
        #set variables
        self.k = h5file["/config"].attrs["k"]
        self.minimumFrequency = h5file["/config"].attrs["minimumFrequency"]
        self.h5file = h5file
        self.maxNumber = 0
        self.debug = debug
        
        #check existence group
        if not "/split" in h5file:
            h5file.create_group("/split")
        
        #create splitting k-mer dataset
        if "/split/ckmer" in h5file:
            self._logger.warning("canonical splitting k-mer dataset already in hdf5 storage")
        else:
        
            try:                                                
                if self.debug:
                    self.tmpDirectory = None
                    pytablesFile = filenameBase+"_tmp_split.h5"
                    if os.path.exists(pytablesFile):
                        os.remove(pytablesFile)
                    self._logger.debug("store temporary in "+pytablesFile)    
                else:
                    self.tmpDirectory = tempfile.TemporaryDirectory()
                    pytablesFile = self.tmpDirectory.name+"/kmer.data_tmp_split.h5"
            
                #create datasets
                with tables.open_file(pytablesFile, mode="w", title="Temporary storage") as self.pytables_storage:
                    #create temporary storage
                    tableKmersDef = {
                        "ckmer": tables.StringCol(self.k,pos=0),
                        "type": tables.StringCol(1,pos=1),
                        "number": tables.UInt32Col(pos=2,dflt=0),
                    }
                    self.tableDumpKmers = self.pytables_storage.create_table(self.pytables_storage.root,
                                                                               "dumpCkmer",tableKmersDef, 
                                                            "Temporary to store dump canonical k-mers")
                    self.tableSortedKmers = self.pytables_storage.create_table(self.pytables_storage.root,
                                                                               "sortedCkmer",tableKmersDef, 
                                                            "Temporary to store sorted canonical k-mers")

                    #get k-mers
                    self._parseIndex(sortedIndexFile)
                    self._sort()
                    self._store()
                    #flush
                    self.h5file.flush()
            #except:
            #    self._logger.error("problem occurred while constructing splits")
            finally:
                if not self.debug:
                    self.tmpDirectory.cleanup()

    
    def _parseIndex(self, filename: str):
        #administration
        rightSplitTrunks = 0
        rightSplitKmers = 0
        
        def saveToDumpStorage(trunk, stored,rightSplitTrunks,rightSplitKmers):   
            rightSplitTrunks+=1
            for key,value in stored.items():
                kmer = trunk+key
                ckmer = haplotyping.General.canonical(kmer)
                #assume right splitting if k-mer equals rc
                #todo: can this (in theory) cause problematic situations?
                ckmerType = "r" if kmer==ckmer else "l"
                ckmerRow = self.tableDumpKmers.row            
                ckmerRow["ckmer"] = ckmer
                ckmerRow["type"] = ckmerType
                ckmerRow["number"] = value
                ckmerRow.append()
                rightSplitKmers+=1
            return (rightSplitTrunks,rightSplitKmers)

        #loop over sorted list with k-mers, detect right splitting k-mers
        try:
            with gzip.open(filename, "rt") as f: 
                self._logger.debug("parse sorted list to detect right splitting k-mers")
                reader = csv.reader(f, delimiter="\t")
                previousTrunk = ""
                previousBranch = ""
                previousKmer = ""
                previousNumber = 0            
                stored={}
                errorNumber = 0
                for line in reader:
                    if len(line)<2:
                        if errorNumber==0:
                            self._logger.warning("unexpected item in sorted list: "+str(line))
                        errorNumber+=1
                    else:
                        currentTrunk = line[0][:-1]
                        currentBranch = line[0][-1]
                        currentKmer = line[0]
                        #always get number
                        try:
                            currentNumber = int(line[1])
                        except ValueError:
                            currentNumber = 0
                        #check k-mer size
                        if not len(currentKmer)==self.k:
                            if errorNumber==0:
                                self._logger.warning("different k-mer size in sorted list ("+
                                                   str(len(currentKmer))+" instead of "+str(self.k)+")")
                            errorNumber+=1
                        #check sorted
                        elif not (previousKmer=="") and not (previousKmer<currentKmer):
                            if errorNumber==0:
                                self._logger.warning("provided list not properly sorted")
                            errorNumber+=1
                        #check frequency
                        elif currentNumber<1:
                            if errorNumber==0:
                                self._logger.warning("unexpected frequency for "+str(currentKmer))
                            errorNumber+=1
                        else:
                            if currentNumber < self.minimumFrequency:
                                pass
                            else:
                                if currentTrunk==previousTrunk:
                                    if currentBranch==previousBranch:
                                        self._logger.warning("detected reoccurrence of same k-mer ("+currentKmer+")")
                                    else:
                                        if not previousBranch in stored.keys():
                                            stored[previousBranch]=previousNumber
                                        stored[currentBranch]=currentNumber
                                else:
                                    if len(stored)>0:
                                        (rightSplitTrunks,rightSplitKmers) = saveToDumpStorage(previousTrunk,stored,
                                                                                              rightSplitTrunks,rightSplitKmers)
                                        stored={}
                                previousTrunk = currentTrunk
                                previousBranch = currentBranch
                                previousKmer = currentKmer
                                previousNumber = currentNumber
                if len(stored)>0:
                    (rightSplitTrunks,rightSplitKmers) = saveToDumpStorage(previousTrunk,stored,
                                                                          rightSplitTrunks,rightSplitKmers)  
                self.tableDumpKmers.flush()
                #warning
                if errorNumber>0:
                    self._logger.warning("skipped "+str(errorNumber)+" items in sorted list")
                #stats
                self._logger.debug("found "+str(rightSplitTrunks)+" rightSplitTrunks")
                self._logger.debug("found "+str(rightSplitKmers)+" rightSplitKmers")
                self.h5file["/config/"].attrs["rightSplitTrunks"]=rightSplitTrunks
                self.h5file["/config/"].attrs["rightSplitKmers"]=rightSplitKmers
        except OSError as ex:
            self._logger.error("problem with sorted list: "+str(ex))


    def _sort(self):
        def saveToSortedStorage(ckmer,ckmerType,value):
            ckmerRow = self.tableSortedKmers.row            
            ckmerRow["ckmer"] = ckmer
            ckmerRow["type"] = ckmerType
            ckmerRow["number"] = value
            ckmerRow.append()
        #create sorted and grouped storage
        self.tableDumpKmers.flush()                
        numberOfSplittingKmers=self.tableDumpKmers.shape[0]
        if numberOfSplittingKmers>0:
            self._logger.debug("sort and group "+str(numberOfSplittingKmers)+" splitting k-mers")
            self.tableDumpKmers.cols.ckmer.create_csindex()
            self.tableDumpKmers.flush()
            previousCkmer=None
            previousType=None
            previousN=None
            self.maxNumber = 0
            for row in self.tableDumpKmers.itersorted("ckmer",checkCSI=True):
                currentCkmer=row["ckmer"]
                currentType=row["type"]
                currentN=row["number"]
                if previousCkmer and not previousCkmer==currentCkmer:
                    saveToSortedStorage(previousCkmer,previousType,previousN)
                    self.maxNumber=max(self.maxNumber,previousN)
                    previousType=currentType
                elif previousCkmer and not previousType==currentType:
                    previousType="b"
                else:
                    previousType=currentType
                previousCkmer=currentCkmer
                previousN=currentN
            if previousCkmer:        
                saveToSortedStorage(previousCkmer,previousType,previousN) 
                self.maxNumber=max(self.maxNumber,previousN)
        else:
            self._logger.warning("no splitting k-mers to sort and group")
            
    def _store(self):
        canonicalSplitKmers = 0
        canonicalSplitKmersLeft = 0
        canonicalSplitKmersRight = 0
        canonicalSplitKmersBoth = 0
        #create storage    
        self.tableSortedKmers.flush()
        numberOfKmers=self.tableSortedKmers.shape[0]
        #don't make the structure unnecessary big
        dtypeList=[("ckmer","S"+str(self.k)),
                   ("type","S1"),
                   ("number",self._getUint(self.maxNumber)),
                   ("direct",[("left",
                               [("distinct","uint8"),
                                ("number",self._getUint(self.maxNumber))]),
                              ("right",
                               [("distinct","uint8"),
                                ("number",self._getUint(self.maxNumber))])
                             ]),
                   ("connected",[("distinct",self._getUint(self.maxNumber))]),
                   ("cycle",[("number",self._getUint(self.maxNumber))]),
                   ("reversal",[("number",self._getUint(self.maxNumber))])]
        dt=np.dtype(dtypeList)
        ds=self.h5file["/split/"].create_dataset("ckmer",(numberOfKmers,), dtype=dt, chunks=None)
        self._logger.info("store "+str(numberOfKmers)+" splitting k-mers")
        #add stored and grouped kmers to the final unchunked storage
        for i in range(0,numberOfKmers,self.stepSizeStorage):
            stepData = self.tableSortedKmers[i:i+self.stepSizeStorage]
            ds[i:i+self.stepSizeStorage] = stepData
            for row in stepData:
                canonicalSplitKmers+=1
                if row[1].decode()=="l":
                    canonicalSplitKmersLeft+=1
                elif row[1].decode()=="b":
                    canonicalSplitKmersBoth+=1
                elif row[1].decode()=="r":
                    canonicalSplitKmersRight+=1
        #store the stats
        self._logger.debug("found "+str(canonicalSplitKmersLeft)+" canonicalSplitKmersLeft")
        self._logger.debug("found "+str(canonicalSplitKmersRight)+" canonicalSplitKmersRight")
        self._logger.debug("found "+str(canonicalSplitKmersBoth)+" canonicalSplitKmersBoth")
        self.h5file["/config/"].attrs["canonicalSplitKmers"]=canonicalSplitKmers
        self.h5file["/config/"].attrs["canonicalSplitKmersLeft"]=canonicalSplitKmersLeft
        self.h5file["/config/"].attrs["canonicalSplitKmersBoth"]=canonicalSplitKmersBoth
        self.h5file["/config/"].attrs["canonicalSplitKmersRight"]=canonicalSplitKmersRight 
        
    def _getUint(self, maximumValue):
        if maximumValue<=np.iinfo(np.uint8).max:
            return "uint8"
        elif maximumValue<=np.iinfo(np.uint16).max:
            return "uint16"
        elif maximumValue<=np.iinfo(np.uint32).max:
            return "uint32"
        else:
            return "uint64"

    