import logging, h5py, tables, gzip, csv
import os, sys, tempfile, numpy as np
import re, haplotyping

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
        self.matchPattern = re.compile(r"^["+"".join(haplotyping.index.Database.letters)+"]+$")
        
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
        
        #create temporary storage
        tableCkmerDef = {
            "ckmer": tables.StringCol(self.k,pos=0),
            "type": tables.StringCol(1,pos=1),
            "number": tables.UInt32Col(pos=2,dflt=0),
        }
        self.tableDumpKmers = self.pytables_storage.create_table(self.pytables_storage.root,
                                                                   "dumpCkmer",tableCkmerDef, 
                                                "Temporary to store dump canonical k-mers")
        
        #administration
        rightSplitBases = 0
        rightSplitKmers = 0
        
        def saveToDumpStorage(base, stored,rightSplitBases,rightSplitKmers):   
            rightSplitBases+=1
            for key,value in stored.items():
                kmer = base+key
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
            return (rightSplitBases,rightSplitKmers)

        #loop over sorted list with k-mers, detect right splitting k-mers
        try:
            with gzip.open(filename, "rt") as f: 
                self._logger.debug("parse sorted list to detect right splitting k-mers")
                reader = csv.reader(f, delimiter="\t")
                previousBase = ""
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
                        currentBase = line[0][:-1]
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
                        #check pattern    
                        elif not self.matchPattern.match(currentKmer):    
                            if errorNumber==0:
                                self._logger.warning("k-mer didn't match pattern")
                            errorNumber+=1
                            pass
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
                                if currentBase==previousBase:
                                    if currentBranch==previousBranch:
                                        self._logger.warning("detected reoccurrence of same k-mer ("+currentKmer+")")
                                    else:
                                        if not previousBranch in stored.keys():
                                            stored[previousBranch]=previousNumber
                                        stored[currentBranch]=currentNumber
                                else:
                                    if len(stored)>0:
                                        (rightSplitBases,rightSplitKmers) = saveToDumpStorage(previousBase,stored,
                                                                                              rightSplitBases,rightSplitKmers)
                                        stored={}
                                previousBase = currentBase
                                previousBranch = currentBranch
                                previousKmer = currentKmer
                                previousNumber = currentNumber
                if len(stored)>0:
                    (rightSplitBases,rightSplitKmers) = saveToDumpStorage(previousBase,stored,
                                                                          rightSplitBases,rightSplitKmers)  
                self.tableDumpKmers.flush()
                #warning
                if errorNumber>0:
                    self._logger.warning("skipped "+str(errorNumber)+" items in sorted list")
                #stats
                self._logger.debug("found "+str(rightSplitBases)+" rightSplitBases")
                self._logger.debug("found "+str(rightSplitKmers)+" rightSplitKmers")
                self.h5file["/config/"].attrs["rightSplitBases"]=rightSplitBases
                self.h5file["/config/"].attrs["rightSplitKmers"]=rightSplitKmers
        except OSError as ex:
            self._logger.error("problem with sorted list: "+str(ex))


    def _sort(self):                        
        
        def saveToSortedCkmerStorage(ckmer,ckmerType,number,numberOfCkmers):
            ckmerRow = self.tableSortedKmers.row            
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
            baseRow = self.tableSortedBases.row    
            baseRow["base"] = base
            baseRow["number"] = number
            for letter in branches.keys():
                baseRow["branches/"+str(letter)+"/number"] = branches[letter]["number"]
                baseRow["branches/"+str(letter)+"/ckmerLink"] = branches[letter]["ckmerLink"]
                ckmerRow = self.tableSortedKmers[branches[letter]["ckmerLink"]]
                ckmer = ckmerRow["ckmer"].decode()
                rkmer=haplotyping.General.reverse_complement(ckmer)
                if ckmer == base+letter:
                    self.tableSortedKmers.modify_column(branches[letter]["ckmerLink"],
                                                        branches[letter]["ckmerLink"]+1,
                                                        column=[numberOfBases],
                                                        colname="rightSplitBaseLink/rightSplit")   
                if rkmer == base+letter:
                    self.tableSortedKmers.modify_column(branches[letter]["ckmerLink"],
                                                        branches[letter]["ckmerLink"]+1,
                                                        column=[numberOfBases],
                                                        colname="rightSplitBaseLink/leftSplit")   
            baseRow.append()
            return numberOfBases+1
            
        def saveToDumpBaseStorage(base,branch,number,ckmerLink):   
            baseRow = self.tableDumpRightSplitBases.row
            baseRow["base"] = base
            baseRow["branch"] = branch
            baseRow["number"] = number
            baseRow["ckmerLink"] = ckmerLink
            baseRow.append()
            
        #create sorted and grouped storage
        self.tableDumpKmers.flush()                
        numberOfSplittingKmers=self.tableDumpKmers.shape[0]
        if numberOfSplittingKmers>0:
            
            #create temporary storage
            tableCkmerDef = {
                "ckmer": tables.StringCol(self.k,pos=0),
                "type": tables.StringCol(1,pos=1),
                "rightSplitBaseLink": {
                    "leftSplit": self._getTablesUint(numberOfSplittingKmers,0),
                    "rightSplit": self._getTablesUint(numberOfSplittingKmers,1),
                },
                "number": tables.UInt32Col(pos=3,dflt=0),
            }
            self.tableSortedKmers = self.pytables_storage.create_table(self.pytables_storage.root,
                                                                       "sortedCkmer",tableCkmerDef, 
                                                    "Temporary to store sorted canonical k-mers")  
            
            #create dump table for bases
            tableDumpBaseDef = {
                "base": tables.StringCol(self.k-1,pos=0),
                "branch": tables.StringCol(1,pos=1),
                "number": tables.UInt32Col(pos=2,dflt=0),
                "ckmerLink": self._getTablesUint(numberOfSplittingKmers,3),
            }
            self.tableDumpRightSplitBases = self.pytables_storage.create_table(self.pytables_storage.root,
                                                                   "dumpRightSplitBase",tableDumpBaseDef, 
                                                                   "Temporary to store dump bases")
            
            #store sorted k-mers
            self._logger.debug("sort and group "+str(numberOfSplittingKmers)+" splitting k-mers")
            self.tableDumpKmers.cols.ckmer.create_csindex()
            self.tableDumpKmers.flush()
            previousCkmer=None
            previousType=None
            previousN=None
            self.maxNumber = 0
            numberOfCkmers = 0
            for row in self.tableDumpKmers.itersorted("ckmer",checkCSI=True):
                currentCkmer=row["ckmer"]
                currentType=row["type"]
                currentN=row["number"]
                if previousCkmer and not previousCkmer==currentCkmer:
                    numberOfCkmers = saveToSortedCkmerStorage(previousCkmer.decode(),previousType,
                                                         previousN,numberOfCkmers)
                    self.maxNumber=max(self.maxNumber,previousN)
                    previousType=currentType.decode()
                elif previousCkmer and not previousType==currentType:
                    previousType="b"
                else:
                    previousType=currentType.decode()
                previousCkmer=currentCkmer
                previousN=currentN
            if previousCkmer:        
                numberOfCkmers = saveToSortedCkmerStorage(previousCkmer.decode(),previousType,
                                                     previousN,numberOfCkmers) 
                self.maxNumber=max(self.maxNumber,previousN)
                
            self.tableDumpRightSplitBases.flush()                
            numberOfBases=self.tableDumpRightSplitBases.shape[0]
                
            #store sorted bases
            self.tableSortedKmers.flush()
            numberOfKmers=self.tableSortedKmers.shape[0]
            tableBaseDef = {
                "base": tables.StringCol(self.k-1,pos=0),
                "number": tables.UInt32Col(pos=1,dflt=0),    
                "branches": {},                            
            }
            for letter in haplotyping.index.Database.letters:
                tableBaseDef["branches"][letter] = {
                    "number": tables.UInt32Col(pos=0,dflt=0),
                    "ckmerLink": self._getTablesUint(numberOfKmers,1),
                }
            self.tableSortedBases = self.pytables_storage.create_table(self.pytables_storage.root,
                                                                   "sortedBase",tableBaseDef, 
                                                                   "Temporary to store dump bases")
            
            #store sorted bases
            self._logger.debug("sort and group "+str(numberOfBases)+" bases")
            self.tableDumpRightSplitBases.cols.base.create_csindex()
            self.tableDumpRightSplitBases.flush()
            previousBase=None
            previousBranches={}
            previousNumber=0
            numberOfBases = 0
            for row in self.tableDumpRightSplitBases.itersorted("base",checkCSI=True):
                currentBase=row["base"]
                currentBranch=row["branch"]
                currentNumber=row["number"]
                currentCkmerLink=row["ckmerLink"]
                if previousBase and not previousBase==currentBase:
                    numberOfBases = saveToSortedBaseStorage(previousBase.decode(),previousBranches,
                                                            previousNumber,numberOfBases)
                    previousBase=None
                    previousBranches={}
                    previousNumber=0
                previousBase=currentBase
                previousBranches[row["branch"].decode()] = {"number": currentNumber, "ckmerLink": currentCkmerLink}
                previousNumber+=currentNumber
            if previousBase:
                    numberOfBases = saveToSortedBaseStorage(previousBase.decode(),previousBranches,
                                                            previousNumber,numberOfBases)    
                
        else:
            self._logger.warning("no splitting k-mers to sort and group")
            
    def _store(self):
        canonicalSplitKmers = 0
        canonicalSplitKmersLeft = 0
        canonicalSplitKmersRight = 0
        canonicalSplitKmersBoth = 0
        #create storage    
        self.tableSortedKmers.flush()
        self.tableSortedBases.flush()
        numberOfKmers=self.tableSortedKmers.shape[0]
        numberOfBases=self.tableSortedBases.shape[0]
        #don't make the structure unnecessary big
        dtypeCkmerList=[("ckmer","S"+str(self.k)),
                   ("type","S1"),
                   ("number",self._getUint(self.maxNumber)),
                   ("rightSplitBaseLink",[("leftSplit",self._getUint(numberOfBases)),
                                          ("rightSplit",self._getUint(numberOfBases))]),
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
        dtCkmer=np.dtype(dtypeCkmerList)
        dsCkmer=self.h5file["/split/"].create_dataset("ckmer",(numberOfKmers,), dtype=dtCkmer, chunks=None)
        self._logger.info("store "+str(numberOfKmers)+" splitting k-mers")
        #add stored and grouped kmers to the final unchunked storage
        for i in range(0,numberOfKmers,self.stepSizeStorage):
            stepData = self.tableSortedKmers[i:i+self.stepSizeStorage]
            dsCkmer[i:i+self.stepSizeStorage] = stepData
            for row in stepData:
                canonicalSplitKmers+=1
                if row[1].decode()=="l":
                    canonicalSplitKmersLeft+=1
                elif row[1].decode()=="b":
                    canonicalSplitKmersBoth+=1
                elif row[1].decode()=="r":
                    canonicalSplitKmersRight+=1
        #don't make the structure unnecessary big
        dtypeBaseList=[("base","S"+str(self.k-1)),
                   ("number",self._getUint(self.maxNumber)),
                   ("branches",[])]
        for i in range(len(haplotyping.index.Database.letters)):
            letter = haplotyping.index.Database.letters[i]
            dtypeBaseList[2][1].append((letter,[("number",self._getUint(self.maxNumber)),
                                                ("ckmerLink",numberOfKmers)]))
        dtBase=np.dtype(dtypeBaseList)
        dsBase=self.h5file["/split/"].create_dataset("base",(numberOfBases,), dtype=dtBase, chunks=None)
        self._logger.info("store "+str(numberOfBases)+" bases")
        #add stored and grouped bases to the final unchunked storage
        for i in range(0,numberOfBases,self.stepSizeStorage):
            stepData = self.tableSortedBases[i:i+self.stepSizeStorage]
            dsBase[i:i+self.stepSizeStorage] = stepData
        #store the stats
        self._logger.debug("found "+str(canonicalSplitKmersLeft)+" canonicalSplitKmersLeft")
        self._logger.debug("found "+str(canonicalSplitKmersRight)+" canonicalSplitKmersRight")
        self._logger.debug("found "+str(canonicalSplitKmersBoth)+" canonicalSplitKmersBoth")
        self.h5file["/config/"].attrs["canonicalSplitKmers"]=canonicalSplitKmers
        self.h5file["/config/"].attrs["canonicalSplitKmersLeft"]=canonicalSplitKmersLeft
        self.h5file["/config/"].attrs["canonicalSplitKmersBoth"]=canonicalSplitKmersBoth
        self.h5file["/config/"].attrs["canonicalSplitKmersRight"]=canonicalSplitKmersRight 
        
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

    