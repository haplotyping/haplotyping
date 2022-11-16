import logging, h5py, tables, gzip, csv
import os, sys, numpy as np
import re, haplotyping, ahocorasick, pickle
import haplotyping.index.database

class Splits:
    
    """
    Internal use, get splitting k-mers from index
    """
    
    stepSizeStorage = 10000
    
    def __init__(self, sortedIndexFile: str, h5file, filenameBase, debug=False, keepTemporaryFiles=False):
        
        """
        Internal use only: initialize
        """
        
        #logger
        self._logger = logging.getLogger(__name__)
        self._logger.info("parse k-mers to find right splitting k-mers")
        
        self.matchPattern = re.compile(r"^["+"".join(haplotyping.index.Database.letters)+"]+$")
        
        #set variables
        self.k = h5file["/config"].attrs["k"]
        self.automatonKmerSize = h5file["/config"].attrs["automatonKmerSize"]
        self.minimumFrequency = h5file["/config"].attrs["minimumFrequency"]
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
                                      title="Temporary storage") as self.pytablesStorage:
                    #get k-mers
                    self._parseIndex(sortedIndexFile)
                    self._sort()
                    self._store()                    
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

    def _parseIndex(self, filename: str):
        
        #create temporary storage
        tableCkmerDef = {
            "ckmer": tables.StringCol(self.k,pos=0),
            "type": tables.StringCol(1,pos=1),
            "number": tables.UInt64Col(pos=2),
        }
        self.tableDumpKmers = self.pytablesStorage.create_table(self.pytablesStorage.root,
                                                                   "dumpCkmer",tableCkmerDef, 
                                                "Temporary to store dump canonical k-mers")
        
        #administration
        rightSplitBases = 0
        rightSplitKmers = 0
        
        def saveToDumpStorage(base, stored, rightSplitBases,rightSplitKmers):   
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
                self.frequencyHistogram["ckmer"][value] = (
                                self.frequencyHistogram["ckmer"].get(value,0)+1)
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
                totalNumberOfKmers = 0
                totalSumOfKmerFrequencies = 0
                minimumAllKmerFrequencies = 2**32
                maximumAllKmerFrequencies = 0
                self.frequencyHistogram = {"kmer": {}, "ckmer": {}, "base": {}}
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
                        totalNumberOfKmers+=1
                        #always get number
                        try:
                            currentNumber = int(line[1])
                            totalSumOfKmerFrequencies += currentNumber
                            minimumAllKmerFrequencies = min(minimumAllKmerFrequencies,currentNumber)
                            maximumAllKmerFrequencies = max(minimumAllKmerFrequencies,currentNumber)
                            self.frequencyHistogram["kmer"][currentNumber] = (
                                self.frequencyHistogram["kmer"].get(currentNumber,0)+1)
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
                                        if (rightSplitBases%1000000)==0:
                                            self._logger.debug("processed "+str(rightSplitBases)
                                                               +" bases and "+str(rightSplitKmers)+" k-mers")
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
                    self._logger.warning("skipped {} items in sorted list".format(errorNumber))
                #stats
                minimumAllKmerFrequencies=min(minimumAllKmerFrequencies,maximumAllKmerFrequencies)
                self._logger.debug("checked {} k-mers with {} total frequency".format(
                    totalNumberOfKmers,totalSumOfKmerFrequencies))
                self._logger.debug("found {} rightSplitBases and {} rightSplitKmers".format(
                    rightSplitBases,rightSplitKmers))
                self._logger.debug("frequency k-mers between {} and {}".format(
                    minimumAllKmerFrequencies,maximumAllKmerFrequencies))
                self.h5file["/config/"].attrs["rightSplitBases"]=rightSplitBases
                self.h5file["/config/"].attrs["rightSplitKmers"]=rightSplitKmers
                self.h5file["/config/"].attrs["totalNumberOfKmers"]=totalNumberOfKmers
                self.h5file["/config/"].attrs["totalSumOfKmerFrequencies"]=totalSumOfKmerFrequencies
                self.h5file["/config/"].attrs["minimumAllKmerFrequencies"]=minimumAllKmerFrequencies
                self.h5file["/config/"].attrs["maximumAllKmerFrequencies"]=maximumAllKmerFrequencies
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
            self.frequencyHistogram["base"][number] = (
                                self.frequencyHistogram["base"].get(number,0)+1)
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
                    "leftSplit": haplotyping.index.Database.getTablesUint(numberOfSplittingKmers,0),
                    "rightSplit": haplotyping.index.Database.getTablesUint(numberOfSplittingKmers,1),
                },
                "number": tables.UInt32Col(pos=3,dflt=0),
            }
            self.tableSortedKmers = self.pytablesStorage.create_table(self.pytablesStorage.root,
                                                                       "sortedCkmer",tableCkmerDef, 
                                                    "Temporary to store sorted canonical k-mers",
                                                                      expectedrows=numberOfSplittingKmers)  
            
            #create dump table for bases
            tableDumpBaseDef = {
                "base": tables.StringCol(self.k-1,pos=0),
                "branch": tables.StringCol(1,pos=1),
                "number": tables.UInt32Col(pos=2,dflt=0),
                "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfSplittingKmers,3),
            }
            self.tableDumpRightSplitBases = self.pytablesStorage.create_table(self.pytablesStorage.root,
                                                                   "dumpRightSplitBase",tableDumpBaseDef, 
                                                                   "Temporary to store dump bases",
                                                                       expectedrows=numberOfSplittingKmers)
            
            #store sorted k-mers
            self._logger.debug("sort and group "+str(numberOfSplittingKmers)+" splitting k-mers")
            self.tableDumpKmers.cols.ckmer.create_csindex()
            self.tableDumpKmers.flush()
            previousCkmer=None
            previousType=None
            previousN=None
            self.maximumNumber = 0
            numberOfCkmers = 0
            for row in self.tableDumpKmers.itersorted("ckmer",checkCSI=True):
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
            if previousCkmer:        
                numberOfCkmers = saveToSortedCkmerStorage(previousCkmer.decode(),previousType,
                                                     previousN,numberOfCkmers) 
                self.maximumNumber=max(self.maximumNumber,previousN)
                
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
                    "ckmerLink": haplotyping.index.Database.getTablesUint(numberOfKmers,1),
                }
            self.tableSortedBases = self.pytablesStorage.create_table(self.pytablesStorage.root,
                                                                   "sortedBase",tableBaseDef, 
                                                                   "Temporary to store dump bases",
                                                                       expectedrows=numberOfKmers)
            
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
                   ("connected",[("link",haplotyping.index.Database.getUint(numberOfKmers)),
                                 ("number",haplotyping.index.Database.getUint(len(haplotyping.index.Database.letters)**4)),
                                 ("connected",
                                      haplotyping.index.Database.getUint(len(haplotyping.index.Database.letters)**4))]),
                   ("paired",[("link",haplotyping.index.Database.getUint(numberOfKmers)),
                              ("number",haplotyping.index.Database.getUint(len(haplotyping.index.Database.letters)**4))]),
                   ("cycle",[("number",haplotyping.index.Database.getUint(self.maximumNumber))]),
                   ("reversal",[("number",haplotyping.index.Database.getUint(self.maximumNumber))])]
        dtCkmer=np.dtype(dtypeCkmerList)
        dsCkmer=self.h5file["/split/"].create_dataset("ckmer",(numberOfKmers,), dtype=dtCkmer, chunks=None)
        self._logger.info("store {} splitting k-mers".format(numberOfKmers))
        #add stored and grouped kmers to the final unchunked storage
        for i in range(0,numberOfKmers,Splits.stepSizeStorage):
            stepData = self.tableSortedKmers[i:i+Splits.stepSizeStorage]
            dsCkmer[i:i+Splits.stepSizeStorage] = stepData
            for row in stepData:
                canonicalSplitKmers+=1
                if row[1].decode()=="l":
                    canonicalSplitKmersLeft+=1
                elif row[1].decode()=="b":
                    canonicalSplitKmersBoth+=1
                elif row[1].decode()=="r":
                    canonicalSplitKmersRight+=1
        # BASE STORAGE - don't make the structure unnecessary big
        dtypeBaseList=[("base","S"+str(self.k-1)),
                   ("number",haplotyping.index.Database.getUint(self.maximumNumber)),
                   ("branches",[])]
        for i in range(len(haplotyping.index.Database.letters)):
            letter = haplotyping.index.Database.letters[i]
            dtypeBaseList[2][1].append((letter,[("number",haplotyping.index.Database.getUint(self.maximumNumber)),
                                                ("ckmerLink",numberOfKmers)]))
        dtBase=np.dtype(dtypeBaseList)
        dsBase=self.h5file["/split/"].create_dataset("base",(numberOfBases,), dtype=dtBase, chunks=None)
        self._logger.info("store {} bases".format(numberOfBases))
        #add stored and grouped bases to the final unchunked storage
        for i in range(0,numberOfBases,Splits.stepSizeStorage):
            stepData = self.tableSortedBases[i:i+Splits.stepSizeStorage]
            dsBase[i:i+Splits.stepSizeStorage] = stepData
        # HISTOGRAM K-MER STORAGE - don't make the structure unnecessary big
        maximumFrequency = max(self.frequencyHistogram["kmer"].keys())
        maximumNumber = max(self.frequencyHistogram["kmer"].values())
        dtypeFrequencyHistogramKmerList=[
                    ("frequency",haplotyping.index.Database.getUint(maximumFrequency)),
                    ("number",haplotyping.index.Database.getUint(maximumNumber)),]
        dtFrequencyHistogramKmer=np.dtype(dtypeFrequencyHistogramKmerList)
        dsFrequencyHistogramKmer=self.h5file["/histogram/"].create_dataset("kmer",(len(self.frequencyHistogram["kmer"]),), 
                                                                  dtype=dtFrequencyHistogramKmer, chunks=None)
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
                                                                  dtype=dtFrequencyHistogramCkmer, chunks=None)
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
                                                                  dtype=dtFrequencyHistogramBase, chunks=None)
        self._logger.info("store {} entries splitting k-mer bases histogram".format(len(self.frequencyHistogram["base"])))
        #store histogram
        dsFrequencyHistogramBase[0:len(self.frequencyHistogram["base"])] = list(
            sorted(self.frequencyHistogram["base"].items()))
        #store the stats
        self._logger.debug("found {} canonicalSplitKmersLeft".format(canonicalSplitKmersLeft))
        self._logger.debug("found {} canonicalSplitKmersRight".format(canonicalSplitKmersRight))
        self._logger.debug("found {} canonicalSplitKmersBoth".format(canonicalSplitKmersBoth))
        self.h5file["/config/"].attrs["maximumFrequency"]=self.maximumNumber
        self.h5file["/config/"].attrs["canonicalSplitKmers"]=canonicalSplitKmers
        self.h5file["/config/"].attrs["canonicalSplitKmersLeft"]=canonicalSplitKmersLeft
        self.h5file["/config/"].attrs["canonicalSplitKmersBoth"]=canonicalSplitKmersBoth
        self.h5file["/config/"].attrs["canonicalSplitKmersRight"]=canonicalSplitKmersRight 
        
    def createAutomatonWithIndex(h5file, indexFile, automatonFile, k):
        k = min(h5file["/config"].attrs["k"],k)
        logger = logging.getLogger(__name__)
        logger.info("create automaton with k' = {}".format(k))
        automatonSplits = ahocorasick.Automaton()
        numberOfKmers = h5file["/split/ckmer"].shape[0]
        kmers = h5file["/split/ckmer"]
        with open(indexFile, "w") as f:
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
        stats = automatonSplits.get_stats()
        logger.debug("automaton with {} words, size {} MB".format(stats["words_count"],
                                                                  round(stats["total_size"]/1048576)))
        with open(automatonFile, "wb") as f:
            pickle.dump(automatonSplits, f)
        
   

    