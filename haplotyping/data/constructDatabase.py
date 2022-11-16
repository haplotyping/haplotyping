#!python

import os,logging,tempfile
from frictionless import Package, extract
import gzip,csv,pandas as pd, numpy as np
import sqlite3, json, time
from sqlite3 import Error
import h5py, random, shutil

class ConstructDatabase:
    
    dbFilename = "db.sqlite"
    identifierBackup = "identifiers.json"
    markerDirectory = "marker"
    kmerDirectory = "kmer"
    
    def __init__(self, basedir:str, exportdir:str):
        #data
        self._baseDir = basedir
        self._exportDir = exportdir        
        #logging
        self._logger = logging.getLogger(__name__)
        #connection
        self._connection = None
        self._cursor = None
        self._clean()
        try:
            self._load_identifiers()
            self._connect_sqlite3()
            self._exportPedigree()
            self._exportResources()
        except Error as e:
            print(e)
        finally:
            self._close_sqlite3()
            self._store_identifiers()
            
    def _load_identifiers(self):
        if os.path.exists(os.path.join(self._exportDir,self.identifierBackup)):
            with open(os.path.join(self._exportDir,self.identifierBackup), mode="r") as f:
                self._identifiers = json.load(f)
                for resource in self._identifiers["collections"]["marker"].keys():
                    for internalId in self._identifiers["collections"]["marker"][resource].keys():
                        self._identifiers["collections"]["marker"][resource][internalId]["active"] = 0
                for resource in self._identifiers["collections"]["kmer"].keys():
                    for internalId in self._identifiers["collections"]["kmer"][resource].keys():
                        self._identifiers["collections"]["kmer"][resource][internalId]["active"] = 0
                for uid in self._identifiers["collections"]["uids"].keys():
                    for id in self._identifiers["collections"]["uids"][uid].keys():
                        self._identifiers["collections"]["uids"][uid][id]["active"] = 0
        else:
            self._identifiers = {"collections": {"uids": {}, "marker": {}, "kmer": {}}, "datasets": []}
            
    def _store_identifiers(self):
        with open(os.path.join(self._exportDir,self.identifierBackup), mode="w") as f:
            json.dump(self._identifiers,f)

    def _connect_sqlite3(self):
        self._logger.info("create {}".format(self.dbFilename))
        self._connection = sqlite3.connect(os.path.join(self._exportDir,self.dbFilename))  
        self._cursor = self._connection.cursor()
        
    def _close_sqlite3(self):
        if self._cursor:
            self._cursor.close()
        if self._connection:            
            self._connection.close()            
            
    def _clean(self):
        if os.path.exists(os.path.join(self._exportDir,self.markerDirectory)):
            self._logger.info("remove previous version {}".format(self.markerDirectory))
            shutil.rmtree(os.path.join(self._exportDir,self.markerDirectory))
        os.mkdir(os.path.join(self._exportDir,self.markerDirectory))
        if not os.path.exists(os.path.join(self._exportDir,self.kmerDirectory)):
            os.mkdir(os.path.join(self._exportDir,self.kmerDirectory))
        if os.path.exists(os.path.join(self._exportDir,self.dbFilename)):
            self._logger.info("remove previous version {}".format(self.dbFilename))
            os.remove(os.path.join(self._exportDir,self.dbFilename))
            
                
    def _exportPedigree(self):  
        #add new
        for filename in sorted(os.listdir(self._baseDir)):
            if os.path.isfile(os.path.join(self._baseDir,filename)):
                if filename=="pedigree.package.json":
                    package = Package(os.path.join(self._baseDir,filename))
                    #create structure
                    self._logger.info("initialise database")
                    sql = """CREATE TABLE IF NOT EXISTS "variety" (
                                "id" INTEGER NOT NULL,
                                "uid" VARCHAR(13) NULL,
                                "name" VARCHAR(255) NULL,
                                "year_min" INTEGER NULL,
                                "year_max" INTEGER NULL,
                                "origin" VARCHAR(4) NULL,
                                "breeder_id" INTEGER NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE INDEX "variety_name" ON "variety" ("name");
                            CREATE INDEX "variety_origin" ON "variety" ("origin");
                            CREATE INDEX "variety_breeder" ON "variety" ("breeder_id");
                            CREATE UNIQUE INDEX "variety_uid" ON "variety" ("uid");
                            CREATE TABLE IF NOT EXISTS "variety_synonym" (
                                "id" INTEGER NOT NULL,
                                "uid" VARCHAR(13) NULL,
                                "synonym" VARCHAR(255) NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE INDEX "synonym_uid" ON "variety_synonym" ("uid");
                            CREATE TABLE IF NOT EXISTS "variety_ancestor" (
                                "id" INTEGER NOT NULL,
                                "variety" VARCHAR(13) NOT NULL,
                                "ancestor" VARCHAR(13) NULL,
                                "type" TEXT NULL,
                                "offspring" INTEGER NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE INDEX "variety_ancestor_ancestor" ON "variety_ancestor" ("ancestor");
                            CREATE INDEX "variety_ancestor_offspring" ON "variety_ancestor" ("offspring");
                            CREATE INDEX "variety_ancestor_variety" ON "variety_ancestor" ("variety");
                            CREATE TABLE IF NOT EXISTS "country" (
                                "id" INTEGER NOT NULL,
                                "uid" VARCHAR(4) NOT NULL,
                                "name" VARCHAR(255) NOT NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE UNIQUE INDEX "country_uid" ON "country" ("uid");
                            CREATE TABLE IF NOT EXISTS "breeder" (
                                "id" INTEGER NOT NULL,
                                "name" VARCHAR(255) NOT NULL,
                                "country" VARCHAR(4) NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE INDEX "breeder_country" ON "breeder" ("country");
                            CREATE TABLE IF NOT EXISTS "collection" (
                                "id" INTEGER NOT NULL,
                                "uid" VARCHAR(13) NULL,
                                "resource" VARCHAR(255) NULL,
                                "internal_id" INTEGER NOT NULL,
                                "type" VARCHAR(6) CHECK("type" IN ("marker","kmer") ) NOT NULL,
                                "name" VARCHAR(255) NULL,
                                "experiment" VARCHAR(255) NULL,
                                "location" VARCHAR(255) NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE UNIQUE INDEX "collection_resource_internal_type" 
                                ON "collection" ("resource","internal_id","type");
                            CREATE UNIQUE INDEX "collection_uid" ON "collection" ("uid");
                            CREATE TABLE IF NOT EXISTS "dataset" (
                                "id" INTEGER NOT NULL,
                                "uid" VARCHAR(13) NULL,
                                "variety" VARCHAR(13) NULL,
                                "collection_id" INTEGER NOT NULL,
                                "internal_id" INTEGER NOT NULL,
                                "location" VARCHAR(255) NULL,
                                "type" VARCHAR(6) CHECK("type" IN ("marker","kmer","split") ) NULL,
                                PRIMARY KEY ("id")
                            );
                            CREATE INDEX "dataset_collection_id" ON "dataset" ("collection_id");
                            CREATE UNIQUE INDEX "dataset_collection_internal" 
                                ON "dataset" ("collection_id","internal_id");
                            CREATE UNIQUE INDEX "dataset_uid" ON "dataset" ("uid");
                            CREATE INDEX "dataset_variety" ON "dataset" ("variety");"""
                    self._connection.executescript(sql)
                    #countries
                    countries = pd.DataFrame(extract(package.get_resource("pedigree_countries")))
                    countries = countries[["code","name"]]
                    countries.columns = ["uid","name"]
                    countries.to_sql("country", con=self._connection, index=False, if_exists="append")
                    self._logger.info("add {} countries".format(len(countries)))
                    #breeders
                    breeders = pd.DataFrame(extract(package.get_resource("pedigree_breeders")))
                    breeders = breeders[["id","name","country"]]
                    breeders.to_sql("breeder", con=self._connection, index=False, if_exists="append")
                    self._logger.info("add {} breeders".format(len(breeders)))
                    #varieties
                    varieties = pd.DataFrame(extract(package.get_resource("pedigree_varieties")))
                    varieties = varieties[["uid", "name", "yearMin", "yearMax", "origin", "breederId"]]
                    varieties.columns = ["uid", "name", "year_min", "year_max", "origin", "breeder_id"]
                    varieties.to_sql("variety", con=self._connection, index=False, if_exists="append")
                    self._logger.info("add {} varieties".format(len(varieties)))
                    #ancestors
                    ancestors = pd.DataFrame(extract(package.get_resource("pedigree_ancestors")))
                    ancestors = ancestors[["id","variety", "ancestor", "type", "offspring"]]
                    ancestors.to_sql("variety_ancestor", con=self._connection, index=False, if_exists="append")
                    self._logger.info("add {} ancestors".format(len(ancestors)))
                    #synonyms
                    synonyms = pd.DataFrame(extract(package.get_resource("pedigree_synonyms")))
                    synonyms = synonyms[["uid","synonym"]]
                    synonyms.to_sql("variety_synonym", con=self._connection, index=False, if_exists="append")
                    self._logger.info("add {} synonyms".format(len(synonyms)))
                    return
        self._logger.error("no pedigree found in {}".format(self._baseDir))
                    
                    
    def _exportResources(self):  
        #add new
        for filename in sorted(os.listdir(self._baseDir)):
            if os.path.isfile(os.path.join(self._baseDir,filename)):
                if filename=="pedigree.package.json":
                    pass
                elif filename.endswith(".package.json"):                    
                    self._exportResource(os.path.join(self._baseDir,filename),os.path.basename(filename)[:-13])
                    self._logger.info("processed package {}".format(filename))
                    
    def _newCollectionUid(self, collectionType, resource, internalId):
        if collectionType in self._identifiers["collections"].keys():
            if resource in self._identifiers["collections"][collectionType]:
                if str(internalId) in self._identifiers["collections"][collectionType][resource]:
                    self._identifiers["collections"][collectionType][resource][str(internalId)]["active"] = 1
                    self._identifiers["collections"][collectionType][resource][str(internalId)]["used"] = time.time()
                    return self._identifiers["collections"][collectionType][resource][str(internalId)]["uid"]
        while True:
            uid = "DC_{}".format(str(random.randint(0,999999)).zfill(6))
            if not collectionType in self._identifiers["collections"].keys():
                self._identifiers["collections"][collectionType] = {}
            if not resource in self._identifiers["collections"][collectionType].keys():
                self._identifiers["collections"][collectionType][resource] = {}            
            if not uid in self._identifiers["collections"]["uids"]:
                now = time.time()
                self._identifiers["collections"][collectionType][resource][str(internalId)] = {
                    "active": 1, "used": now, "created": now, "uid": uid
                }
                self._identifiers["collections"]["uids"][uid] = {}
                return uid             
                        
    def _newDatasetUid(self, collectionUid, internalId):
        if collectionUid in self._identifiers["collections"]["uids"].keys():
            if str(internalId) in self._identifiers["collections"]["uids"][collectionUid]:
                self._identifiers["collections"]["uids"][collectionUid][str(internalId)]["active"] = 1
                self._identifiers["collections"]["uids"][collectionUid][str(internalId)]["used"] = time.time()
                return self._identifiers["collections"]["uids"][collectionUid][str(internalId)]["uid"]
            else:
                while True:
                    uid = "DS_{}".format(str(random.randint(0,999999)).zfill(6))
                    if not uid in self._identifiers["datasets"]:
                        now = time.time()
                        self._identifiers["collections"]["uids"][collectionUid][str(internalId)] = {
                            "active": 1, "used": now, "created": now, "uid": uid
                        }
                        self._identifiers["datasets"].append(uid)
                        return uid 
        else:
            raise Exception("unknown collection uid {}".format(collectionUid))
                        
    def _createCollection(self, resource, internal_id, collection_type, resource_name, experiment_name, location):
        uid = self._newCollectionUid(collection_type, resource, internal_id)
        query = """INSERT OR IGNORE INTO `collection` (`uid`,`resource`,`internal_id`,`type`,`name`,`experiment`,`location`) 
        VALUES (?,?,?,?,?,?,?)"""
        self._connection.execute(query, (uid,resource, int(internal_id), collection_type, 
                                         resource_name, experiment_name, location))
        self._connection.commit()
        query = """SELECT `id` FROM `collection` WHERE `resource` = ? AND `internal_id` = ? AND `type` = ?"""
        self._cursor.execute(query, (resource, int(internal_id),collection_type,))  
        collection = self._cursor.fetchone()
        if collection:
            return (collection[0],uid,)
        else:
            raise Exception("could not create collection")
        
    def _createDataset(self,variety,collection_id,collection_uid,internal_id,location,type=None):
        uid = self._newDatasetUid(collection_uid, internal_id)
        query = """INSERT OR IGNORE INTO `dataset` (`uid`,`variety`,`collection_id`,`internal_id`,`location`,`type`) 
        VALUES (?,?,?,?,?,?)"""
        self._connection.execute(query, (uid,variety, int(collection_id), int(internal_id), location, type, ))
        self._connection.commit()
        
    def _getUint(maximumValue):
        if maximumValue<=np.iinfo(np.uint8).max:
            return "uint8"
        elif maximumValue<=np.iinfo(np.uint16).max:
            return "uint16"
        elif maximumValue<=np.iinfo(np.uint32).max:
            return "uint32"
        else:
            return "uint64"
    
    def _getMaxLengthColumn(df,name):
        return max(1,max([len(item) if item else 0 for item in df[name]]))
    
    def _exportResource(self, fileName:str, resourceName:str):
        package = Package(fileName)        
        if package.has_resource("varieties") and package.has_resource("metadata") and package.has_resource("experiments"):
            if (package.has_resource("sequences") or 
                    (package.has_resource("scores") 
                     and package.has_resource("markers") 
                     and package.has_resource("mappings"))):
                #get shared data
                metadata = pd.DataFrame(extract(package.get_resource("metadata")))
                metadata = {x:y for (x,y,) in zip(metadata["label"],metadata["value"])}
                varieties = pd.DataFrame(extract(package.get_resource("varieties")))
                experiments = pd.DataFrame(extract(package.get_resource("experiments")))
                if len(varieties)>0 and len(experiments)>0:
                    varieties = varieties.set_index("id")
                    experiments = experiments.set_index("id")
                    self._logger.info("store data for '{}'".format(metadata.get("name",resourceName)))
                else:
                    return
                #process sequence data
                if package.has_resource("sequences"):
                    sequences = pd.DataFrame(extract(package.get_resource("sequences")))                    
                    if len(sequences)>0:
                        for experiment_id in set(sequences["experiment_id"]):
                            collectionId=None
                            for id,row in sequences[sequences["experiment_id"]==experiment_id].iterrows():
                                variety = varieties.loc[row["variety_id"]]
                                if variety["uid"]:
                                    if collectionId==None:
                                        experiment_name = experiments.loc[experiment_id]["name"]
                                        locationResource = metadata.get("location","")
                                        locationResource = "" if locationResource==None else locationResource
                                        locationExperiment = experiments.loc[experiment_id]["location"]
                                        locationExperiment = "" if locationExperiment==None else locationExperiment
                                        location = os.path.join(locationResource,locationExperiment)
                                        (collectionId,collectionUid) = self._createCollection(
                                             resourceName, experiment_id, "kmer",
                                             metadata.get("name",resourceName), experiment_name, location)
                                    self._createDataset(variety["uid"],collectionId,collectionUid,
                                                        row["variety_id"],row["location"])

                #process marker data
                if (package.has_resource("scores") and package.has_resource("markers") and package.has_resource("mappings")):
                    scores = pd.DataFrame(extract(package.get_resource("scores")))
                    markers = pd.DataFrame(extract(package.get_resource("markers")))
                    if len(markers)>0:
                        markers = markers.set_index("id")
                    mappings = pd.DataFrame(extract(package.get_resource("mappings")))
                    #process scores
                    for score_id,row in scores.iterrows():
                        varietyList = []
                        markerList = []
                        resourceData = {"data":[],"lg":[],"marker":[],"variety":[], "varietyIndex": {}}
                        with tempfile.NamedTemporaryFile() as tmpMarkerFile:
                            sourceFilename = os.path.join(package.basepath,row["source"])
                            open_fn = gzip.open if sourceFilename.endswith(".gz") else open
                            #get all data from external file
                            with open_fn(sourceFilename, "rt") as f: 
                                reader = csv.reader(f, delimiter=",")
                                line = list(next(reader))
                                markerList = [int(item) for item in line[1:]]
                                for line in reader:
                                    items = list(line)
                                    variety_id = int(items[0])
                                    resourceData["data"].append([int(item) if item.isnumeric() else -1 for item in line[1:]])
                                    varietyList.append(variety_id)
                            if len(varietyList)>0 and len(markerList)>0:
                                experiment_id = row["experiment_id"]
                                experiment_location = experiments.loc[experiment_id]["location"]
                                experiment_location = "" if experiment_location==None else experiment_location
                                markerFilename = "markers_{}_{}_{}.h5".format(resourceName,experiment_id,score_id)
                                os.makedirs(os.path.join(self._exportDir,self.markerDirectory,
                                                         experiment_location), exist_ok=True)
                                self._logger.info("create marker database '{}'".format(markerFilename))
                                hf = h5py.File(os.path.join(self._exportDir,self.markerDirectory,
                                                            experiment_location,markerFilename), "w")
                                #data
                                hf.create_dataset("data", data=resourceData["data"], dtype=h5py.h5t.NATIVE_INT8)
                                #store linkage groups (chromosomes) and make an index to get assigned identifier later
                                lgs = [item for item in set(mappings["lg"]) if item]
                                lgIndex = {}
                                resourceData["lg"].append(("",))
                                for item in sorted(lgs):
                                    lgIndex[item]=len(resourceData["lg"])
                                    resourceData["lg"].append((item,))                                    
                                dtypeLgsList=[
                                    ("name","S{}".format(
                                        ConstructDatabase._getMaxLengthColumn(mappings,"lg")))
                                 ]
                                dtLgs=np.dtype(dtypeLgsList)
                                dsLgs=hf.create_dataset("lg", data=resourceData["lg"], dtype=dtLgs, chunks=None)
                                #mapping
                                mappingList = list(set(mappings["mapping"]))
                                dtypeMappingsList=[
                                ]
                                for mapping in mappingList:
                                    dtypeMappingsList.append(
                                        (mapping,[("lg",ConstructDatabase._getUint(len(resourceData["lg"])),),
                                                  ("pos","uint32"),]))                                    
                                dtMappings=np.dtype(dtypeMappingsList)
                                dsMappings=hf.create_dataset("mapping",(len(markerList),), 
                                                      dtype=dtMappings, chunks=None)
                                for mapping in mappingList:
                                    mappingsSubset = mappings[mappings["mapping"]==mapping]
                                    mappingsSubset = mappingsSubset[mappingsSubset["marker_id"].isin(markerList)]
                                    if len(mappingsSubset)>0:
                                        mappingsSubset = mappingsSubset.set_index("marker_id")
                                        for mappingId,mappingRow in mappingsSubset.iterrows():
                                            lgId = lgIndex.get(mappingRow["lg"],0)
                                            if not np.isnan(mappingRow["pos"]):
                                                dsMappings[markerList.index(mappingId),mapping] = (
                                                    lgId,mappingRow["pos"],)
                                #markers
                                for id,row in markers.loc[markerList].iterrows():
                                    entry = (row["name"],row["ref"],row["alt"],
                                                   row["sequence_ref"],row["sequence_alt"],
                                                   row["sequence_left"],row["sequence_right"],)
                                    entry = tuple([item if item else "" for item in entry])
                                    resourceData["marker"].append(entry)
                                dtypeMarkersList=[
                                    ("name","S{}".format(
                                        ConstructDatabase._getMaxLengthColumn(markers,"name"))),
                                    ("ref","S1"),
                                    ("alt","S1"),
                                    ("sequenceRef","S{}".format(
                                        ConstructDatabase._getMaxLengthColumn(markers,"sequence_ref"))),
                                    ("sequenceAlt","S{}".format(
                                        ConstructDatabase._getMaxLengthColumn(markers,"sequence_alt"))),
                                    ("sequenceLeft","S{}".format(
                                        ConstructDatabase._getMaxLengthColumn(markers,"sequence_left"))),
                                    ("sequenceRight","S{}".format(
                                        ConstructDatabase._getMaxLengthColumn(markers,"sequence_right")))
                                 ]
                                dtMarkers=np.dtype(dtypeMarkersList)
                                dsMarkers=hf.create_dataset("marker",data=resourceData["marker"], 
                                                      dtype=dtMarkers, chunks=None)
                                #varieties
                                for id,row in varieties.loc[varietyList].iterrows():
                                    resourceData["varietyIndex"][id] = len(resourceData["variety"])
                                    resourceData["variety"].append((row["name"] if row["name"] else "",))
                                varietyNameLength = max(1,max([len(item) for item in varieties["name"]]))
                                dtypeVarietiesList=[
                                    ("name","S{}".format(varietyNameLength))
                                 ]
                                dtVarieties=np.dtype(dtypeVarietiesList)
                                dsVarieties=hf.create_dataset("variety",data=resourceData["variety"], 
                                                      dtype=dtVarieties, chunks=None)
                                hf.close()
                                self._logger.info("register entries database '{}'".format(markerFilename))
                                #register markers in sqlite database
                                experiment_name = experiments.loc[experiment_id]["name"]
                                experiment_location = experiments.loc[experiment_id]["location"]
                                experiment_location = "" if experiment_location==None else experiment_location
                                (collectionId,collectionUid) = self._createCollection(resourceName, experiment_id, "marker",
                                                     metadata.get("name",resourceName), 
                                                     experiment_name, experiment_location)
                                for id,row in varieties.loc[varietyList].iterrows():
                                    internal_id = resourceData["varietyIndex"][id]
                                    self._createDataset(row["uid"], collectionId, collectionUid,
                                                            internal_id, markerFilename, "marker")
                                    
                            
                            
                             
                    

            
                    
                    

