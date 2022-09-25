#!python

import os,logging
import sqlite3
from sqlite3 import Error
import h5py

class CheckDatabase:
    
    dbFilename = "db.sqlite"
    markerDirectory = "marker"
    kmerDirectory = "kmer"
    
    def __init__(self, dbdir:str, markerdir:str = None, kmerdir:str = None):
        #data
        self._dbDir = dbdir
        self._markerDir = os.path.join(self._dbDir,self.markerDirectory) if markerdir==None else markerdir
        self._kmerDir = os.path.join(self._dbDir,self.kmerDirectory) if kmerdir==None else kmerdir
        #logging
        self._logger = logging.getLogger(__name__)
        #check directories
        if not os.path.exists(self._markerDir):
            self._logger.error("marker directory {} not found".format(self._markerDir))
            return
        if not os.path.exists(self._kmerDir):
            self._logger.error("kmer directory {} not found".format(self._kmerDir))
            return
        #connection
        self._connection = None
        self._cursor = None
        try:
            self._connect_sqlite3()
            self._reset_dataset_types()
            self._check_markers()
            self._check_kmers()
        except Error as e:
            print(e)
        finally:
            self._close_sqlite3()

    def _connect_sqlite3(self):
        if os.path.exists(os.path.join(self._dbDir,self.dbFilename)):
            self._logger.info("open connection {}".format(self.dbFilename))
            self._connection = sqlite3.connect(os.path.join(self._dbDir,self.dbFilename))  
            self._connection.row_factory = sqlite3.Row
            self._cursor = self._connection.cursor()
        else:
            raise Exception("no {} found".format(self.dbFilename))
        
    def _close_sqlite3(self):
        if self._cursor:
            self._cursor.close()
        if self._connection:            
            self._connection.close()
            
    def _reset_dataset_types(self):
        self._logger.debug("reset dataset types")
        query = """UPDATE `dataset` SET `type` = NULL"""
        self._connection.execute(query)
        self._connection.commit()
            
    def _check_markers(self):
        query = """SELECT * FROM `collection` WHERE `type` = 'marker' ORDER BY `id`"""
        self._cursor.execute(query)
        collections = self._cursor.fetchall()        
        for row in collections:
            collection = dict(zip(row.keys(), row))
            if collection["location"]!=None:
                if not os.path.exists(os.path.join(self._markerDir,collection["location"])):
                    self._logger.error("marker collection {}: location {} not found".format(
                        collection["name"], collection["location"]))
                    continue
                else:
                    self._logger.info("marker collection {}: check datasets".format(collection["name"]))
            #get datasets
            query = """SELECT * FROM `dataset` WHERE `collection_id` = ? ORDER BY `id`"""
            self._cursor.execute(query, (collection["id"],))  
            datasets = self._cursor.fetchall()
            markerDbs = {}
            for row in datasets:
                dataset = dict(zip(row.keys(), row))
                if dataset["location"]!=None:
                    if not dataset["location"] in markerDbs.keys():
                        if not collection["location"]==None:
                            dbLocation = os.path.join(self._markerDir,collection["location"],dataset["location"])
                        else:
                            dbLocation = os.path.join(self._markerDir,dataset["location"])
                        if not os.path.exists(dbLocation):   
                            self._logger.error("marker collection {}: database {} not found".format(
                                collection["name"], dbLocation))
                            markerDbs[dataset["location"]] = []
                        else:
                            with h5py.File(dbLocation, "r") as hf:
                                markerDbs[dataset["location"]] = list(range(len(hf["variety"])))
                    if dataset["internal_id"] in markerDbs[dataset["location"]]:
                        query = """UPDATE `dataset` SET `type` = 'marker' WHERE `id` = ? AND `collection_id` = ?"""
                        self._cursor.execute(query, (dataset["id"], collection["id"],))  
                    else:
                        self._logger.error("marker collection {}: markers for {} not found in database {}".format(
                            collection["name"], dataset["variety"], dataset["location"]))
            self._connection.commit()
            
    def _check_kmers(self):
        query = """SELECT * FROM `collection` WHERE `type` = 'kmer' ORDER BY `id`"""
        self._cursor.execute(query)
        collections = self._cursor.fetchall()        
        for row in collections:
            collection = dict(zip(row.keys(), row))
            if collection["location"]!=None:
                if not os.path.exists(os.path.join(self._kmerDir,collection["location"])):
                    self._logger.error("kmer collection {}: location {} not found".format(
                        collection["name"], collection["location"]))
                    continue
                else:
                    self._logger.info("kmer collection {}: check datasets".format(collection["name"]))
            #get datasets
            query = """SELECT * FROM `dataset` WHERE `collection_id` = ? ORDER BY `id`"""
            self._cursor.execute(query, (collection["id"],))  
            datasets = self._cursor.fetchall()
            locationsNotFound = 0
            kmerNotFound = 0
            splitNotFound = 0
            for row in datasets:
                dataset = dict(zip(row.keys(), row))
                if dataset["location"]!=None:
                    if not collection["location"]==None:
                        kmerLocation = os.path.join(self._kmerDir,collection["location"],dataset["location"])
                    else:
                        kmerLocation = os.path.join(self._kmerDir,dataset["location"])
                    if not os.path.exists(kmerLocation):   
                        self._logger.debug("kmer collection {}: location {} not found".format(
                                collection["name"], kmerLocation))
                        locationsNotFound+=1
                    else:
                        if not os.path.exists(os.path.join(kmerLocation,"kmer.kmc.kmc_pre")):  
                            kmerNotFound+=1
                        elif not os.path.exists(os.path.join(kmerLocation,"kmer.kmc.kmc_suf")):  
                            kmerNotFound+=1
                        else:
                            datasetType = "kmer"
                            if not os.path.exists(os.path.join(kmerLocation,"kmer.data.h5")): 
                                splitNotFound+=1
                            else:
                                datasetType = "split"
                            query = """UPDATE `dataset` SET `type` = ? WHERE `id` = ? AND `collection_id` = ?"""
                            self._cursor.execute(query, (datasetType, dataset["id"], collection["id"],))  
            self._connection.commit()            
            if locationsNotFound>0:
                self._logger.error("kmer collection {}: {} location(s) not found".format(
                                collection["name"], locationsNotFound))
            if kmerNotFound>0:
                self._logger.error("kmer collection {}: {} k-mer databases not found".format(
                                collection["name"], locationsNotFound))
            if splitNotFound>0:
                self._logger.error("kmer collection {}: {} splitting k-mer databases not found".format(
                                collection["name"], locationsNotFound))
                            
                    

                
            
            
    

