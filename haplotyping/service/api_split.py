from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
import json, haplotyping, sqlite3

from haplotyping.service.split import Split

namespace = Namespace("split", description="Splitting k-mer information for a dataset (WORK IN PROGRESS)", path="/split")
parser = namespace.parser()
            
@namespace.route("/<uid>/kmer/<kmer>")
class SplitKmer(Resource):
    
    @namespace.doc(description="Get splitting k-mer information from dataset defined by uid for specified k-mer")
    @namespace.doc(params={"uid": "unique identifier dataset","split": "splitting k-mer"})
    def get(self,uid,split):
        try:
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT `dataset`.`uid`, \
                            `dataset`.`location_split`, \
                            `collection`.`location` AS `location` \
                            FROM `dataset` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            WHERE `dataset`.`uid` = ? \
                            AND NOT `dataset`.`location_split` IS NULL \
                            AND NOT `collection`.`location` IS NULL",(str(uid),))  
            data = cursor.fetchone()
            if data:
                location_split = haplotyping.service.API.get_data_location() + data["location"] + data["location_split"]
                response = {}
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/base/<base>")
class SplitBase(Resource):
    
    @namespace.doc(description="Get splitting k-mer base information from dataset defined by uid for specified base")
    @namespace.doc(params={"uid": "unique identifier dataset","base": "splitting base"})
    def get(self,uid,base):
        try:
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT `dataset`.`uid`, \
                            `dataset`.`location_split`, \
                            `collection`.`location` AS `location` \
                            FROM `dataset` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            WHERE `dataset`.`uid` = ? \
                            AND NOT `dataset`.`location_split` IS NULL \
                            AND NOT `collection`.`location` IS NULL",(str(uid),))  
            data = cursor.fetchone()
            if data:
                location_split = haplotyping.service.API.get_data_location() + data["location"] + data["location_split"]
                response = {}
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/connected/<kmer>")
class SplitConnected(Resource):
    
    @namespace.doc(description="Get connected k-mer information from dataset defined by uid for specified k-mer")
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "splitting k-mer"})
    def get(self,uid,kmer):
        try:
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT `dataset`.`uid`, \
                            `dataset`.`location_split`, \
                            `collection`.`location` AS `location` \
                            FROM `dataset` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            WHERE `dataset`.`uid` = ? \
                            AND NOT `dataset`.`location_split` IS NULL \
                            AND NOT `collection`.`location` IS NULL",(str(uid),))  
            data = cursor.fetchone()
            if data:
                location_split = haplotyping.service.API.get_data_location() + data["location"] + data["location_split"]
                response = {}
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
