from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
from flask_caching import Cache
import json, haplotyping, sqlite3, os

from haplotyping.service.split import Split

def _make_cache_key(*args, **kwargs):
    #todo: find better solution
    try:
        cacheKey = "%s_%s_%s_%s" % (args[0].__class__.__name__, str(request.path), str(request.args), str(namespace.payload))
    except:
        cacheKey = "%s_%s_%s" % (args[0].__class__.__name__, str(request.path), str(request.args))
    return cacheKey
   
namespace = Namespace("split", description="Splitting k-mer information for a dataset", path="/split")
cache = Cache()
parser = namespace.parser()

def _getDataset(uid):
    db_connection = haplotyping.service.API.get_db_connection()
    db_connection.row_factory = sqlite3.Row
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                    `dataset`.`type`, \
                    `dataset`.`location` AS `dataset_location`, \
                    `collection`.`location` AS `collection_location` \
                    FROM `dataset` \
                    LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                    WHERE `dataset`.`uid` = ? \
                    AND (`dataset`.`type` = 'kmer' OR `dataset`.`type` = 'split') \
                    AND NOT `collection`.`location` IS NULL",(str(uid),))  
    data = cursor.fetchone()
    return data
            
@namespace.route("/<uid>/info")
class SplitKmerInfo(Resource):
    
    @namespace.doc(description="Get information describing dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid):
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.info(location_split)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/kmer/<kmer>")
class SplitKmerSingle(Resource):
    
    @namespace.doc(description="Get splitting k-mer information from dataset defined by uid for specified k-mer")
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "splitting k-mer"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,kmer):
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_info(location_split, kmer)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/kmer")
class SplitKmerMultiple(Resource):
    
    dataset_kmers = namespace.model("k-mer list to get splitting k-mer information", {
        "kmers": fields.List(fields.String, attribute="items", required=True, description="list of k-mers")
    })
    
    @namespace.doc(description="Get splitting k-mer information for a list of k-mers from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(dataset_kmers)
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):
        kmers = namespace.payload.get("kmers",[])
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_list_info(location_split, kmers)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/kmer/direct/<kmer>")
class SplitKmerDirectSingle(Resource):
    
    @namespace.doc(description="Get direct connections splitting k-mer from dataset defined by uid for specified k-mer")
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "splitting k-mer"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,kmer):
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_direct(location_split, kmer)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))    
            
@namespace.route("/<uid>/kmer/direct")
class SplitKmerDirectMultiple(Resource):
    
    dataset_kmers = namespace.model("k-mer list to get direct connections splitting k-mer", {
        "kmers": fields.List(fields.String, attribute="items", required=True, description="list of k-mers")
    })
    
    @namespace.doc(description="Get direct connections splitting k-mer for a list of k-mers from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(dataset_kmers)
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):
        kmers = namespace.payload.get("kmers",[])
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_list_direct(location_split, kmers)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))  
            
@namespace.route("/<uid>/kmer/connected/<kmer>")
class SplitKmerConnectedSingle(Resource):
    
    @namespace.doc(description="Get connected splitting k-mers from dataset defined by uid containing specified k-mer")
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "splitting k-mer"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,kmer):
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_connected(location_split, kmer)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))      
            
@namespace.route("/<uid>/kmer/connected")
class SplitKmerConnectedMultiple(Resource):
    
    dataset_kmers = namespace.model("k-mer list to get connected splitting k-mers", {
        "kmers": fields.List(fields.String, attribute="items", required=True, description="list of k-mers")
    })
    
    @namespace.doc(description="Get connected splitting k-mers for a list of k-mers from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(dataset_kmers)
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):
        kmers = namespace.payload.get("kmers",[])
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_list_connected(location_split, kmers)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))              
            
@namespace.route("/<uid>/kmer/sequence")
class SplitKmerSequence(Resource):
    
    sequence_data = namespace.model("sequence to get splitting k-mer information", {"sequence": 
                                            fields.String(required=True, description="sequence")})
    
    @namespace.doc(description="Get splitting k-mer information for a sequence from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(sequence_data)
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):
        sequence = namespace.payload.get("sequence","")
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.kmer_sequence_info(location_split, sequence)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))      
            
@namespace.route("/<uid>/base/<base>")
class SplitBaseSingle(Resource):
    
    @namespace.doc(description="Get splitting k-mer base information from dataset defined by uid for specified k-mer base")
    @namespace.doc(params={"uid": "unique identifier dataset","base": "splitting k-mer base"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,base):
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.base_info(location_split, base)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/base")
class SplitBaseMultiple(Resource):
    
    dataset_bases = namespace.model("k-mer base list to get splitting k-mer base information", {
        "bases": fields.List(fields.String, attribute="items", required=True, description="list of bases")
    })
    
    @namespace.doc(description="Get splitting k-mer base information for a list of k-mer bases from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(dataset_bases)
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):
        bases = namespace.payload.get("bases",[])
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.data.h5")
                else:
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.data.h5")
                response = Split.base_list_info(location_split, bases)
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with splitting k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))   
            
