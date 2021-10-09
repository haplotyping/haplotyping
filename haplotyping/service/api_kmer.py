from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
from flask_caching import Cache
import json, haplotyping, sqlite3

from haplotyping.service.kmer_kmc import Kmer as KmerKMC

def _make_cache_key(*args, **kwargs):
    cacheKey = "%s_%s_%s_%s" % (args[0].__class__.__name__, str(request.path), str(request.args), str(namespace.payload))
    return cacheKey
    
namespace = Namespace("kmer", description="K-mer frequencies for a dataset", path="/kmer")
cache = Cache()
parser = namespace.parser()
            
def _getDataset(uid):
    db_connection = haplotyping.service.API.get_db_connection()
    db_connection.row_factory = sqlite3.Row
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                    `dataset`.`location_kmer`, \
                    `collection`.`location` AS `location` \
                    FROM `dataset` \
                    LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                    WHERE `dataset`.`uid` = ? \
                    AND NOT `dataset`.`location_kmer` IS NULL \
                    AND NOT `collection`.`location` IS NULL",(str(uid),))  
    data = cursor.fetchone()
    return data

@namespace.route("/<uid>/<kmer>")
class KmerSingle(Resource):
    
    dataset_kmers = parser.copy()
    dataset_kmers.add_argument("mismatches", type=int, required=False, location="values", 
                              help="maximum number of mismatches", choices=[0,1,2])
    
    @namespace.doc(description="Get k-mer frequencies from dataset defined by uid")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "k-mer"})
    @cache.cached(key_prefix="kmer_", make_cache_key=_make_cache_key)
    def get(self,uid,kmer):
        mm = min(2,max(0,int(request.args.get("mismatches",0))))
        try:
            data = _getDataset(uid)
            if data:
                location_kmc = haplotyping.service.API.get_data_location() + data["location"] + data["location_kmer"]
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,[kmer],mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmer))
                elif kmc_query_binary_location:
                    response = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location,location_kmc,[kmer],mm)
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmer))
                    response["info"] = KmerKMC.kmc_binary_info(kmc_query_binary_location,location_kmc)
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>")
class KmerMultiple(Resource):            
                
    dataset_kmers = namespace.model("k-mer list to get frequencies", {
        "kmers": fields.List(fields.String, attribute="items", required=True, description="list of k-mers"),
        "mismatches": fields.Integer(min=0,max=2, required=False, description="maximum number of mismatches")
    })
                
    @namespace.doc(description="Get k-mer frequencies from dataset defined by uid")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @cache.cached(key_prefix="kmer_", make_cache_key=_make_cache_key)
    def post(self,uid):
        mm = min(2,max(0,namespace.payload.get("mismatches",0)))
        kmers = namespace.payload.get("kmers",[])
        try:
            data = _getDataset(uid)
            if data:
                location_kmc = haplotyping.service.API.get_data_location() + data["location"] + data["location_kmer"]
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmers))
                elif kmc_query_binary_location:
                    response = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmers))
                    response["info"] = KmerKMC.kmc_binary_info(kmc_query_binary_location,location_kmc)
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
        
            
@namespace.route("/<uid>/sequence")
class KmerSequence(Resource):            
                
    sequence_data = namespace.model("sequence to get k-mer frequencies", {"sequence": 
                                            fields.String(required=True, description="sequence"),
                                    "mismatches": fields.Integer(min=0,max=2, required=False, 
                                                                 description="maximum number of mismatches")})
    
    @namespace.doc(description="Get k-mer frequencies for a sequence from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(sequence_data)
    @cache.cached(key_prefix="kmer_", make_cache_key=_make_cache_key)
    def post(self,uid):  
        mm = min(2,max(0,namespace.payload.get("mismatches",0)))
        sequence = namespace.payload.get("sequence","")
        kmers = []
        try:
            data = _getDataset(uid)
            if data:
                location_kmc = haplotyping.service.API.get_data_location() + data["location"] + data["location_kmer"]
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,[],0)
                    if not response:
                        abort(500,"no response using kmc query library")
                    k = response["info"].get("kmer_length",0)
                    kmers = [sequence[i:i+k] for i in range(len(sequence)-(k-1)) if k>0]
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmers))
                elif kmc_query_binary_location:
                    info = KmerKMC.kmc_binary_info(kmc_query_binary_location,location_kmc)
                    k = info.get("kmer_length",0)
                    kmers = [sequence[i:i+k] for i in range(len(sequence)-(k-1)) if k>0]
                    response = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location,location_kmc,kmers,mm)
                    response["info"] = info;
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmers))
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))