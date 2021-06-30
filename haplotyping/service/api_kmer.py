from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
import json, haplotyping, sqlite3, subprocess

from haplotyping.service.kmer_kmc import Kmer as KmerKMC

namespace = Namespace("kmer", description="K-mer frequencies for a dataset", path="/kmer")
parser = namespace.parser()
            
@namespace.route("/<uid>/<kmer>")
class KmerSingle(Resource):
    
    dataset_kmers = parser.copy()
    dataset_kmers.add_argument("mismatches", type=int, required=False, location="values", 
                              help="maximum number of mismatches", choices=[0,1,2])
    
    @namespace.doc(description="Get k-mer frequencies from dataset defined by uid")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "k-mer"})
    def get(self,uid,kmer):
        mm = min(2,max(0,int(request.args.get("mismatches",0))))
        try:
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
            if data:
                location_kmc = data["location"] + data["location_kmer"]
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary = haplotyping.service.API.get_kmc_query_binary()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,[kmer],mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmer))
                elif kmc_query_binary:
                    response = KmerKMC.kmc_binary(kmc_query_binary,location_kmc,[kmer],mm)
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmer))
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>")
class KmerMultiple(Resource):            
                
    dataset_kmers = namespace.model("k-mer list to get frequencies from dataset", {
        "kmers": fields.List(fields.String, attribute="items", required=True, description="list of k-mers"),
        "mismatches": fields.Integer(min=0,max=2, required=False, description="maximum number of mismatches")
    })
                
    @namespace.doc(description="Get k-mer frequencies from dataset defined by uid")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset"})
    def post(self,uid):
        mm = min(2,max(0,namespace.payload.get("mismatches",0)))
        kmers = namespace.payload.get("kmers",[])
        try:
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
            if data:
                location_kmc = data["location"] + data["location_kmer"]
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary = haplotyping.service.API.get_kmc_query_binary()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmers))
                elif kmc_query_binary:
                    response = KmerKMC.kmc_binary(kmc_query_binary,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmers))
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
        
            
            