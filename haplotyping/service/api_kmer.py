from flask import Response, request
from flask_restx import Namespace, Resource, abort, fields
import json, haplotyping, sqlite3

namespace = Namespace("kmer", description="K-mer frequencies for a dataset", path="/kmer")
parser = namespace.parser()

            
@namespace.route("/<uid>/<kmer>")
@namespace.doc(params={"uid": "unique identifier dataset","kmer": "k-mer"})
class Kmer(Resource):
    
    @namespace.doc(description="Get k-mer frequency from dataset defined by uid")
    def get(self,uid,kmer):
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
                            AND NOT `collection`.`location` IS NULL",(uid,))  
            data = cursor.fetchone()
            if data:
                response = adjust_dataset_response(dict(data), db_connection)
                return Response(json.dumps(response), mimetype="application/json")
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))