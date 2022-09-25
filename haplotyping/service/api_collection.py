from flask import Response, request, abort
from flask_restx import Namespace, Resource
import json, haplotyping, sqlite3

namespace = Namespace("collection", description="Datasets all belong to a collection", path="/collection")

def adjust_collection_response(item):
    item["datasets"] = {
        "total": item["datasets_total"], 
        "kmer": item["datasets_kmer"],
        "split": item["datasets_split"],
        "marker": item["datasets_marker"]
    }
    del item["datasets_total"]
    del item["datasets_kmer"]
    del item["datasets_split"]
    del item["datasets_marker"]
    #finished
    return item


@namespace.route("/")
class CollectionList(Resource):
    
    @namespace.doc(description="Get collections")
    @namespace.param("start", "start of listed results", type="int", required=False)
    @namespace.param("number", "(maximum) number of listed results", type="int", required=False)
    def get(self):
        try:
            start = int(request.args.get("start",0))
            number = int(request.args.get("number",10))        
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT COUNT(*) AS `number` FROM `collection`")  
            total = cursor.fetchone()[0]
            if start<total:
                cursor.execute("SELECT `collection`.`uid`, `collection`.`name`,\
                                       `collection`.`experiment`, `collection`.`type`,\
                                       COUNT(DISTINCT(`dataset`.`id`)) AS `datasets`\
                                FROM `collection`\
                                INNER JOIN `dataset` ON `collection`.`id` = `dataset`.`collection_id`\
                                AND NOT `dataset`.`uid` IS NULL\
                                AND NOT `dataset`.`type` IS NULL\
                                GROUP BY `collection`.`id` ORDER BY `name`,\
                                `collection`.`type`, `collection`.`id` LIMIT ?,?",(start,number,))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                resultList = []
            response = {"start": start, "number": number, "total": total, "list": resultList}
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))

