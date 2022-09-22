from flask import Response, request, abort
from flask_restx import Namespace, Resource
import json, haplotyping, sqlite3

namespace = Namespace("country", description="Varieties originate from these countries", path="/country")

@namespace.route("/")
class CountryList(Resource):
    
    @namespace.doc(description="Get countries")
    @namespace.param("start", "start of listed results", type="int", required=False)
    @namespace.param("number", "(maximum) number of listed results", type="int", required=False)
    def get(self):
        try:
            start = int(request.args.get("start",0))
            number = int(request.args.get("number",10))        
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT COUNT(*) AS `number` FROM `country`")  
            total = cursor.fetchone()[0]
            if start<total:
                cursor.execute("SELECT `country`.`uid`, `country`.`name` FROM `country` \
                                INNER JOIN `variety` ON `country`.`uid` = `variety`.`origin` \
                                GROUP BY `country`.`id` ORDER BY `country`.`uid` LIMIT ?,? ",(start,number,))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                resultList = []
            response = {"start": start, "number": number, "total": total, "list": resultList}
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))

