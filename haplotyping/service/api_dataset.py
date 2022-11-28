from flask import Response, request, abort
from flask_restx import Namespace, Resource, abort, fields
import json, haplotyping, sqlite3

namespace = Namespace("dataset", description="All available datasets", path="/dataset")
parser = namespace.parser()

def _make_bool(text):
    return str(text).lower() in ("yes", "true", "t", "1")

def year_description(year_min,year_max):
    if year_min==None and year_max==None:
        return None
    elif year_min==None:
        return "<"+str(year_max+1)
    elif year_max==None:
        return ">"+str(year_min-1)
    elif year_min==year_max:
        return str(year_min)
    else:
        return str(year_min)+"-"+str(year_max)

def adjust_dataset_response(item):
    #variety
    if item["variety_uid"]:
        item["variety"] = {
            "uid": item["variety_uid"], 
            "name": item["variety_name"]
        }
        if item["variety_origin"]:
            item["variety"]["origin"] = {"uid": item["variety_origin"], "country": item["variety_country"]}                
        if item["variety_year_min"] or item["variety_year_max"]:
            item["variety"]["year"] = {"description": year_description(item["variety_year_min"],item["variety_year_max"]), 
                                       "min": item["variety_year_min"], 
                                       "max": item["variety_year_max"]}
    else:
        item["variety"] = None
    del item["variety_uid"]
    del item["variety_name"]
    del item["variety_origin"]
    del item["variety_year_min"]
    del item["variety_year_max"]
    #collection
    if item["collection_uid"]:
        item["collection"] = {
            "uid": item["collection_uid"], 
            "name": item["collection_name"]
        }
    else:
        item["collection"] = None
    del item["collection_uid"]
    del item["collection_name"]    
    #finished
    return item


    
@namespace.route("/")
class DatasetList(Resource):

    dataset_list = parser.copy()
    dataset_list.add_argument("start", type=int, required=False, location="values", 
                              help="paging")       
    dataset_list.add_argument("number", type=int, required=False, location="values", 
                              help="paging")
    dataset_list.add_argument("collection", type=str, required=False, location="values", 
                              help="variety has dataset from comma separated list of collection uids")
    dataset_list.add_argument("hasVariety", type=bool, required=False, location="values", 
                              help="linked to a variety")
    dataset_list.add_argument("dataType", type=str, required=False, location="values", choices=["marker","kmer","split"], 
                              help="type of dataset")
    
    dataset_set = namespace.model("uid list to get datasets", {"uids": fields.List(fields.String, attribute="items", 
                          required=True, description="list of uids")})
    
    @namespace.doc(description="Get datasets")    
    @namespace.expect(dataset_list)
    def get(self):
        try:
            start = int(request.args.get("start",0))
            number = int(request.args.get("number",10))        
            name = request.args.get("name",None)
            collection = request.args.get("collection",None)
            hasVariety = request.args.get("hasVariety",None)
            dataType = request.args.get("dataType",None)
            condition_sql = "NOT (`dataset`.`uid` IS NULL OR `dataset`.`type` IS NULL)"
            condition_variables = []
            if not collection==None:
                collection_list = collection.split(",")
                condition_sql = condition_sql + " AND (`collection`.`uid` IN ("+",".join(["?"]*len(collection_list))+"))"
                condition_variables.extend(collection_list)
            if not hasVariety==None:
                if _make_bool(hasVariety):
                    condition_sql = condition_sql + " AND NOT (`dataset`.`variety` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`dataset`.`variety` IS NULL)"       
            if not dataType==None:
                if dataType=="kmer":
                    condition_sql = condition_sql + " AND ((`dataset`.`type` IS 'kmer') OR (`dataset`.`type` IS 'split'))"  
                elif dataType=="split":
                    condition_sql = condition_sql + " AND (`dataset`.`type` IS 'split')" 
                elif dataType=="marker":
                    condition_sql = condition_sql + " AND (`dataset`.`type` IS 'marker')"
                else:
                    abort(422, "incorrect dataType condition "+str(dataType))   
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT COUNT(DISTINCT `dataset`.`id`) AS `number` \
                            FROM `dataset` \
                            LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
                            LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            WHERE "+condition_sql, tuple(condition_variables))  
            total = cursor.fetchone()[0]
            if start<total:
                cursor.execute("SELECT `dataset`.`uid`, \
                                `dataset`.`type`, \
                                `collection`.`uid` AS `collection_uid`, \
                                `collection`.`name` AS `collection_name`, \
                                `variety`.`uid` AS `variety_uid`, \
                                `variety`.`name` AS `variety_name`, \
                                `variety`.`origin` AS `variety_origin`, \
                                `country`.`name` AS `variety_country`, \
                                `variety`.`year_min` AS `variety_year_min`, \
                                `variety`.`year_max` AS `variety_year_max` \
                                FROM `dataset` \
                                LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                                WHERE "+condition_sql+" \
                                GROUP BY `dataset`.`id` \
                                ORDER BY `dataset`.`uid` \
                                LIMIT ?,?",tuple(condition_variables + [start,number]))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                resultList = []
            for i in range(len(resultList)):
                resultList[i] = adjust_dataset_response(resultList[i])                
            response = {"start": start, "number": number, "total": total, "list": resultList}
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
    @namespace.doc(description="Get multiple datasets by uid")    
    @namespace.expect(dataset_set)
    def post(self):
        uids = namespace.payload.get("uids",[])
        try:
            if len(uids)>0:
                condition_sql = "`dataset`.`uid` IN ("+(",".join(["?"]*len(uids)))+")"
                condition_variables = [] + uids;
                db_connection = haplotyping.service.API.get_db_connection()
                db_connection.row_factory = sqlite3.Row
                cursor = db_connection.cursor()
                cursor.execute("SELECT COUNT(DISTINCT `dataset`.`id`) AS `number` \
                                FROM `dataset` \
                                LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                                WHERE "+condition_sql, tuple(condition_variables))  
                total = cursor.fetchone()[0]
                cursor.execute("SELECT `dataset`.`uid`, \
                                `dataset`.`type`, \
                                `collection`.`uid` AS `collection_uid`, \
                                `collection`.`name` AS `collection_name`, \
                                `variety`.`uid` AS `variety_uid`, \
                                `variety`.`name` AS `variety_name`, \
                                `variety`.`origin` AS `variety_origin`, \
                                `country`.`name` AS `variety_country`, \
                                `variety`.`year_min` AS `variety_year_min`, \
                                `variety`.`year_max` AS `variety_year_max` \
                                FROM `dataset` \
                                LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                                WHERE "+condition_sql+" \
                                GROUP BY `dataset`.`id` \
                                ORDER BY `dataset`.`uid`",tuple(condition_variables))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                total = 0
                resultList = []
            for i in range(len(resultList)):
                resultList[i] = adjust_dataset_response(resultList[i])                
            response = {"total": total, "list": resultList}
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>")
@namespace.doc(params={"uid": "unique identifier dataset"})
class DatasetId(Resource):
    
    @namespace.doc(description="Get dataset by uid")
    def get(self,uid):
        try:
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT `dataset`.`uid`, \
                            `dataset`.`type`, \
                            `collection`.`uid` AS `collection_uid`, \
                            `collection`.`name` AS `collection_name`, \
                            `variety`.`uid` AS `variety_uid`, \
                            `variety`.`name` AS `variety_name`, \
                            `variety`.`origin` AS `variety_origin`, \
                            `country`.`name` AS `variety_country`, \
                            `variety`.`year_min` AS `variety_year_min`, \
                            `variety`.`year_max` AS `variety_year_max` \
                            FROM `dataset` \
                            LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
                            LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            WHERE `dataset`.`uid` = ? \
                            GROUP BY `dataset`.`id`",(uid,))  
            data = cursor.fetchone()
            if data:
                response = adjust_dataset_response(dict(data))
                return Response(json.dumps(response), mimetype="application/json")
            else:
                abort(404, "no dataset with uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))