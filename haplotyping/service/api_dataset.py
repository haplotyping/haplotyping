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
    if item["variety_uid"]:
        item["variety"] = {
            "uid": item["variety_uid"], 
            "name": item["variety_name"],
            "origin": item["variety_origin"],
            "year": year_description(item["variety_year_min"],item["variety_year_max"])
        }
    else:
        item["variety"] = None
    del item["variety_uid"]
    del item["variety_name"]
    del item["variety_origin"]
    del item["variety_year_min"]
    del item["variety_year_max"]
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
                              help="variety has dataset from comma separated list of collections")
    dataset_list.add_argument("variety", type=bool, required=False, location="values", 
                              help="linked to a variety")
    dataset_list.add_argument("kmer", type=bool, required=False, location="values", 
                              help="k-mer information available")
    dataset_list.add_argument("split", type=bool, required=False, location="values", 
                              help="splitting k-mer information available")
    dataset_list.add_argument("marker", type=bool, required=False, location="values", 
                              help="marker information available")
    
    @namespace.doc(description="Get datasets")    
    @namespace.expect(dataset_list)
    def get(self):
        try:
            start = int(request.args.get("start",0))
            number = int(request.args.get("number",10))        
            name = request.args.get("name",None)
            collection = request.args.get("collection",None)
            variety = request.args.get("variety",None)
            kmer = request.args.get("kmer",None)
            split = request.args.get("split",None)
            marker = request.args.get("marker",None)
            condition_sql = "1"
            condition_variables = []
            if not collection==None:
                collection_list = collection.split(",")
                condition_sql = condition_sql + " AND (`collection`.`name` IN ("+",".join(["?"]*len(collection_list))+"))"
                condition_variables.extend(collection_list)
            if not variety==None:
                if _make_bool(variety):
                    condition_sql = condition_sql + " AND NOT (`dataset`.`variety` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`dataset`.`variety` IS NULL)"       
            if not kmer==None:
                if _make_bool(kmer):
                    condition_sql = condition_sql + " AND NOT (`dataset`.`location_kmer` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`dataset`.`location_kmer` IS NULL)"
            if not split==None:
                if _make_bool(split):
                    condition_sql = condition_sql + " AND NOT (`dataset`.`location_split` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`dataset`.`location_split` IS NULL)"   
            if not marker==None:
                if _make_bool(marker):
                    condition_sql = condition_sql + " AND NOT (`dataset`.`location_marker` IS NULL \
                                                            OR `dataset`.`marker_id` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`dataset`.`location_marker` IS NULL \
                                                            OR `dataset`.`marker_id` IS NULL)"   
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT COUNT(DISTINCT `dataset`.`id`) AS `number` \
                            FROM `dataset` \
                            LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            WHERE "+condition_sql, tuple(condition_variables))  
            total = cursor.fetchone()[0]
            if start<total:
                cursor.execute("SELECT `dataset`.`uid`, \
                                (CASE WHEN `dataset`.`location_kmer` IS NULL THEN 0 ELSE 1 END) AS `kmer`, \
                                (CASE WHEN `dataset`.`location_split` IS NULL THEN 0 ELSE 1 END) AS `split`, \
                                (CASE WHEN `dataset`.`location_marker` IS NULL OR \
                                    `dataset`.`marker_id` IS NULL THEN 0 ELSE 1 END) AS `marker`, \
                                `collection`.`name` AS `collection`, \
                                `variety`.`uid` AS `variety_uid`, \
                                `variety`.`name` AS `variety_name`, \
                                `variety`.`origin` AS `variety_origin`, \
                                `variety`.`year_min` AS `variety_year_min`, \
                                `variety`.`year_max` AS `variety_year_max` \
                                FROM `dataset` \
                                LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
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
                            (CASE WHEN `dataset`.`location_kmer` IS NULL THEN 0 ELSE 1 END) AS `kmer`, \
                            (CASE WHEN `dataset`.`location_split` IS NULL THEN 0 ELSE 1 END) AS `split`, \
                            (CASE WHEN `dataset`.`location_marker` IS NULL OR \
                                    `dataset`.`marker_id` IS NULL THEN 0 ELSE 1 END) AS `marker`, \
                            `collection`.`name` AS `collection`, \
                            `variety`.`uid` AS `variety_uid`, \
                            `variety`.`name` AS `variety_name`, \
                            `variety`.`origin` AS `variety_origin`, \
                            `variety`.`year_min` AS `variety_year_min`, \
                            `variety`.`year_max` AS `variety_year_max` \
                            FROM `dataset` \
                            LEFT JOIN `variety` ON `dataset`.`variety` = `variety`.`uid` \
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