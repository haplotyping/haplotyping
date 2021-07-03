from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
import json, haplotyping, sqlite3, re

namespace = Namespace("variety", description="All available varieties", path="/variety")
parser = namespace.parser()

yearPattern = re.compile(r"^[\d]{4}$")
yearMaxPattern = re.compile(r"^\<[\d]{4}$")
yearMinPattern = re.compile(r"^\>[\d]{4}$")
yearRangePattern = re.compile(r"^[\d]{4}\-[\d]{4}$")

def _make_bool(text):
    return str(text).lower() in ("yes", "true", "t", "1")

def get_variety_datasets(uid, db_connection):
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                        (CASE WHEN `dataset`.`location_kmer` IS NULL THEN 0 ELSE 1 END) AS `kmer`, \
                        (CASE WHEN `dataset`.`location_split` IS NULL THEN 0 ELSE 1 END) AS `split`, \
                        `collection`.`name` AS `collection` \
                        FROM `dataset` \
                        LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                        WHERE `dataset`.`variety` = ? \
                        ORDER BY `collection`.`name`, `dataset`.`uid`",(uid,)) 
    return [dict(row) for row in cursor.fetchall()]  

def parents_tree(id,parentData):
    parents = []   
    for item in parentData:
        if item["offspring"]==id:
            if item["ancestor"]==None:
                parent = {"type": item["type"],
                          "parents": parents_tree(item["id"],parentData)}
                parents.append(parent)
            else:
                parent = {"uid": item["ancestor"],
                          "name": item["name"],
                          "relation": item["type"]}
                parents.append(parent)
    return parents

def offspring_list(offspringData):
    offspring = []
    for item in offspringData:
        offspring_item = {
            "uid": item["variety"],
            "name": item["name"],
            "relation": "unknown" if item["type"]==None else item["type"],
            "direct": True if item["offspring"]==None else False
        }
        offspring.append(offspring_item)
    return offspring
        
def get_variety_parents(uid, db_connection):
    cursor = db_connection.cursor()
    cursor.execute("SELECT `variety_ancestor`.*, `variety`.`name` FROM `variety_ancestor` \
                    LEFT JOIN `variety` ON `variety_ancestor`.`ancestor` = `variety`.`uid` \
                    WHERE `variety_ancestor`.`variety` = ? \
                    ORDER BY `variety_ancestor`.`type` DESC",(uid,)) 
    parentData = [dict(row) for row in cursor.fetchall()]
    return parents_tree(None,parentData)

def get_variety_offspring(uid, db_connection):
    cursor = db_connection.cursor()
    cursor.execute("SELECT `variety_ancestor`.*, `variety`.`name` FROM `variety_ancestor` \
                    LEFT JOIN `variety` ON `variety_ancestor`.`variety` = `variety`.`uid` \
                    WHERE `variety_ancestor`.`ancestor` = ? \
                    ORDER BY `variety`.`name`, `variety`.`uid`",(uid,)) 
    offspringData = [dict(row) for row in cursor.fetchall()]
    return offspring_list(offspringData)

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
    
def adjust_variety_response(item, db_connection):
    #origin
    item["origin"] = {"uid": item["origin"], "country": item["country"], 
                        "intermediate-region": item["intermediate-region"],
                        "sub-region": item["sub-region"],
                        "region": item["region"]
                        }
    del item["country"]
    del item["region"]
    del item["sub-region"]
    del item["intermediate-region"]
    #year
    item["year"] = {"description": year_description(item["year_min"],item["year_max"]), 
                    "min": item["year_min"], 
                    "max": item["year_max"]}            
    del item["year_min"]
    del item["year_max"]
    #datasets
    item["datasets"] = get_variety_datasets(item["uid"], db_connection)
    #parents
    item["parents"] = get_variety_parents(item["uid"], db_connection)
    #offspring
    item["offspring"] = get_variety_offspring(item["uid"], db_connection)
    #finished
    return item
    

variety_list = parser.copy()
variety_list.add_argument("start", type=int, required=False, location="args", 
                          help="paging")       
variety_list.add_argument("number", type=int, required=False, location="args", 
                          help="paging")
variety_list.add_argument("name", type=str, required=False, location="args", 
                          help="name of variety")
variety_list.add_argument("origin", type=str, required=False, location="args", 
                          help="comma separated list of country codes")
variety_list.add_argument("year", type=str, required=False, location="args", 
                          help="year of variety (e.g. '1995', '<1995', '>1995', '1990-1995')")
variety_list.add_argument("collection", type=str, required=False, location="args", 
                          help="variety has dataset from comma separated list of collections")
variety_list.add_argument("dataset", type=str, required=False, location="args", 
                          help="variety has dataset (of specific type)", choices=["any","none","kmer","split"])
variety_list.add_argument("parents", type=bool, required=False, location="args", 
                          help="parent(s) known for this variety")
variety_list.add_argument("offspring", type=bool, required=False, location="args", 
                          help="offspring known for this variety")
    
@namespace.route("/")
class VarietyList(Resource):
    
    @namespace.doc(description="Get varieties")    
    @namespace.expect(variety_list)
    def get(self):
        try:
            start = int(request.args.get("start",0))
            number = int(request.args.get("number",10))        
            name = request.args.get("name",None)
            origin = request.args.get("origin",None)
            year = request.args.get("year",None)
            collection = request.args.get("collection",None)
            dataset = request.args.get("dataset",None)
            parents = request.args.get("parents",None)
            offspring = request.args.get("offspring",None)
            condition_sql = "1"
            condition_variables = []
            if not name==None:
                condition_sql = condition_sql + " AND (`variety`.`name` LIKE ?)"
                condition_variables.append(name)
            if not origin==None:
                origin_list = origin.split(",")
                condition_sql = condition_sql + " AND (`variety`.`origin` IN ("+",".join(["?"]*len(origin_list))+"))"
                condition_variables.extend(origin_list)
            if not year==None:
                condition_sql = condition_sql + " AND NOT (`variety`.`year_max` IS NULL AND `variety`.`year_min` IS NULL)"
                if yearPattern.match(year):
                    condition_sql = condition_sql + " AND ((`variety`.`year_max` >= ?) OR (`variety`.`year_max` IS NULL))"
                    condition_sql = condition_sql + " AND ((`variety`.`year_min` <= ?) OR (`variety`.`year_min` IS NULL))"
                    condition_variables.extend([year,year])
                elif yearMaxPattern.match(year):
                    condition_sql = condition_sql + " AND (`variety`.`year_max` < ?)"
                    condition_variables.append(year[1:5])
                elif yearMinPattern.match(year):
                    condition_sql = condition_sql + " AND (`variety`.`year_min` > ?)"
                    condition_variables.append(year[1:5])
                elif yearRangePattern.match(year):
                    year_min = year[0:4]
                    year_max = year[5:9]
                    condition_sql = condition_sql + " AND ((`variety`.`year_min` >= ?) OR "
                    condition_sql = condition_sql + "((`variety`.`year_min` IS NULL) AND (`variety`.`year_max` >= ?)))"
                    condition_sql = condition_sql + " AND ((`variety`.`year_max` <= ?) OR "
                    condition_sql = condition_sql + "((`variety`.`year_max` IS NULL) AND (`variety`.`year_min` <= ?)))"
                    condition_variables.extend([year_min,year_min,year_max,year_max])
                else:
                    abort(422, "incorrect year condition "+str(year))
            if not collection==None:
                collection_list = collection.split(",")
                condition_sql = condition_sql + " AND (`collection`.`name` IN ("+",".join(["?"]*len(collection_list))+"))"
                condition_variables.extend(collection_list)
            if not dataset==None:
                if dataset=="any":
                    condition_sql = condition_sql + " AND NOT (`dataset`.`uid` IS NULL)"
                elif dataset=="none":
                    condition_sql = condition_sql + " AND (`dataset`.`id` IS NULL)"  
                elif dataset=="kmer":
                    condition_sql = condition_sql + " AND NOT (`dataset`.`location_kmer` IS NULL)"  
                elif dataset=="split":
                    condition_sql = condition_sql + " AND NOT (`dataset`.`location_split` IS NULL)"  
                else:
                    abort(422, "incorrect dataset condition "+str(dataset))
            if not parents==None:
                if _make_bool(parents):
                    condition_sql = condition_sql + " AND NOT (`ancestor`.`id` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`ancestor`.`id` IS NULL)"       
            if not offspring==None:
                if _make_bool(offspring):
                    condition_sql = condition_sql + " AND NOT (`offspring`.`id` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`offspring`.`id` IS NULL)"       
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT COUNT(DISTINCT `variety`.`id`) AS `number` \
                            FROM `variety` \
                            LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                            LEFT JOIN `variety_ancestor` AS `ancestor` ON `variety`.`uid` = `ancestor`.`variety` \
                            LEFT JOIN `variety_ancestor` AS `offspring` ON `variety`.`uid` = `offspring`.`ancestor` \
                            WHERE "+condition_sql, tuple(condition_variables))  
            total = cursor.fetchone()[0]
            if start<total:
                cursor.execute("SELECT `variety`.`uid`, `variety`.`name`, `variety`.`origin`, \
                                `country`.`name` AS `country`, `country`.`region`, \
                                `country`.`sub-region`, `country`.`intermediate-region`, \
                                NULL AS `year`, `variety`.`year_min`, `variety`.`year_max`, \
                                COUNT(DISTINCT `dataset`.`id`) as `datasets` \
                                FROM `variety` \
                                LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                                LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `variety_ancestor` AS `ancestor` ON `variety`.`uid` = `ancestor`.`variety` \
                                LEFT JOIN `variety_ancestor` AS `offspring` ON `variety`.`uid` = `offspring`.`ancestor` \
                                WHERE "+condition_sql+" \
                                GROUP BY `variety`.`id` \
                                ORDER BY `variety`.`name`, `variety`.`uid` \
                                LIMIT ?,?",tuple(condition_variables + [start,number]))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                resultList = []
            for i in range(len(resultList)):
                resultList[i] = adjust_variety_response(resultList[i], db_connection)                
            response = {"start": start, "number": number, "total": total, "list": resultList}
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>")
@namespace.doc(params={"uid": "unique identifier variety"})
class VarietyId(Resource):
    
    @namespace.doc(description="Get variety by uid")
    def get(self,uid):
        try:
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT `variety`.`uid`, `variety`.`name`, `variety`.`origin`, \
                            `country`.`name` AS `country`, `country`.`region`, \
                            `country`.`sub-region`, `country`.`intermediate-region`, \
                            NULL AS `year`, `variety`.`year_min`, `variety`.`year_max`, \
                            NULL as `datasets` \
                            FROM `variety` \
                            LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                            LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                            WHERE `variety`.`uid` = ? \
                            GROUP BY `variety`.`id`",(uid,))  
            data = cursor.fetchone()
            if data:
                response = adjust_variety_response(dict(data), db_connection)
                return Response(json.dumps(response), mimetype="application/json")
            else:
                abort(404, "no variety with uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))