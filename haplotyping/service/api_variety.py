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

def get_variety_datasets(uid, collection, dataset, db_connection):
    condition_variables = [uid]
    condition_sql = "`dataset`.`variety` = ? \
                        AND (NOT `dataset`.`uid` IS NULL) AND (NOT `dataset`.`type` IS NULL)"
    if not collection==None:
        collection_list = collection.split(",")
        condition_sql = condition_sql + " AND (`collection`.`uid` IN ("+",".join(["?"]*len(collection_list))+"))"
        condition_variables.extend(collection_list)
    if not dataset==None:
        if dataset=="any":
            condition_sql = condition_sql + " AND NOT (`dataset`.`id` IS NULL)"
        elif dataset=="none":
            condition_sql = condition_sql + " AND (`dataset`.`id` IS NULL)"  
        elif dataset=="kmer":
            condition_sql = condition_sql + " AND ((`dataset`.`type` IS 'kmer') OR (`dataset`.`type` IS 'split'))"  
        elif dataset=="split":
            condition_sql = condition_sql + " AND (`dataset`.`type` IS 'split')" 
        elif dataset=="marker":
            condition_sql = condition_sql + " AND (`dataset`.`type` IS 'marker')" 
        else:
            abort(422, "incorrect dataset condition "+str(dataset))
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                           `dataset`.`type`, \
                           `collection`.`uid` AS `collection_uid`, \
                           `collection`.`name` AS `collection_name`, \
                           `collection`.`type` AS `collection_type`, \
                           `collection`.`experiment` AS `collection_experiment` \
                        FROM `dataset` \
                        LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                        WHERE "+condition_sql+" ORDER BY `collection`.`name`, `dataset`.`uid`", 
                   tuple(condition_variables)) 
    return [dict(row) for row in cursor.fetchall()]  

def parents_tree(id,parentData):
    parents = []   
    for item in parentData:
        if item["offspring"]==id:
            if item["ancestor"]==None:
                parent = {"relation": item["type"],
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
    
def adjust_variety_response(item, collection, dataType, db_connection):
    #origin
    if item["origin"]:
        item["origin"] = {"uid": item["origin"], "country": item["country"]}
    else:
        del item["origin"]
    del item["country"]
    #year
    if item["year_min"] or item["year_max"]:
        item["year"] = {"description": year_description(item["year_min"],item["year_max"]), 
                        "min": item["year_min"], 
                        "max": item["year_max"]}
    else:
        del item["year"]
    del item["year_min"]
    del item["year_max"]
    if "synonyms" in item and not item["synonyms"]:
        del item["synonyms"]
    #datasets
    item["datasets"] = get_variety_datasets(item["uid"], collection, dataType, db_connection)
    for i in range(len(item["datasets"])):
        item["datasets"][i]["collection"] = {
            "uid": item["datasets"][i]["collection_uid"],
            "name": item["datasets"][i]["collection_name"],
            "type": item["datasets"][i]["collection_type"],
            "experiment": item["datasets"][i]["collection_experiment"]
        }
        del item["datasets"][i]["collection_uid"]
        del item["datasets"][i]["collection_name"]
        del item["datasets"][i]["collection_type"]
        del item["datasets"][i]["collection_experiment"]
    #parents
    item["parents"] = get_variety_parents(item["uid"], db_connection)
    #offspring
    item["offspring"] = get_variety_offspring(item["uid"], db_connection)
    #finished
    return item
        
@namespace.route("/")
class VarietyList(Resource):
    
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
                              help="variety has dataset from comma separated list of collection uids")
    variety_list.add_argument("dataType", type=str, required=False, location="args", 
                              help="variety has dataset (of specific type)", choices=["any","none","kmer","split","marker"])
    variety_list.add_argument("hasParents", type=bool, required=False, location="args", 
                              help="parent(s) known for this variety")
    variety_list.add_argument("hasOffspring", type=bool, required=False, location="args", 
                              help="offspring known for this variety")

    variety_set = namespace.model("uid list to get varieties", {"uids": fields.List(fields.String, attribute="items", 
                              required=True, description="list of uids")})

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
            dataType = request.args.get("dataType",None)
            hasParents = request.args.get("hasParents",None)
            hasOffspring = request.args.get("hasOffspring",None)
            condition_sql = "1"
            condition_variables = []
            if not name==None:
                condition_sql = condition_sql + " AND ((`variety`.`name` LIKE ?) OR (`condition_synonym`.`synonym` LIKE ?))"
                condition_variables.extend([name,name])
            if not origin==None:
                origin_list = origin.split(",")
                condition_sql = condition_sql + " AND (`variety`.`origin` IN ("+",".join(["?"]*len(origin_list))+"))"
                condition_variables.extend(origin_list)
            if not year==None:
                condition_sql = condition_sql + " AND (NOT (`variety`.`year_max` IS NULL AND `variety`.`year_min` IS NULL))"
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
                condition_sql = condition_sql + " AND (`collection`.`uid` IN ("+",".join(["?"]*len(collection_list))+"))"
                condition_variables.extend(collection_list)
            if not dataType==None:
                if dataType=="any":
                    condition_sql = condition_sql + " AND NOT (`dataset`.`id` IS NULL)"  
                elif dataType=="none":
                    condition_sql = condition_sql + " AND (`dataset`.`id` IS NULL)"  
                elif dataType=="kmer":
                    condition_sql = condition_sql + " AND ((`dataset`.`type` IS 'kmer') OR (`dataset`.`type` IS 'split'))"  
                elif dataType=="split":
                    condition_sql = condition_sql + " AND (`dataset`.`type` IS 'split')" 
                elif dataType=="marker":
                    condition_sql = condition_sql + " AND (`dataset`.`type` IS 'marker')"
                else:
                    abort(422, "incorrect dataType condition "+str(dataset))
            if not hasParents==None:
                if _make_bool(hasParents):
                    condition_sql = condition_sql + " AND NOT (`ancestor`.`id` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`ancestor`.`id` IS NULL)"       
            if not hasOffspring==None:
                if _make_bool(hasOffspring):
                    condition_sql = condition_sql + " AND NOT (`offspring`.`id` IS NULL)"
                else:
                    condition_sql = condition_sql + " AND (`offspring`.`id` IS NULL)"     
            db_connection = haplotyping.service.API.get_db_connection()
            db_connection.row_factory = sqlite3.Row
            cursor = db_connection.cursor()
            cursor.execute("SELECT COUNT(DISTINCT `variety`.`id`) AS `number` \
                            FROM `variety` \
                            LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                            AND NOT `dataset`.`uid` IS NULL AND NOT `dataset`.`type` IS NULL \
                            LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                            LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                            LEFT JOIN `variety_ancestor` AS `ancestor` ON `variety`.`uid` = `ancestor`.`variety` \
                            LEFT JOIN `variety_ancestor` AS `offspring` ON `variety`.`uid` = `offspring`.`ancestor` \
                            LEFT JOIN `variety_synonym` AS `synonym` ON `variety`.`uid` = `synonym`.`uid` \
                            LEFT JOIN `variety_synonym` AS `condition_synonym` ON `variety`.`uid` = \
                                        `condition_synonym`.`uid` \
                            WHERE "+condition_sql, tuple(condition_variables))  
            total = cursor.fetchone()[0]
            if start<total:
                cursor.execute("SELECT `variety`.`uid`, `variety`.`name`, `variety`.`origin`, \
                                `country`.`name` AS `country`, \
                                NULL AS `year`, `variety`.`year_min`, `variety`.`year_max`, \
                                GROUP_CONCAT(DISTINCT `synonym`.`synonym`) AS `synonyms`, \
                                COUNT(DISTINCT `dataset`.`id`) as `datasets` \
                                FROM `variety` \
                                LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                                AND NOT `dataset`.`uid` IS NULL AND NOT `dataset`.`type` IS NULL\
                                LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `variety_ancestor` AS `ancestor` ON `variety`.`uid` = `ancestor`.`variety` \
                                LEFT JOIN `variety_ancestor` AS `offspring` ON `variety`.`uid` = `offspring`.`ancestor` \
                                LEFT JOIN `variety_synonym` AS `synonym` ON `variety`.`uid` = `synonym`.`uid` \
                                LEFT JOIN `variety_synonym` AS `condition_synonym` ON `variety`.`uid` = \
                                            `condition_synonym`.`uid` \
                                WHERE "+condition_sql+" \
                                GROUP BY `variety`.`id` \
                                ORDER BY `variety`.`name`, `variety`.`uid` \
                                LIMIT ?,?",tuple(condition_variables + [start,number]))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                resultList = []
            for i in range(len(resultList)):
                resultList[i] = adjust_variety_response(resultList[i], collection, dataType, db_connection)                
            response = {"start": start, "number": number, "total": total, "list": resultList}
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
    @namespace.doc(description="Get multiple varieties by uid")    
    @namespace.expect(variety_set)
    def post(self):
        uids = namespace.payload.get("uids",[])
        try:
            if len(uids)>0:
                condition_sql = "`variety`.`uid` IN ("+(",".join(["?"]*len(uids)))+")"
                condition_variables = [] + uids;
                db_connection = haplotyping.service.API.get_db_connection()
                db_connection.row_factory = sqlite3.Row
                cursor = db_connection.cursor()
                cursor.execute("SELECT COUNT(DISTINCT `variety`.`id`) AS `number` \
                                FROM `variety` \
                                LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                                AND NOT (`dataset`.`uid` IS NULL) AND NOT (`dataset`.`type` IS NULL) \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `variety_synonym` AS `synonym` ON `variety`.`uid` = `synonym`.`uid` \
                                WHERE "+condition_sql, tuple(condition_variables))  
                total = cursor.fetchone()[0]
                cursor.execute("SELECT `variety`.`uid`, `variety`.`name`, `variety`.`origin`, \
                                `country`.`name` AS `country`, \
                                NULL AS `year`, `variety`.`year_min`, `variety`.`year_max`, \
                                GROUP_CONCAT(DISTINCT `synonym`.`synonym`) AS `synonyms`, \
                                COUNT(DISTINCT `dataset`.`id`) as `datasets` \
                                FROM `variety` \
                                LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                                AND NOT (`dataset`.`uid` IS NULL) AND NOT (`dataset`.`type` IS NULL) \
                                LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                                LEFT JOIN `variety_synonym` AS `synonym` ON `variety`.`uid` = `synonym`.`uid` \
                                WHERE "+condition_sql+" \
                                GROUP BY `variety`.`id` \
                                ORDER BY `variety`.`name`, `variety`.`uid`",tuple(condition_variables))  
                resultList = [dict(row) for row in cursor.fetchall()]
            else:
                total = 0
                resultList = []
            for i in range(len(resultList)):
                resultList[i] = adjust_variety_response(resultList[i], None, None, db_connection)                
            response = {"total": total, "list": resultList}
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
                            `country`.`name` AS `country`, \
                            GROUP_CONCAT(DISTINCT `synonym`.`synonym`) AS `synonyms`, \
                            NULL AS `year`, `variety`.`year_min`, `variety`.`year_max`, \
                            NULL as `datasets` \
                            FROM `variety` \
                            LEFT JOIN `dataset` ON `variety`.`uid` = `dataset`.`variety` \
                            AND NOT (`dataset`.`uid` IS NULL) AND NOT (`dataset`.`type` IS NULL) \
                            LEFT JOIN `country` ON `variety`.`origin` = `country`.`uid` \
                            LEFT JOIN `variety_synonym` AS `synonym` ON `variety`.`uid` = `synonym`.`uid` \
                            WHERE `variety`.`uid` = ? \
                            GROUP BY `variety`.`id`",(uid,))  
            data = cursor.fetchone()
            if data:
                response = adjust_variety_response(dict(data), None, None, db_connection)
                return Response(json.dumps(response), mimetype="application/json")
            else:
                abort(404, "no variety with uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))