from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
import json, haplotyping, sqlite3

from haplotyping.service.marker import Marker

namespace = Namespace("marker", description="Marker information for a dataset", path="/marker")
parser = namespace.parser()

def getDataset(uid):
    db_connection = haplotyping.service.API.get_db_connection()
    db_connection.row_factory = sqlite3.Row
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                    `collection`.`name` AS `collection`, \
                    `dataset`.`location_marker`, \
                    `dataset`.`marker_id`, \
                    `collection`.`location` AS `location` \
                    FROM `dataset` \
                    LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                    WHERE `dataset`.`uid` = ? \
                    AND NOT `dataset`.`location_marker` IS NULL \
                    AND NOT `dataset`.`marker_id` IS NULL \
                    AND NOT `collection`.`location` IS NULL",(str(uid),))  
    return cursor.fetchone()
    
            
@namespace.route("/<uid>/data")
class MarkerData(Resource):
    
    @namespace.doc(description="Get marker data from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    def get(self,uid):
        try:
            data = getDataset(uid)
            if data:
                location_marker = haplotyping.service.API.get_data_location() + data["location"] + data["location_marker"]
                marker_id = data["marker_id"]
                (ploidy,markerList,) = Marker.getMarkerData(location_marker,marker_id)
                response={"collection": data["collection"], "ploidy": ploidy, "data": markerList}
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with markers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/mapping")
class MarkerMapping(Resource):
    
    @namespace.doc(description="Get marker mapping for the collection containing the dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    def get(self,uid):
        try:
            data = getDataset(uid)
            if data:
                location_marker = haplotyping.service.API.get_data_location() + data["location"] + data["location_marker"]
                markerMapping = Marker.getMarkerMapping(location_marker)
                response={"collection": data["collection"], "mapping": markerMapping}
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with markers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            

@namespace.route("/<uid>/info")
class MarkerInfo(Resource):
    
    @namespace.doc(description="Get marker information like reference/alternative "+
                   "and flanking sequences for the collection containing the dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    def get(self,uid):
        try:
            data = getDataset(uid)
            if data:
                location_marker = haplotyping.service.API.get_data_location() + data["location"] + data["location_marker"]
                markerInfo = Marker.getMarkerInfo(location_marker)
                response={"collection": data["collection"], "info": markerInfo}
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with markers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
            