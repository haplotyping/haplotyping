from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
import json, haplotyping, sqlite3, os

from haplotyping.service.marker import Marker

namespace = Namespace("marker", description="Marker information for a dataset", path="/marker")
parser = namespace.parser()

def getDataset(uid):
    db_connection = haplotyping.service.API.get_db_connection()
    db_connection.row_factory = sqlite3.Row
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                    `dataset`.`type`, \
                    `dataset`.`internal_id`, \
                    `dataset`.`location` AS `dataset_location`, \
                    `collection`.`location` AS `collection_location` \
                    FROM `dataset` \
                    LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                    WHERE `dataset`.`uid` = ? \
                    AND `dataset`.`type` = 'marker' \
                    AND NOT `collection`.`location` IS NULL",(str(uid),))  
    data = cursor.fetchone()
    return data
    
            
@namespace.route("/<uid>/data")
class MarkerData(Resource):
    
    @namespace.doc(description="Get marker data from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    def get(self,uid):
        try:
            data = getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_marker = os.path.join(haplotyping.service.API.get_data_marker_location(),
                                                data["collection_location"],data["dataset_location"])
                else:
                    location_marker = os.path.join(haplotyping.service.API.get_data_marker_location(),
                                                data["dataset_location"])
                marker_id = data["internal_id"]
                (name,markerList,) = Marker.getMarkerData(location_marker,marker_id)
                response=markerList
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
                if not data["collection_location"]==None:
                    location_marker = os.path.join(haplotyping.service.API.get_data_marker_location(),
                                                data["collection_location"],data["dataset_location"])
                else:
                    location_marker = os.path.join(haplotyping.service.API.get_data_marker_location(),
                                                data["dataset_location"])
                markerMapping = Marker.getMarkerMapping(location_marker)
                response=markerMapping
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
                if not data["collection_location"]==None:
                    location_marker = os.path.join(haplotyping.service.API.get_data_marker_location(),
                                                data["collection_location"],data["dataset_location"])
                else:
                    location_marker = os.path.join(haplotyping.service.API.get_data_marker_location(),
                                                data["dataset_location"])
                markerInfo = Marker.getMarkerInfo(location_marker)
                response=markerInfo
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with markers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
            