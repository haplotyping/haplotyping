import h5py, numpy as np

class Marker:
    
    def getMarkerData(filename: str, varietyId: int):
        with h5py.File(filename, mode="r") as h5file:
            varietyTable = h5file.get("/variety")
            dataTable = h5file.get("/data")
            if varietyId>=0 and varietyId<dataTable.shape[0]:
                markerList = []
                name = varietyTable["name"][varietyId].decode("utf8")
                values = dataTable[varietyId]
                for i in range(len(values)):
                    if values[i]>=0:
                        markerList.append(int(values[i]))
                    else:
                        markerList.append(None)
                return (name,markerList,)
            else:
                raise Exception("invalid identifier '"+str(varietyId)+"', data not available"+str(dataTable.shape))
                
    def getMarkerMapping(filename: str):
        with h5py.File(filename, mode="r") as h5file:
            markerTable = h5file.get("/marker")
            mappingTable = h5file.get("/mapping")
            lgTable = h5file.get("/lg")
            markerMapping = [] 
            markerNames = markerTable["name"]
            lgNames = lgTable["name"]
            mappings = {}
            for mapping in mappingTable.dtype.fields.keys():
                mappings[mapping] = [mappingTable[mapping]["lg"],mappingTable[mapping]["pos"]]
            for i in range(len(markerNames)):
                item = {
                    "name": markerNames[i].decode("utf-8"),
                    "mappings": {}
                }
                for mapping in mappings.keys():
                    lg = lgNames[mappings[mapping][0][i]].decode("utf-8")
                    pos = mappings[mapping][1][i]
                    if isinstance(pos, np.integer):
                        pos = int(pos)
                    elif isinstance(pos, np.floating):
                        pos = float(pos)
                    else:
                        pos = None
                    item["mappings"][mapping] = {"lg": lg,"pos": pos}
                markerMapping.append(item)            
            return markerMapping
        
    def getMarkerInfo(filename: str):
        with h5py.File(filename, mode="r") as h5file:
            markerTable = h5file.get("/marker")
            markerInfo = [] 
            infoData = {}
            markerNames = markerTable["name"]
            for key in markerTable.dtype.fields.keys():
                if not key=="name":
                    infoData[key] = markerTable[key]
            for i in range(len(markerNames)):
                item = {
                    "name": markerNames[i].decode("utf-8"),
                    "info": {}
                }
                for key in infoData.keys():
                    item["info"][key] = infoData[key][i].decode("utf-8")
                markerInfo.append(item)            
            return markerInfo
                
        