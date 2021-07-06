import h5py

class Marker:
    
    def getMarkerData(filename: str, varietyId: int):
        with h5py.File(filename, mode="r") as h5file:
            markerTable = h5file.get("/marker")
            scoreTable = h5file.get("/value")
            dataTable = h5file.get("/data")
            mappingTable = h5file.get("/mapping")
            mappingDataTable = h5file.get("/mappingData")
            if varietyId>=0 and varietyId<dataTable.shape[0]:
                response = []
                for i in range(dataTable.shape[1]):
                    score = dataTable[varietyId][i]
                    if score>=0:
                        markerName = markerTable[i][0].decode("utf-8")
                        scoreValue = scoreTable[score][0].decode("utf-8")
                        entry = {"name": markerName, "score": scoreValue, "mappings": []}
                        for j in range(mappingTable.shape[0]):
                            mappingEntry = mappingDataTable[i,j]
                            if mappingEntry[1]>=0:
                                entry["mappings"].append({
                                  "name": mappingTable[j][0].decode("utf-8"),
                                  "type": ("genetic" if mappingTable[j][1]=="g" else "physical"),
                                  "lg": mappingEntry[0].decode("utf-8"),
                                  "pos": (mappingEntry[1] if mappingTable[j][1]=="g" else int(mappingEntry[1]))
                                })
                        response.append(entry)
                return response
            else:
                raise Exception("invalid identifier '"+str(varietyId)+"', data not available"+str(dataTable.shape))
        