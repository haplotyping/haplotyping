from flask import Response, request, abort
from flask_restx import Namespace, Resource, fields
from flask_caching import Cache
import json, haplotyping, sqlite3, os

from haplotyping.service.kmer_kmc import Kmer as KmerKMC
from haplotyping.service.split import Split

def _make_cache_key(*args, **kwargs):
    #todo: find better solution    
    try:
        cacheKey = "%s_%s_%s_%s" % (args[0].__class__.__name__, str(request.path), str(request.args), str(namespace.payload))
    except:
        cacheKey = "%s_%s_%s" % (args[0].__class__.__name__, str(request.path), str(request.args))
    return cacheKey
    
namespace = Namespace("kmer", description="K-mer frequencies for a dataset", path="/kmer")
cache = Cache()
parser = namespace.parser()
            
def _getDataset(uid):
    db_connection = haplotyping.service.API.get_db_connection()
    db_connection.row_factory = sqlite3.Row
    cursor = db_connection.cursor()
    cursor.execute("SELECT `dataset`.`uid`, \
                    `dataset`.`type`, \
                    `dataset`.`location` AS `dataset_location`, \
                    `collection`.`location` AS `collection_location` \
                    FROM `dataset` \
                    LEFT JOIN `collection` ON `dataset`.`collection_id` = `collection`.`id` \
                    WHERE `dataset`.`uid` = ? \
                    AND (`dataset`.`type` = 'kmer' OR `dataset`.`type` = 'split') \
                    AND NOT `collection`.`location` IS NULL",(str(uid),))  
    data = cursor.fetchone()
    return data

@namespace.route("/<uid>/<kmer>")
class KmerSingle(Resource):
    
    dataset_kmers = parser.copy()
    dataset_kmers.add_argument("mismatches", type=int, required=False, location="values", 
                              help="maximum number of mismatches", choices=[0,1,2])
    
    @namespace.doc(description="Get k-mer frequencies from dataset defined by uid")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "k-mer"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,kmer):
        mm = min(2,max(0,int(request.args.get("mismatches",0))))
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.kmc")
                else:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.kmc")
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,[kmer],mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmer))
                elif kmc_query_binary_location:
                    kmc_query_binary_location_query = os.path.join(kmc_query_binary_location,"kmc_query")
                    kmc_query_binary_location_analysis = os.path.join(kmc_query_binary_location,"kmc_analysis")
                    response = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location_query,location_kmc,[kmer],mm)
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmer))
                    response["info"] = KmerKMC.kmc_binary_info(kmc_query_binary_location_analysis,location_kmc)
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/<kmer>/split")
class KmerSplit(Resource):
    
    dataset_kmers = parser.copy()
    dataset_kmers.add_argument("distance", type=int, required=False, location="values", 
                              help="maximum distance, default 1000")
    dataset_kmers.add_argument("minimumFrequency", type=int, required=False, location="values", 
                              help="minimum frequency, default 1")
    
    @namespace.doc(description="Get path from k-mer to first splitting k-mer")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset","kmer": "k-mer"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,kmer):
        distance = max(1,int(request.args.get("distance",1000)))
        minimumFrequency = max(1,int(request.args.get("minimumFrequency",1)))
        try:
            dsData = _getDataset(uid)
            if dsData:
                if not dsData["collection_location"]==None:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                dsData["collection_location"],dsData["dataset_location"],"kmer.kmc")
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                dsData["collection_location"],dsData["dataset_location"],"kmer.data.h5")
                else:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                dsData["dataset_location"],"kmer.kmc")
                    location_split = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                dsData["dataset_location"],"kmer.data.h5")
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                
                newKmer = kmer
                result = {"distance": None, "splittingKmers":[], "pathKmers": [newKmer]}
                for i in range(distance):
                    rightNeighbours = set()
                    leftSplitters = set()
                    for rb in ["A","C","G","T"]:
                        rightNeighbours.add(newKmer[1:]+rb)
                    if i>0:
                        for lb in ["A","C","G","T"]:
                            if not lb==newKmer[0]:
                                leftSplitters.add(lb+newKmer[1:])
                    kmerList = list(rightNeighbours.union(leftSplitters))
                    if kmc_query_library:
                        data = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmerList,0)
                        if not data:
                            abort(500,"no response using kmc query library for "+str(kmerList))
                    elif kmc_query_binary_location:
                        kmc_query_binary_location_query = os.path.join(kmc_query_binary_location,"kmc_query")
                        kmc_query_binary_location_analysis = os.path.join(kmc_query_binary_location,"kmc_analysis")
                        data = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location_query,location_kmc,kmerList,0)
                        if not data:
                            abort(500,"no response using kmc query binary for "+str(kmerList))
                    else:
                        abort(500,"no kmc binary or library configured")
                    #analyse result
                    kmerFound = [k for k in data["kmers"].keys() if data["kmers"][k]>=minimumFrequency]
                    if len([k for k in leftSplitters if k in kmerFound])>0:
                        result["splittingKmers"] = [newKmer]
                        result["distance"] = i
                        result["pathKmers"] = result["pathKmers"][:-1]
                        break               
                    else:
                        rightSplitters = [k for k in rightNeighbours if k in kmerFound or 
                                          haplotyping.General.reverse_complement(k) in kmerFound]
                        if len(rightSplitters)>1:
                            result["splittingKmers"] = list(rightSplitters)
                            result["distance"] = i+1
                            break
                        elif len(rightSplitters)==0:
                            return result
                        else:
                            newKmer = rightSplitters[0] 
                            result["pathKmers"].append(newKmer)  
                if len(result["splittingKmers"])>0:
                    checkedSplittingKmers = []
                    for splittingKmer in result["splittingKmers"]:
                        response = Split.kmer_info(location_split, splittingKmer)
                        if response:
                            checkedSplittingKmers.append(splittingKmer)
                    result["splittingKmers"] = checkedSplittingKmers
                return Response(json.dumps(result), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
            
@namespace.route("/<uid>/<kmer1>/path/<kmer2>")
class KmerPath(Resource):
    
    dataset_kmers = parser.copy()
    dataset_kmers.add_argument("distance", type=int, required=False, location="values", 
                              help="maximum distance, default 1000")
    dataset_kmers.add_argument("minimumFrequency", type=int, required=False, location="values", 
                              help="minimum frequency, default 1")
    
    @namespace.doc(description="Get path from k-mer 1 to k-mer 2")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset","kmer1": "k-mer","kmer2": "k-mer"})
    @cache.cached(make_cache_key=_make_cache_key)
    def get(self,uid,kmer1,kmer2):
        distance = max(1,int(request.args.get("distance",1000)))
        minimumFrequency = max(1,int(request.args.get("minimumFrequency",1)))
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.kmc")
                else:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.kmc")
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                #construct path
                newKmer = kmer1
                path = kmer1
                for i in range(distance):
                    rightNeighbours = set()
                    for rb in ["A","C","G","T"]:
                        rightNeighbours.add(newKmer[1:]+rb)
                    kmerList = list(rightNeighbours)                    
                    if kmc_query_library:
                        data = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmerList,0)
                        if not data:
                            abort(500,"no response using kmc query library for "+str(kmerList))
                    elif kmc_query_binary_location:
                        kmc_query_binary_location_query = os.path.join(kmc_query_binary_location,"kmc_query")
                        kmc_query_binary_location_analysis = os.path.join(kmc_query_binary_location,"kmc_analysis")
                        data = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location_query,location_kmc,kmerList,0)
                        if not data:
                            abort(500,"no response using kmc query binary for "+str(kmerList))
                    else:
                        abort(500,"no kmc binary or library configured")
                    #analyse result
                    kmerFound = [k for k in data["kmers"].keys() if data["kmers"][k]>=minimumFrequency]
                    rightSplitters = [k for k in rightNeighbours if k in kmerFound or 
                                          haplotyping.General.reverse_complement(k) in kmerFound]
                    if len(rightSplitters)==1:
                        newKmer = rightSplitters[0]
                        path = path + newKmer[-1]
                    elif kmer2 in rightSplitters:
                        for rightSplitter in rightSplitters:
                            if rightSplitter==kmer2:
                                path = path + rightSplitter[-1]                                
                    else:
                        return None
                    if kmer2==path[-len(kmer2):]:
                        return Response(json.dumps(path), mimetype="application/json")                    
                else:
                    return None
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
        
            
@namespace.route("/<uid>")
class KmerMultiple(Resource):            
                
    dataset_kmers = namespace.model("k-mer list to get frequencies", {
        "kmers": fields.List(fields.String, attribute="items", required=True, description="list of k-mers"),
        "mismatches": fields.Integer(min=0,max=2, required=False, description="maximum number of mismatches")
    })
                
    @namespace.doc(description="Get k-mer frequencies from dataset defined by uid")
    @namespace.expect(dataset_kmers)
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):
        mm = min(2,max(0,namespace.payload.get("mismatches",0)))
        kmers = namespace.payload.get("kmers",[])
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.kmc")
                else:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.kmc")
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmers))
                elif kmc_query_binary_location:
                    kmc_query_binary_location_query = os.path.join(kmc_query_binary_location,"kmc_query")
                    kmc_query_binary_location_analysis = os.path.join(kmc_query_binary_location,"kmc_analysis")
                    response = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location_query,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmers))
                    response["info"] = KmerKMC.kmc_binary_info(kmc_query_binary_location_analysis,location_kmc)
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
        
            
@namespace.route("/<uid>/sequence")
class KmerSequence(Resource):            
                
    sequence_data = namespace.model("sequence to get k-mer frequencies", {"sequence": 
                                            fields.String(required=True, description="sequence"),
                                    "mismatches": fields.Integer(min=0,max=2, required=False, 
                                                                 description="maximum number of mismatches")})
    
    @namespace.doc(description="Get k-mer frequencies for a sequence from dataset defined by uid")
    @namespace.doc(params={"uid": "unique identifier dataset"})
    @namespace.expect(sequence_data)
    @cache.cached(make_cache_key=_make_cache_key)
    def post(self,uid):  
        mm = min(2,max(0,namespace.payload.get("mismatches",0)))
        sequence = namespace.payload.get("sequence","")
        kmers = []
        try:
            data = _getDataset(uid)
            if data:
                if not data["collection_location"]==None:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["collection_location"],data["dataset_location"],"kmer.kmc")
                else:
                    location_kmc = os.path.join(haplotyping.service.API.get_data_kmer_location(),
                                                data["dataset_location"],"kmer.kmc")
                kmc_query_library = haplotyping.service.API.get_kmc_query_library()
                kmc_query_binary_location = haplotyping.service.API.get_kmc_query_binary_location()
                if kmc_query_library:
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,[],0)
                    if not response:
                        abort(500,"no response using kmc query library")
                    k = response["info"].get("kmer_length",0)
                    kmers = [sequence[i:i+k] for i in range(len(sequence)-(k-1)) if k>0]
                    response = KmerKMC.kmc_library(kmc_query_library,location_kmc,kmers,mm)
                    if not response:
                        abort(500,"no response using kmc query library for "+str(kmers))
                elif kmc_query_binary_location:
                    kmc_query_binary_location_query = os.path.join(kmc_query_binary_location,"kmc_query")
                    kmc_query_binary_location_analysis = os.path.join(kmc_query_binary_location,"kmc_analysis")
                    info = KmerKMC.kmc_binary_info(kmc_query_binary_location_analysis,location_kmc)
                    k = info.get("kmer_length",0)
                    kmers = [sequence[i:i+k] for i in range(len(sequence)-(k-1)) if k>0]
                    response = KmerKMC.kmc_binary_frequencies(kmc_query_binary_location_query,location_kmc,kmers,mm)
                    response["info"] = info;
                    if not response:
                        abort(500,"no response using kmc query binary for "+str(kmers))
                else:
                    abort(500,"no kmc binary or library configured")
                return Response(json.dumps(response), mimetype="application/json")                
            else:
                abort(404, "no dataset with k-mers for uid "+str(uid))
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))