import ctypes, tempfile, subprocess

class Kmer:
    
    def kmc_library(library: str, filename: str, kmers: list = [], mm: int = 0):
        def get_status(status):
            response = {}
            response["size_pre"] = status[0];
            response["size_suf"] = status[1];
            response["kmer_length"] = status[2];
            response["mode"] = status[3];
            response["suffix_counter_size"] = status[4];
            response["prefix_length"] = status[5];
            response["signature_length"] = status[6];
            response["min_count"] = status[7];
            response["max_count"] = status[8];
            response["total_kmers"] = status[9];
            response["both_strands"] = status[10];
            response["kmc_version"] = status[11];
            response["signature_map_size"] = status[12];
            response["signature_map_position"] = status[13];
            response["prefixes_list_size"] = status[14];
            response["prefixes_size"] = status[15];
            response["prefixes_position"] = status[16];
            response["suffix_size"] = status[17];
            response["suffix_record_size"] = status[18];
            response["suffixes_position"] = status[19];
            response["suffixes_size"] = status[20];
            return response
        
        def get_stats(stats):
            response = {}
            response["checked"] = stats[0];
            response["positive"] = stats[1];
            response["minimum"] = stats[2];
            response["maximum"] = stats[3];
            return response
        
        outputFile = None
        try:            
            lib = ctypes.cdll.LoadLibrary(library)
            n=len(kmers)
            kmers_bytes = [bytes(kmer, "utf-8") for kmer in kmers]
            kmers_array = (ctypes.c_char_p * (n+1))()
            kmers_array[:-1] = kmers_bytes
            stats = (ctypes.c_uint32*4)()
            status = (ctypes.c_uint64*21)()
            if mm==0:
                frequencies = (ctypes.c_uint32*n)()
                if lib.kmer_frequencies(bytes(filename, "utf-8"),kmers_array,n,
                                        frequencies,stats,status):
                    response = {"info": get_status(status), "stats": get_stats(stats), "kmers": {}}                
                    for i in range(n):
                        response["kmers"][kmers[i]] = frequencies[i]
                    return response
                else:
                    return None   
            else:
                outputFile = tempfile.NamedTemporaryFile()                
                if lib.kmer_frequencies_mm(bytes(filename, "utf-8"),kmers_array,n,mm,bytes(outputFile.name, "utf-8"),
                                           stats,status):
                    result = [item.decode("ascii").strip().split("\t") for item in outputFile.readlines()]
                    response = {"info": get_status(status), "stats": get_stats(stats), "kmers": {}}
                    for item in result:
                        response["kmers"][item[0]]=item[1]
                    return response
                else:
                    return None   
        except Exception as e:
            return None   
        finally:
            if outputFile:
                outputFile.close()

    
    def kmc_binary(binary: str, filename: str, kmers: list = [], mm: int = 0):
        inputFile = None
        try:        
            response = {"info":{}, "stats":{}, "kmers": {}}
            if len(kmers)>0:
                inputFile = tempfile.NamedTemporaryFile()
                with open(inputFile.name, "wb") as file:
                    for kmer in kmers:
                        line = str(kmer).strip()+"\n"
                        file.write(line.encode("utf-8"))               
                args = [binary, filename, "-mm", str(mm), "-f", inputFile.name]
                p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                result = p.communicate()
                if result[1]:
                    #raise Exception(result[1].decode("UTF-8").strip())
                    return None
                else:
                    #kmers            
                    answers = result[0].decode("UTF-8").strip().split("\n")
                    for i in range(len(answers)):
                        items = answers[i].strip().split("\t")
                        if len(items)>1:
                            number = int(items[1])
                            response["kmers"][items[0]] = number  
            #stats
            response["stats"]["checked"] = None
            response["stats"]["positive"] = len(response["kmers"])
            if response["stats"]["positive"]>0:
                response["stats"]["minimum"] = min(response["kmers"].values())
                response["stats"]["maximum"] = max(response["kmers"].values())
            else:
                response["stats"]["minimum"] = None
                response["stats"]["maximum"] = None
            return response
        except Exception as e:
            return None   
        finally:
            if inputFile:
                inputFile.close()
        