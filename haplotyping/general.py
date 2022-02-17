class General:
    
    map={"A" : "T", "C" : "G", "T" : "A", "G" : "C", "N" : "N"}
    
    def reverse_complement(kmer: str) -> str: 
        """Return the reverse complement of the provided k-mer"""
        try:
            rckmer = ""
            for base in kmer:
                rckmer += General.map[base]
            rckmer = rckmer[::-1]
            return rckmer[:]
        except Exception as e:
            raise Exception("invalid k-mer: "+str(kmer))

    def canonical_old(kmer: str) -> str:
        """Return the canonical form of the provided k-mer"""
        rckmer = General.reverse_complement(kmer)
        if rckmer<kmer:
            return rckmer
        else:
            return kmer
        
    def canonical(kmer: str) -> str:
        """Return the canonical form of the provided k-mer"""
        try:
            rckmer = ""
            do_check=True
            for base in kmer:
                rcbase = General.map[base]
                if do_check:
                    if base<rcbase:
                        return kmer
                    elif rcbase<base:
                        do_check = False
                rckmer += rcbase
            rckmer = rckmer[::-1]
            return rckmer[:]
        except Exception as e:
            raise Exception("invalid k-mer: "+str(kmer))
        
    def is_canonical(kmer: str) -> str:
        """Test if provided k-mer is canonical"""
        try:
            for base in kmer:                
                if base<map[base]:
                    return True
                elif map[base]<base:
                    return False
            return True
        except Exception as e:
            raise Exception("invalid k-mer: "+str(kmer))