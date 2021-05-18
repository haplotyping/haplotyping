class General:
    
    def reverse_complement(kmer: str) -> str: 
        """Return the reverse complement of the provided k-mer"""
        map={"A" : "T", "C" : "G", "T" : "A", "G" : "C", "N" : "N"}
        try:
            rckmer = ""
            for base in kmer:
                rckmer += map[base]
            rckmer = rckmer[::-1]
            return rckmer[:]
        except:
            raise Exception("invalid k-mer")

    def canonical(kmer: str) -> str:
        """Return the canonical form of the provided k-mer"""
        rckmer = General.reverse_complement(kmer)
        if rckmer<kmer:
            return rckmer
        else:
            return kmer