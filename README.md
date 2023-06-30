# haplotyping

## Indexation process

See [haplotyping.index](haplotyping/index) for more details

```mermaid
graph TD;

    fastq["<b>READS</b><br/>fastq files"]; 
    detect["detect splitting k-mer in reads"];
    direct["store direct connections reads"];
    
    fastq --> kmc;
    kmerlist --> splitlist;

    fastq --> detect;

    subgraph KMC Software
        k["k-mer<br/>size"]
        minfreq["minimum<br/>k-mer frequency"]
        maxfreq["maximum<br/>k-mer frequency"]
        kmc(<b>kmc database</b>)  
        k.->kmc
        minfreq.->kmc 
        maxfreq.->kmc  
    end

    subgraph scripts
        kmcanalysis["<i>kmc_analysis</i>"]
        kmerlist["<b>k-mer list</b>"]
        kmcanalysis-->sort
        sort-->gzip
        gzip-->kmerlist
        kmc-->kmcanalysis
    end

    subgraph HDF5 Database
        splits["splitting k-mers"]
        bases["splitting bases"]
        direct["direct connections"]
        partitioning["partitioning"]
        indirect["indirect connections"]
        splits-.->partitioning
        partitioning-.->indirect
    end

    subgraph haplotyping.index
      subgraph Splits
        splitlist["get splitting k-mers and<br/>bases from sorted list"]
        automaton["automaton<br/>splitting k-mers"]
        splitlist-->automaton
        splitlist--"storage"-->splits
        splitlist--"storage"-->bases
      end
      subgraph Connections
        detect["detect splitting k-mer in reads"]
        partition["partition the graph"]
        read["sort indirect connections reads"]
        automaton-->detect
        detect--"optional"-->partition
        partition --> read
        detect--"storage"-->direct
        partition--"storage"-->partitioning
        read--"storage"-->indirect
      end
    end

```


