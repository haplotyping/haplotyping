# Indexation process

The number of workers is computed automatically 
and depends on available memory and the 
maximum number of parallel processes.


```mermaid
graph TD;

  reads(("<b>reads</b><br/>fastq<br/>files"));
  reads-->queueAutomaton;

  subgraph Splits
    k["k"]
    k1["k' < k"]
    automaton(("<b>automaton</b><br/>splitting<br/>k-mers"))
    index(("<b>index</b><br/>splitting<br/>k-mers"))
    k1.->automaton
    k.->index
    automaton.->index
  end

  subgraph HDF5 Database
    hdf5Kmer["splitting k-mers"]
    hdf5Base["splitting bases"]
  end
  subgraph HDF5 Database
    hdf5Direct["direct connections"]
  end
  subgraph HDF5 Database
    hdf5Partition["partitioning"]
    hdf5ReadInfo["read info"]
    hdf5ReadData["read data"]
    hdf5Partition.->hdf5ReadInfo
    hdf5ReadInfo.->hdf5ReadData
  end

  mergeDirect-->hdf5Direct
  hdf5Kmer-->automaton
  hdf5Kmer-->index
  hdf5Direct--"optional"-->partition
  partition-->hdf5Partition
  mergeReads-->hdf5ReadInfo
  mergeReads-->hdf5ReadData

  subgraph Connections
    queueAutomaton(("queue<br/>automaton"))
    workerAutomaton["<b>workers</b> automaton"]
    queueIndex(("queue<br/>index"))
    workerIndex["<b>workers</b> index"]
    queueMatches(("queue<br/>matches"))
    workerMatches["<b>workers</b> matches"]
    storageMatches{"temporary<br/>storage read<br/>matches"}
    storageDirect{"temporary<br/>storage direct<br/>connections"}
    workerMergeDirect["<b>workers</b> merge direct"]
    storageMergedDirect{"temporary<br/>storage merged<br/>direct"}
    mergeDirect["merge and store direct"]
    partition["partition graph"]
    workerReads["<b>workers</b> reads"]
    storageReads{"temporary<br/>storage partitioned<br/>read matches"}
    workerMergeReads["<b>workers</b> merge reads"]
    storageMergedReads{"temporary<br/>storage merged<br/>partitioned reads"}
    mergeReads["merge and store reads"]
    automaton-->workerAutomaton
    index-->workerIndex
    queueAutomaton.->workerAutomaton
    workerAutomaton.->queueIndex
    queueIndex.->workerIndex
    workerIndex.->queueMatches
    queueMatches.->workerMatches
    workerIndex--"optional"-->storageMatches
    workerMatches-->storageDirect
    storageDirect-->workerMergeDirect
    workerMergeDirect-->storageMergedDirect
    storageMergedDirect-->mergeDirect
    partition-->workerReads
    storageMatches-->workerReads
    workerReads-->storageReads
    storageReads-->workerMergeReads
    workerMergeReads-->storageMergedReads
    storageMergedReads-->mergeReads
  end

```

