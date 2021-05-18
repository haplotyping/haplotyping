#created kmc files and sorted list from reads
find . -name "*.fastq.gz"  > kmer.lst 
kmc -k31 -ci1 -cx100000000 -cs1000000000 -t20 @"./kmer.lst" "./kmer.kmc" tmp
kmc_analysis dump "./kmer.kmc" "./kmer.list" -min 2 -max 999999999 -rc
sort "./kmer.list" > "./kmer.list.sorted"
gzip "./kmer.list.sorted"
rm "./kmer.list"

