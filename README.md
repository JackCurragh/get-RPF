# get-RPF
Extract just the Ribosome Protected Fragment (RPF) from Ribo-Seq Reads 

This method aims to return just the RPFs (minus adpaters/UMIs/Barcodes) from a given FASTQ file. 

The goal of the method is to find the number of bases that can be trimmed from the 5' and 3' ends of each read producing an alignable RPF. 

## Concept 
1. The method starts by searching for known adapters and trimming those that are detected from the reads. 
2. A subset of the resulting reads are aligned against the reference 
3. The resulting alignment stats and mapping qualities are used to determine how many bases can be trimmed to remove other ligated sequences such as barcodes or UMIs. 
