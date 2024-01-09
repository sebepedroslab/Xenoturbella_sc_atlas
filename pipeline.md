# Load R libraries and source code
```
library(metacell)
source("xboc_functions.R")
```

# Run metacell clustering
Optional. You can also start from the final metacell clustering object (provided in the "input_data/" folder, see below).

1. Direct mapping of MARS-seq reads, using FASTX-toolkit and STAR. Consisder using only a subset of reads. 

```
zcat scdb_path/raw_reads/*/orig_files/*_R1*.fastq.gz > all_reads.fastq

##Strip plate barcodes
fastx_trimmer -Q33 -f 10 -l 49 -i all_reads.fastq -o all_reads_trim.fastq

##Generate STAR genome index and map trimmed reads
#you may need to adjust --genomeChrBinNbits and --genomeSAindexNbases  if you genome is split in >5,000 scaffolds (check STAR manual)
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./spsXSTARindex/ --genomeFastaFiles spsX_scaffolds.fasta --sjdbGTFfile spsX_annotation.gtf --sjdbOverhang 45

##Mapping
STAR --runThreadN 64 --genomeDir spsXSTARindex/ --readFilesIn all_reads_trim.fastq --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outFilterMultimapNmax 20 --limitBAMsortRAM 50000000000 --outFilterMismatchNmax 5 --alignIntronMax 5000 --genomeLoad LoadAndRemove --readNameSeparator ' ' --outFileNamePrefix SpsX_R1_

##Generate coverage bigWig for IGV analysis
$HOME/bin/ucsc_tools/wigToBigWig SpsX_R1_Signal.Unique.str1.out.wig Nvec.chrom.sizes SpsX_MARS_plus.bw
$HOME/bin/ucsc_tools/wigToBigWig SpsX_R1_Signal.Unique.str2.out.wig Nvec.chrom.sizes SpsX_MARS_minus.bw
```
    
