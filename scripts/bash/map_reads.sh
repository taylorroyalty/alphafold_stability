#!/bin/bash

mapdir=/home/troyalty/Documents/projects/alphafold_stability/data/map_counts/
seqfile=/home/troyalty/Documents/projects/alphafold_stability/data/dbCAN.fna
gtffile=/home/troyalty/Documents/projects/alphafold_stability/data/dbCAN.gtf
fastqdir=/home/troyalty/Documents/projects/protein_expression/data/fastq_transcripts/

rm -fr $mapdir 2> /dev/null
mkdir $mapdir

#file extensions
sam_ex=.sam
bam_ex=.bam
bam_sort_ex=.srt.bam
bed_ex=.bed
abundance_ex=_abundance.tsv
tmp=.tmp

fP_ex=_1P.fastq
bP_ex=_2P.fastq

basename=dbCAN_annotated

tmp_list=tmp_list.txt
tmp_id_list=tmp_ids.txt
tmp_fastq=$mapdir$tmp_list
tmp_id=$mapdir$tmp_id_list

grep -v "^#" $gtffile | cut -f 1,9 | cut -f 1 -d ";" | sed 's/ID=//' > $tmp_id
ls $fastqdir | cut -d '_' -f 1 | sort | uniq > $tmp_fastq

bowtie2-build -f $seqfile $mapdir$basename --threads 30

while read line; do

	bowtie2 -x $mapdir$basename -1 $fastqdir$line$fP_ex -2 $fastqdir$line$bP_ex -S $mapdir$line$sam_ex -p 30

	samtools view -@ 30 -bS $mapdir$line$sam_ex > $mapdir$line$bam_ex
	samtools sort -@ 30 -o $mapdir$line$bam_sort_ex $mapdir$line$bam_ex

	featureCounts -p -T 30 -s 0 -t CDS -a $gtffile -o $mapdir$line$abundance_ex $mapdir$line$bam_sort_ex
	sed '1,2d' $mapdir$line$abundance_ex | cut -f 2,7 > $mapdir$line$abundance_ex$tmp
	mv $mapdir$line$abundance_ex$tmp $mapdir$line$abundance_ex

done<$tmp_fastq

#rm $tmp_fastq $tmp_id
