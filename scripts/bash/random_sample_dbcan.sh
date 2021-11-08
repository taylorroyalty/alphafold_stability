#! /bin/bash

target_annotation=GH109.hmm
sample_n=40

seqfile=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/SRR7066493.faa
blastfile=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/SRR7066492.faa
seqdir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/
sigdir_aa=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/signalp/aa/
randir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/random_sample/
blastdir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/blast/

rm -fr $blastdir 2> /dev/null
mkdir $blastdir
#rm -fr $randir 2> /dev/null
#mkdir $randir


aa_ex=.faa
cleaved_ex=_cleaved.faa
dbcan_ex=_dbcan_p.tsv
ran_ex=_random
db_ex=_db
gap=_
blast_ex=_results.tsv

#for f in $seqdir*$aa_ex; do
	file=$(echo $seqfile | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
	blastfile=$(echo $blastfile | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
#	file=$(echo $f | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
	
	grep $target_annotation $seqdir$file$dbcan_ex | shuf | head -n $sample_n > $blastdir$file$ran_ex$dbcan_ex
	grep $target_annotation $seqdir$blastfile$dbcan_ex > $blastdir$blastfile$ran_ex$dbcan_ex
#	grep $target_annotation $seqdir$file$dbcan_ex | shuf | head -n $sample_n > $randir$file$ran_ex$dbcan_ex

	
	#extract sequences that have signal peptide
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 2 $blastdir$file$ran_ex$dbcan_ex) $sigdir_aa$file$cleaved_ex |\
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' |\
	sed '1d' > $blastdir$file$ran_ex$cleaved_ex
#	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 2 $randdir$file$ran_ex$dbcan_ex) $sigdir_aa$file$cleaved_ex |\
#	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' |\
#	sed '1d' > $randir$file$ran_ex$cleaved_ex

	makeblastdb -in $sigdir_aa$blastfile$cleaved_ex -dbtype prot -out $blastdir$blastfile$db_ex
	blastp -query $blastdir$file$ran_ex$cleaved_ex -db $blastdir$blastfile$db_ex -max_target_seqs 1 -outfmt 6 -out $blastdir$file$gap$blastfile$blast_ex
	
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 2 $blastdir$file$gap$blastfile$blast_ex | sort | uniq) <(perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 2 $blastdir$blastfile$ran_ex$dbcan_ex)) $sigdir_aa$blastfile$cleaved_ex |\
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' |\
	sed '1d' > $blastdir$blastfile$ran_ex$cleaved_ex
#done
