#! /bin/bash
#annotate_dbcan.sh

aadir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/genes/aa/
nucdir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/genes/nuc/
results_dir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/
dbcan_db=/srv/data/dbCAN_v10/db/dbCAN-HMMdb-V10.txt

parser=/srv/data/dbCAN_v10/scripts/hmmscan-parser.sh

tmpid=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/tmpid.txt

rm -fr $results_dir 2> /dev/null
mkdir $results_dir

aa_ex=.faa
nuc_ex=.fna
dbcan_ex=_dbcan.tsv
tmp_ex=.tmp

for f in $aadir*.faa; do
	file=$(echo $f | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
	
	hmmscan --cpu 30 --domtblout $results_dir$file$dbcan_ex$tmp_ex -E 1e-10 $dbcan_db $aadir$file$aa_ex
	bash $parser $results_dir$file$dbcan_ex$tmp_ex | awk '$5<1e-18&&$10>0.35' > $results_dir$file$dbcan_ex
	rm $results_dir$file$dbcan_ex$tmp_ex
	
	cut -f 3 $results_dir$file$dbcan_ex > $tmpid

	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $tmpid $aadir$file$aa_ex > $results_dir$file$aa_ex
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $tmpid $nucdir$file$nuc_ex > $results_dir$file$nuc_ex

done

rm $tmpid



