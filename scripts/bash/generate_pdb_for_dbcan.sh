#! /bin/bash

cpus=30
alphafold_db=/srv/data/alphafold/
alphafold_path=/home/troyalty/Software/alphafold/
aa_path=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/blast/
tmp_path=/home/troyalty/Documents/projects/alphafold_stability/tmp/
results_path=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/alphafold_results
pdb_path=$results_path/pdb/

fa_ex=.faa
tmp_file=tmp$fa_ex

rm -fr $tmp_path 2> /dev/null
mkdir $tmp_path

rm -fr $results_path 2> /dev/null
mkdir $results_path
mkdir $pdb_path

o_pwd=$(pwd)
cd $alphafold_path

for f in $aa_path*.faa; do
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $f | sed '1d' > $tmp_path$tmp_file
	filename=$(echo $f | cut -f 1 -d '.' | rev | cut -f 1 -d '/' | rev)
	mkdir $pdb_path$filename
	while read -r header; do
		read -r seq
		gene_name=$(echo $header | cut -f 1 -d ' ' | cut -f 2 -d '>')
		seq_file=$tmp_path/$gene_name$fa_ex
		echo $header > $seq_file
		echo $seq >> $seq_file
		bash $alphafold_path/run_alphafold.sh -d $alphafold_db -o $pdb_path$filename  -m model_1 -f $seq_file -t 2025-05-14 -n $cpus
		mv $seq_file $pdb_path$filename/$gene_name/
	done<$tmp_path$tmp_file

done

cd $o_pwd
