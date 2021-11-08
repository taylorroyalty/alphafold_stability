#! /bin/bash

cpus=30
alphafold_db=/srv/data/alphafold
alphafold_path=/home/troyalty/Software/alphafold
proj_dir=/home/troyalty/Documents/projects/alphafold_stability
foldx_path=/home/troyalty/Software/foldx/foldx

aa_path=$proj_dir/data/dbCAN
tmp_path=$proj_dir/tmp/
foldx_results_path=$proj_dir/data/foldx/
alphafold_results_path=$proj_dir/data/alphafold_results

pdb_alphafold_path=$alphafold_results_path/pdb/

pdb_foldx_path=$foldx_results_path/pdb/

tmp_fa=$tmp_path/tmp.fa_ex
tmp_pdb=$tmp_path/tmp_pdb.txt

export $cpus
export $alphafold_db
export $alphafold_path
export $proj_dir
export $aa_path
export $tmp_path
export $alphafold_results_path
export $pdb_alphafold_path
export $pdb_foldx_path
export $tmp_file
export $tmp_pdb

rm -fr $tmp_path 2> /dev/null
mkdir $tmp_path

rm -fr $alphafold_results_path 2> /dev/null
mkdir $alphafold_results_path
mkdir $pdb_alphafold_path

rm -fr $foldx_results_path 2> /dev/null
mkdir $pdb_foldx_path
