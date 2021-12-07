#! /bin/bash

#filepath containing fasta files
faa_filepath=/home/troyalty/Documents/projects/alphafold_stability/data/sequence/structures/genes/richness_analysis
#filepath housing alphafold2 run_docker.py script
alphafold2_script=/home/troyalty/Software/alphafold2_docker/alphafold/docker/run_docker.py
#directory to output structure results
output_path=/home/troyalty/Documents/projects/alphafold_stability/data/sequence/structures/genes/richness_analysis/structures
#a temporary directory for writing temporary fasta files
tmp_path=/home/troyalty/Documents/projects/alphafold_stability/data/sequence/structures/tmp/
#directory containing alphafold2 reference data
reference_data=/srv/data/alphafold/

#fasta extension for temporary files
fa_ex=.faa

#this allows for using conda's activate command. Note, that conda's path is added to the PATH in ~/.profile. There is apparently some variability in how
#rstudio executes .bashrc, .bash_profile, .profile, etc. among different OS. GL...
source ~/.profile
source activate
conda activate alphafold_docker

rm -fr $tmp_path 2> /dev/null
mkdir $tmp_path

rm -fr $output_path 2> /dev/null
mkdir $output_path


#loop through list 
for f in $faa_filepath/*.faa; do
cat $f | while read -r header; do #the while loop iterates through headers and sequences--note this requires the fasta file to be single-line fasta format
read -r seq #this reads the sequence under the header
gene_name=$(echo $header | cut -f 1 -d ' ' | cut -f 2 -d '>')
seq_file=$tmp_path/$gene_name$fa_ex

#alphafold only takes one sequence at at time; this requires isolating sequences in multifasta files into temporary files
echo $header > $seq_file
echo $seq >> $seq_file

#run the alphafold2 run_docker.py script
python3 $alphafold2_script \
--fasta_paths=$seq_file \
--max_template_date=2020-05-14 \
--data_dir=$reference_data \
--output_dir=$output_path
done
done

rm -fr $tmp_path 2> /dev/null
