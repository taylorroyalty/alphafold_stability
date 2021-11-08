#! /bin/bash

pdb_ex=.pdb
tsv_ex=.tsv
gap=_

proj_dir=/home/troyalty/Documents/projects/alphafold_stability
depth_chemistry=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/geochemistry/chemistry_for_stability_predictions_0m_16m.tsv
foldx_path=/home/troyalty/Software/foldx/foldx
tmp_path=$proj_dir/tmp
foldx_results_path=$proj_dir/data/M0059E/foldx
alphafold_results_path=$proj_dir/data/M0059E/alphafold_results

pdb_alphafold_path=$alphafold_results_path/pdb/

pdb_foldx_path=$foldx_results_path/pdb

tmp_fa=$tmp_path/tmp.fa_ex
tmp_pdb=$tmp_path/tmp_pdb.txt

rm -fr $foldx_results_path 2> /dev/null
mkdir $foldx_results_path
mkdir $pdb_foldx_path $foldx_results_path/stability_results 

find $pdb_alphafold_path -type f -name "ranked_0.pdb" | sed '/^$/d' > $tmp_pdb

pwd_o=$(pwd)
cd $pdb_foldx_path

while read line; do
	seq=$(echo $line | rev  | cut -f 2 -d '/' | rev)
	cp $line $seq$pdb_ex
done<$tmp_pdb


>$foldx_results_path/stability_results/all_stability.tsv

#while read -r dataset depth ionic ph temp; do
while read -r dataset depth ionic ph temp; do
#	stability_calculation=$foldx_results_path/stability_results/$dataset$tsv_ex
	stability_calculation=$foldx_results_path/stability_results/$dataset$gap$depth$tsv_ex
	echo -e gene'\t'total_energy'\t'backbone_Hbond'\t'sidechain_Hbond'\t'van_der_waals'\t' \
        	electrostatics'\t'solvation_polar'\t'solvation_hydrophobic'\t'van_der_walls_clashes'\t'entropy_side_chain'\t' \
		entropy_main_chain'\t'sloop_entropy'\t'mloop_entropy'\t'cis_bond'\t'torsional_clash'\t' \
		backbone_clash'\t'helix_dipole'\t'water_bridge'\t'disulfide'\t'electrostatic_kon'\t' \
		partial_covalent_bonds'\t'energy_ionisation'\t'entropy_complex'\t'residue_number > $stability_calculation


	for f in ./*.pdb; do
		seq=$(echo $f | rev | cut -f 1 -d '/' | rev)
		$foldx_path -c Stability --temperature $temp --pH $ph --ionStrength $ionic --pdb $seq --output-dir ../ --output-file tmp
		cat ../tmp_ST.fxout | sed 's/\.pdb//' | sed 's/\.\///' >> $stability_calculation
		rm ../tmp_ST.fxout
	done
#	cat $stability_calculation | sed '1d' | cut -f 1,2 |  sed "s/^/$dataset\t$depth\t/" >> $foldx_results_path/stability_results/all_stability.tsv
	cat $stability_calculation | sed '1d' | sed "s/^/$dataset\t$depth\t/" >> $foldx_results_path/stability_results/all_stability.tsv

done< <(sed '1d' $depth_chemistry)

sed -i '1i\event\tdepth_m\tgene\ttotal_energy\tbackbone_Hbond\tsidechain_Hbond\tvan_der_waals\telectrostatics\tsolvation_polar\tsolvation_hydrophobic\tvan_der_walls_clashes\tentropy_side_chain\tentropy_main_chain\tsloop_entropy\tmloop_entropy\tcis_bond\ttorsional_clash\tbackbone_clash\thelix_dipole\twater_bridge\tdisulfide\telectrostatic_kon\tpartial_covalent_bonds\tenergy_ionisation\tentropy_complex\tresidue_number' $foldx_results_path/stability_results/all_stability.tsv

cd $pwd_o

