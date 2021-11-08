#! /bin/bash

seqdir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/
sigdir=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/signalp/
sigdir_aa=/home/troyalty/Documents/projects/alphafold_stability/data/M0059E/dbcan_annotation/signalp/aa/

rm -fr $sigdir 2> /dev/null
mkdir $sigdir
mkdir $sigdir_aa


signal_ex=_summary.signalp5
aa_ex=.faa
sig_ex=.signalp5
cleaved_ex=_cleaved.faa
tmp_ex=.tmp
dbcan_ex=_dbcan.tsv
dbcan_p_ex=_dbcan_p.tsv

for f in $seqdir*$aa_ex; do
	tmp1="${sigdir}tmp1"
	tmp2="${sigdir}tmp2"
	tmp3="${sigdir}tmp3"
	tmp4="${sigdir}tmp4"

	file=$(echo $f | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
	
	#predict presence of signal peptide with signalp
	signalp -fasta $f -org gram+ -format short -prefix $tmp1
	signalp -fasta $f -org gram- -format short -prefix $tmp2
	signalp -fasta $f -org arch -format short -prefix $tmp3
	signalp -fasta $f -org euk -format short -prefix $tmp4
	
	#delete header lines from signalp output files
	sed -i '1,2d' $tmp1$signal_ex
	sed -i '1,2d' $tmp2$signal_ex
	sed -i '1,2d' $tmp3$signal_ex
	sed -i '1,2d' $tmp4$signal_ex
	
	#select signal peptide prediction with highest probability
	cat $tmp1$signal_ex $tmp2$signal_ex $tmp3$signal_ex $tmp4$signal_ex |\
	grep -v "OTHER" |\
	awk -F '\t' 'BEGIN{OFS="\t"}{if (NF==5) m_i=5; else m_i=7; m=$3; for(i=3;i<m_i;i++) {if ($i>m) {m=$i}}; print $1,$2,m,$(NF)}' |\
	sort -r -k 1,3 |\
	rev |\
	uniq -f 8 |\
	rev > $sigdir$file$sig_ex
	
	#filter sequence ids for unambiguous dbcan annotations
	grep -f <(cut -f 1 $sigdir$file$sig_ex) $seqdir$file$dbcan_ex |\
	cut -f 1,3 |\
	sort -k 1,2 |\
	uniq |\
	sort -k 2 |\
	uniq -f 1 > $seqdir$file$dbcan_p_ex

	#extract sequences that have signal peptide
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 1 $sigdir$file$sig_ex) $f |\
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' |\
	sed '1d' > $sigdir_aa$file$cleaved_ex$tmp_ex

	#truncate signal peptides from sequences using signalp cleavage site predictions
	while read -r gene type pr cleave; do
		trn_pos=$(echo $cleave | cut -f 3 -d ' ' | cut -f 1 -d '-')
		grep -A 1 $gene $sigdir_aa$file$cleaved_ex$tmp_ex | sed "2s/^.\{$trn_pos\}//" >> $sigdir_aa$file$cleaved_ex
	
	done < $sigdir$file$sig_ex

	#remove temporary files
	rm $tmp1$signal_ex $tmp2$signal_ex $tmp3$signal_ex $tmp4$signal_ex $sigdir_aa$file$cleaved_ex$tmp_ex
done
