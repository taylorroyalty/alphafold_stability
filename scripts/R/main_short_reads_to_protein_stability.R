library(tidyverse)
library(foreach)
library(doParallel)

source("scripts/R/assembly_annotation_functions.R")

conda.path <- "/home/troyalty/Software/bioinformatics/miniconda3/bin/conda"
conda.env <- "assembly_annotation"
accession.list <- c("SRR7066493")
dbcan.db <- 'data/database/dbCAN-HMMdb-V10.txt'
override <- TRUE

#Identify conda bin path
conda.list <- reticulate::conda_list(conda = conda.path) #locate all conda env
bin.path <- conda.list[conda.list[,1] %in% conda.env, 2] %>%
  gsub("/python","",.)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), bin.path, sep = ":")) #add specified conda.env to PATH

#create standardized directory
create_sequence_directories(override = override)

#download specified short read datasets
download_SRR(accession.list = accession.list, numCores = 2)

#trim adaptars and drop low quality reads
trim_reads(accession.list = accession.list)

#assemble short reads into contigs
assemble_reads(accession.list, megahit.cores=30, min.contig.length = 1000)

#predict contig ORFs
predict_ORFs(accession.list)

#make a dbcan database
make_database(dbcan.db)

#annotate genes with dbcan annotations
dbcan_annotation(accession.list=accession.list, db=dbcan.db)

#filter and isolate genes annotated with dbcan annotations  
filter_fasta_with_gene_list(accession.list = accession.list, dir.gene.list = "data/sequence/annotation/hmm/", ex_gene_list = "_dbCAN-HMMdb-V10.txt_parsed.tsv", col = 3)

#separate dbCAN-positive annotations in signal peptide-positive and negative sequences--cleave signal peptide
filter_signal_peptide(accession.list = accession.list, faa_ex="_filtered.faa", cleave = TRUE)
