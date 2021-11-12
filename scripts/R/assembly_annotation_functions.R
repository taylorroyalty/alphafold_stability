#create standardized directory format for downstream processing.
create_subdirectories <- function(base.layer="all",override=FALSE){
  
  if (!base.layer %in% c("all","fastq.raw","fastq.trim","sequence","genes","hmm","database")){
    stop("The specified base.layer value is not an option.\n\nPlease select from: fastq.raw, fastq.trim, sequence, genes, hmm, or database.")
  }
  
  if (override == TRUE){
    unlink("data/sequence",recursive = TRUE)
  }
  
  if (base.layer %in% "all"){
    create <- TRUE 
  } else {
    create <- FALSE
  }
  
  dir.fastq.raw <- "data/sequence/fastq/raw/"
  dir.fastq.trim <- "data/sequence/fastq/trimmed/"
  dir.sequence <- "data/sequence/sequences/"
  dir.genes.aa <- "data/sequence/genes/aa/"
  dir.intracellular <- "data/sequence/genes/aa/signalp/intracellular/"
  dir.extracellular <- "data/sequence/genes/aa/signalp/extracellular/"
  dir.genes.nuc <- "data/sequence/genes/nuc/"
  dir.genes.gff <- "data/sequence/genes/gff/"
  dir.annotation.hmm <- "data/sequence/annotation/hmm/"
  dir.database <- "data/database/"
  
  warning('If not using fastq-dump, genomic data will need to be manually added to whatever the user-defined base layer is.')
  
  if ((!dir.exists(dir.fastq.raw) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "fastq.raw")){
    create <- TRUE
    dir.create(dir.fastq.raw,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.fastq.trim) & create == TRUE) | (!dir.exists(dir.fastq.trim) & base.layer %in% "fastq.trim")){
    create <- TRUE
    dir.create(dir.fastq.trim,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.sequence) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "sequence")){
    create <- TRUE
    dir.create(dir.sequence,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.genes.aa) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "genes")){
    create <- TRUE
    dir.create(dir.intracellular,recursive = TRUE)
    dir.create(dir.extracellular,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.genes.nuc) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "genes")){
    create <- TRUE
    dir.create(dir.genes.nuc,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.genes.gff) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "genes")){
    create <- TRUE
    dir.create(dir.genes.gff,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.annotation.hmm) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "hmm")){
    create <- TRUE
    dir.create(dir.annotation.hmm,recursive = TRUE)
  }
  
  if ((!dir.exists(dir.database) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "database")){
    create <- TRUE
    dir.create(dir.database,recursive = TRUE)
  }
}

#fastq-dump
#download sequence short read sequences from NCBI's SRA
download_SRR <- function(accession.list, dir.fastq.raw="", numCores=1){
  
  registerDoParallel(numCores)
  
  if (dir.fastq.raw %in% ""){
    dir.fastq.raw <- "./data/sequence/fastq/raw/"
  }
  
  if (!dir.exists(dir.fastq.raw)){
    stop("The output directory for fastq-dump does not exist.\n\nConsider executing the create_subdirectories function first.")
  }
  
  n.accession <- length(accession.list)
  foreach(i=1:n.accession) %dopar% {
    cmd.download <- sprintf("fastq-dump -O %s --split-files %s",dir.fastq.raw,accession.list[i])
    system(cmd.download)
  }
}

#fastp
#trim and quality control short read sequences
trim_reads <- function(accession.list, dir.fastq.trim="", dir.fastq.raw="", fastp.cores=16){
  
  #check to make sure user did exceed fastp's core limit 
  if (fastp.cores > 16){
    fastp.cores <- 6
    warning('Fastp caps the number of threads at 16.\n\n fastp.cores has been set to 16.')
  }
  
  if (dir.fastq.trim %in% ""){
    dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
  }
  
  if (!dir.exists(dir.fastq.trim)){
    stop("The output directory for trimmomatic does not exist.\n\nConsider executing the create_subdirectories function first.")
  }
  
  
  if (dir.fastq.raw %in% ""){
    dir.fastq.raw <- "./data/sequence/fastq/raw/"
  }
  
  if (!dir.exists(dir.fastq.raw) | (length(list.files(dir.fastq.raw)) == 0 )){
    stop("The directory containing raw short read sequences either does not exist or has no files.")
  }
  
  outpath.list <- paste0(dir.fastq.trim,accession.list)
  inpath.list <- paste0(dir.fastq.raw,accession.list)
  
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    cmd.trim <- sprintf("fastp -j %s.json -h %s.html -i %s_1.fastq -I %s_2.fastq -o %s_1.fastq -O %s_2.fastq -w %s",
                        outpath.list[i],
                        outpath.list[i],
                        inpath.list[i],
                        inpath.list[i],
                        outpath.list[i],
                        outpath.list[i],
                        fastp.cores)
    system(cmd.trim)
  }
}

#megahit
#assembly high quality reads into contigs
#note this is for paired reads only
assemble_reads <- function(accession.list, dir.fastq.trim="", dir.assembly="", megahit.cores=30, min.contig.length = 1000){
  
  if (dir.assembly %in% ""){
    dir.assembly <- "./data/sequence/assemblies/"
  }
  
  if (!dir.exists(dir.assembly)){
    stop("The output directory for megahit does not exist.\n\nConsider executing the create_subdirectories function first.")
  }
  
  if (dir.fastq.trim %in% ""){
    dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
  }
  
  if (!dir.exists(dir.fastq.trim) | (length(list.files(dir.fastq.trim)) == 0 )){
    stop("The directory containing raw short read sequences either does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.fastq.trim,accession.list)
  outdirs <- paste0(dir.assembly,accession.list)
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    if (!file.exists(paste0(inpath.list[i],'/',accession.list[i],'_1.fastq')) | !file.exists(paste0(inpath.list[i],'/',accession.list[i],'_2.fastq'))){
      warning(sprintf("One or both paired fastq files for: '%s' is missing",accession.list[i]))
      next
    }
    
    cmd.assemble <- sprintf("megahit -1 %s_1.fastq -2 %s_2.fastq --presets meta-sensitive -t %s -o %s --out-prefix %s --min-contig-len %s",
                            inpath.list[i],
                            inpath.list[i],
                            megahit.cores,
                            outdirs[i],
                            accession.list[i],
                            min.contig.length)
    system(cmd.assemble)
  }
}

#prodigal
#predict open reading frames; extension corresponds to genome extension. 
#The default corresponds to the megahit's default extension
predict_ORFs <- function(accession.list, dir.assembly="", dir.genes.aa="", dir.genes.nuc="", dir.genes.gff="", extension=".contigs.fa"){
  
  if (dir.genes.aa %in% ""){
    dir.genes.aa <- "data/sequence/genes/aa/"
  }
  
  if (dir.genes.nuc %in% ""){
    dir.genes.nuc <- "data/sequence/genes/nuc/"
  }
  
  if (dir.genes.gff %in% ""){
    dir.genes.gff <- "data/sequence/genes/gff/"
  }
  
  if (!dir.exists(dir.genes.aa) | !dir.exists(dir.genes.nuc) | !dir.exists(dir.genes.gff)){
    stop("The output directory for genes does not exist.\n\nConsider executing the create_subdirectories function first.")
  }
  
  if (dir.assembly %in% ""){
    dir.assembly <- "./data/sequence/assemblies/"
  }
  
  if (!dir.exists(dir.assembly) | (length(list.files(dir.assembly)) == 0 )){
    stop("The directory containing assemblies/genomes either does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.assembly,accession.list)
  outdir.aa <- paste0(dir.genes.aa,accession.list)
  outdir.nuc <- paste0(dir.genes.nuc,accession.list)
  outdir.gff <- paste0(dir.genes.gff,accession.list)
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    inpath.tmp <- paste0(inpath.list[i],'/',accession.list[i],extension)
    
    if (!file.exists(inpath.tmp)){
      warning(sprintf("The file: '%s' does not exist.",inpath.tmp))
      next
    }
    
    cmd.ORF <- sprintf("prodigal -q -f gff -i %s -a %s.faa -d %s.fna -o %s.gff",
                       inpath.tmp, #genomes/assemblies
                       outdir.aa[i], #amino acid sequences
                       outdir.nuc[i], #nucleic acid sequences
                       outdir.gff[i]) #generic feature format
    system(cmd.ORF)
    
    
    cmd.clean.aa <- sprintf("sed -i 's/\\*//g;s/^>/>%s\\_/' %s",
                            accession.list[i],
                            outdir.aa[i]) #remove astrieks and replace contig name with accession using sed
    
    cmd.clean.nuc <- sprintf("sed -i 's/^>/>%s\\_/' %s",
                             accession.list[i],
                             outdir.nuc[i]) #replace contig name with accession using sed
    system(cmd.clean.aa)
    system(cmd.clean.nuc)
  }
}

#hmmpress or makeblastdb
#Make databases for annotation. new functions (e.g., hmmbuild) could be added in the future
make_database <- function(file,fun='hmmpress',type='prot') {
  
  # dir.database <- "./data/database"
  # if (!dir.exists(dir.database)){
  #   stop("The output directory for databases does not exist.\n\nConsider executing the create_subdirectories function first.")
  # }
  
  warning("This generates the database in the same directory as the file.\n\nIf you do not have write permissions, please considering copying 'file' to a new location.")
  
  cmd.makedb <- switch(fun,
                       'hmmpress' = sprintf("hmmpress -f %s",file),
                       'makeblastdb' = sprintf("makeblastdb -dbtype %s -in %s -out %s",type,file,file))
  
  system(cmd.makedb)
  
}

#hmmscan
#annotate genes
annotate_hmmscan <- function(accession.list, dir.genes.aa="", dir.annotation.hmm="", db, hmm.cores=30, evalue=1e-10){
  
  if (dir.annotation.hmm %in% ""){
    dir.annotation.hmm <- "./data/sequence/annotation/hmm/"
  }
  
  if (!dir.exists(dir.annotation.hmm)){
    stop("The output directory for hmm annotations does not exist.\n\nConsider executing the create_subdirectories function first.")
  }
  
  if (dir.genes.aa %in% ""){
    dir.genes.aa <- "./data/sequence/genes/aa/"
  }
  
  if (!dir.exists(dir.genes.aa) | (length(list.files(dir.genes.aa)) == 0 )){
    stop("The output directory for amino acid fasta files does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.genes.aa,accession.list,".faa")
  outdir.list <- paste0(dir.annotation.hmm,accession.list,'_',basename(db),'.tsv')
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    if (!file.exists(inpath.list[i])){
      warning(sprintf("The amino acid fasta file for: '%s' is missing.",accession.list[i]))
      next
    }
    
    cmd.hmmscan <- sprintf("hmmscan --cpu %s --domtblout %s -E %s %s %s",
                           hmm.cores,
                           outdir.list[i],
                           evalue,
                           db,
                           inpath.list[i])
    system(cmd.hmmscan) #can't suppress terminal output from hmmscan
  }
  
  
}

#uses annotate_hmmscan and applies hmm-parser.sh function on the output as is done on dbcan's server
#annotate genes
dbcan_annotation <- function(accession.list, dir.genes.aa="", dir.annotation.hmm="", db, hmm.cores=30, evalue=1e-10){
  
  if (dir.annotation.hmm %in% ""){
    dir.annotation.hmm <- "./data/sequence/annotation/hmm/"
  }
  
  annotate_hmmscan(accession.list, dir.genes.aa=dir.genes.aa, dir.annotation.hmm=dir.annotation.hmm, db, hmm.cores, evalue)
  
  hmm.parser <- './scripts/bash/hmmscan-parser.sh'
  if (!file.exists(hmm.parser)){
    stop('This function uses a parsing file published for dbCAN. The filepath for this script should be ./scripts/bash/hmm-parser.sh\n\nThis script can be downloaded from: https://bcb.unl.edu/dbCAN2/download/Databases/V10/hmmscan-parser.sh\n\nNote, this verion corresponds to dbCAN v10')
  }
  
  inpath.list <- paste0(dir.annotation.hmm,accession.list,'_',basename(db),'.tsv')
  outpath.list <- paste0(dir.annotation.hmm,accession.list,'_',basename(db),'_parsed','.tsv')
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    cmd.clean.dbcan <- sprintf("bash %s %s > %s",hmm.parser,inpath.list[i],outpath.list[i])
    
    system(cmd.clean.dbcan)
  }
  
}

#signalp v5.0 
#identifies putative signal peptides and provides the option to cleave at predicted cleavage sites
#Note that signalp is stored in the local /usr/local/bin and is in the usr PATH (for marie). 
filter_signal_peptide <- function(accession.list, dir.genes.aa="", dir.extracellular="", dir.intracellular="", dir.tmp="", faa_ex=".faa", cleave = FALSE){
  
  if (dir.extracellular %in% ""){
    dir.extracellular <- "./data/sequence/genes/aa/signalp/extracellular/"
  }
  
  if (dir.intracellular %in% ""){
    dir.intracellular <- "./data/sequence/genes/aa/signalp/intracellular/"
  }
  
  if (!dir.exists(dir.extracellular) | !dir.exists(dir.intracellular)){
    stop("The output directories for signalp do not exist.\n\nConsider executing the create_subdirectories function first.")
  }
  
  if (dir.genes.aa %in% ""){
    dir.genes.aa <- "./data/sequence/genes/aa/"
  }
  
  if (!dir.exists(dir.genes.aa) | (length(list.files(dir.genes.aa)) == 0 )){
    stop("The output directory for amino acid fasta files does not exist or has no files.")
  }
  
  
  #make a temporary directory to house all temporary files generated from signalp
  if (dir.tmp %in% ""){
    dir.tmp <- "./data/sequence/genes/aa/signalp/tmp/"
  }
  
  unlink(dir.tmp,recursive = TRUE)
  dir.create(dir.tmp,recursive = TRUE)
  outdir.tmp <- paste0(dir.tmp,c("tmp_gram_pos","tmp_gram_neg","tmp_arc","tmp_euk")) #temporary files for each organism in sigalp
  
  inpath.list <- paste0(dir.genes.aa,accession.list,".faa")
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    cmd.signalp_gram_pos <- sprintf("signalp -fasta %s -org gram+ -format short -prefix %s",
                                    inpath.list[i],
                                    outdir.tmp[1])
    cmd.signalp_gram_neg <- sprintf("signalp -fasta %s -org gram- -format short -prefix %s",
                                    inpath.list[i],
                                    outdir.tmp[2])
    cmd.signalp_gram_arc <- sprintf("signalp -fasta %s -org arch -format short -prefix %s",
                                    inpath.list[i],
                                    outdir.tmp[3])
    cmd.signalp_gram_euk <- sprintf("signalp -fasta %s -org euk -format short -prefix %s",
                                    inpath.list[i],
                                    outdir.tmp[4])
    
    system(cmd.signalp_gram_pos)
    system(cmd.signalp_gram_neg)
    system(cmd.signalp_gram_arc)
    system(cmd.signalp_gram_euk)
  }
  
  unlink(dir.tmp,recursive = TRUE)
}


filter_fasta_with_gene_list <- function(accession.list, dir.gene.list, dir.genes.aa="", ex_gene_list, ex_gene=".faa", col=1, header.n=0, sep='\t', min.sep=1, max.sep=1){
  
  if (dir.genes.aa %in% ""){
    dir.genes.aa <- "./data/sequence/genes/aa/"
  }
  
  if (!dir.exists(dir.genes.aa) | (length(list.files(dir.genes.aa)) == 0 )){
    stop("The amino acid fasta file directory does not exist or has no files.")
  }
  
  if (!dir.exists(dir.gene.list) | (length(list.files(dir.gene.list)) == 0 )){
    stop("The gene list directory does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.genes.aa,accession.list,ex_gene)
  gene.list <- paste0(dir.gene.list,accession.list,ex_gene_list)
  outpath.list <- paste0(dir.genes.aa,accession.list,"_filtered",ex_gene)
  tmp.list <- paste0(dir.gene.list,"tmp_list.txt")
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    if (header.n > 0){
      cmd.preprocess <- sprintf("sed '1,%sd;",header.n)
    } else {
      cmd.preprocess <- sprintf("")
    }
    
    if (max.sep == 1){
      cmd.preprocess <- sprintf("%ss/%s/\t/g' %s",
                                cmd.preprocess,
                                sep,
                                gene.list[i])
    } else if (is.infinite(max.sep)) {
      cmd.preprocess <- sprintf("%ss/%s\\{%s,\\}/\t/g' %s",
                                cmd.preprocess,
                                sep,
                                min.sep,
                                gene.list[i])
    } else {
      cmd.preprocess <- sprintf("%ss/%s\\{%s,%s\\}/\t/g' %s",
                                cmd.preprocess,
                                sep,
                                min.sep,
                                max.sep,
                                gene.list[i])
    }
    
    cmd.preprocess <- sprintf("%s | cut -f %s > %s",
                              cmd.preprocess,
                              col,
                              tmp.list)
    
    system(cmd.preprocess)
    
    cmd.extract <- sprintf("seqtk subseq %s %s > %s",
                           inpath.list[i],
                           tmp.list,
                           outpath.list[i])
    
    system(cmd.extract)
  }
  
  unlink(tmp.list)
}