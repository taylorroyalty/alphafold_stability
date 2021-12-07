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
  dir.genes.abundance <- "data/sequence/genes/abundance/"
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
  
  if ((!dir.exists(dir.genes.abundance) & create == TRUE) | (!dir.exists(dir.fastq.raw) & base.layer %in% "genes")){
    create <- TRUE
    dir.create(dir.genes.abundance,recursive = TRUE)
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
    dir.assembly <- "./data/sequence/sequences/"
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

    if (!file.exists(paste0(inpath.list[i],'_1.fastq')) | !file.exists(paste0(inpath.list[i],'_2.fastq'))){
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
predict_ORFs <- function(accession.list, dir.assembly="", dir.genes.aa="", dir.genes.nuc="", dir.genes.gff="", in.ex = ".contigs.fa",out.ex="", gff.only=FALSE, extension=".contigs.fa"){
  
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
    dir.assembly <- "./data/sequence/sequences/"
  }
  
  if (!dir.exists(dir.assembly) | (length(list.files(dir.assembly)) == 0 )){
    stop("The directory containing assemblies/genomes either does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.assembly,accession.list,in.ex)
  outdir.aa <- paste0(dir.genes.aa,accession.list,out.ex)
  outdir.nuc <- paste0(dir.genes.nuc,accession.list,out.ex)
  outdir.gff <- paste0(dir.genes.gff,accession.list,out.ex)
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    # inpath.tmp <- paste0(inpath.list[i],'/',accession.list[i],extension)
    if (!file.exists(inpath.list[i])){
      warning(sprintf("The file: '%s' does not exist.",inpath.list[i]))
      next
    }
    
    if (gff.only ==TRUE){
      
      cmd.ORF <- sprintf("prodigal -q -f gff -i %s -o %s.gff",
                         inpath.list[i], #genomes/assemblies
                         outdir.gff[i]) #generic feature format
      
      system(cmd.ORF)
      
    } else {
      
      cmd.ORF <- sprintf("prodigal -q -f gff -i %s -a %s.faa -d %s.fna -o %s.gff",
                         inpath.list[i], #genomes/assemblies
                         outdir.aa[i], #amino acid sequences
                         outdir.nuc[i], #nucleic acid sequences
                         outdir.gff[i]) #generic feature format
      
      system(cmd.ORF)
      
      cmd.clean.aa <- sprintf("sed -i 's/\\*//g;s/^>/>%s\\_/' %s",
                              accession.list[i],
                              paste0(outdir.aa[i],'.faa')) #remove astrieks and replace contig name with accession using sed
      
      cmd.clean.nuc <- sprintf("sed -i 's/^>/>%s\\_/' %s",
                               accession.list[i],
                               paste0(outdir.nuc[i],'.fna')) #replace contig name with accession using sed
      system(cmd.clean.aa)
      system(cmd.clean.nuc)
    }
    

    
    
    
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
annotate_hmmscan <- function(accession.list, dir.genes.aa="", dir.annotation.hmm="", db, top.hit = TRUE, hmm.cores=30, evalue=1e-10){
  
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
  
  annotate_hmmscan(accession.list, dir.genes.aa=dir.genes.aa, dir.annotation.hmm=dir.annotation.hmm, db, hmm.cores=hmm.cores, evalue=evalue)
  
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
    
    #here I filter for the annotation with the lowest evalue--in some cases a gene is annotated multiple times with either different annotations or the same annotation but shifted AA frames
    read.table(outpath.list[i], sep='\t') %>%
      group_by(V3) %>%
      filter(V5 == min(V5)) %>%
      write.table(file=outpath.list[i],sep='\t',quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

#signalp v5.0 
#identifies putative signal peptides and provides the option to cleave at predicted cleavage sites
#Note that signalp is stored in the local /usr/local/bin and is in the usr PATH (for marie). 
filter_signal_peptide <- function(accession.list, dir.genes.aa="", dir.extracellular="", dir.intracellular="", dir.signalp="", faa_ex=".faa", cleave = FALSE){
  
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
  if (dir.signalp %in% ""){
    dir.signalp <- "./data/sequence/genes/aa/signalp/signalp/"
  }
  
  unlink(dir.signalp,recursive = TRUE)
  dir.create(dir.signalp,recursive = TRUE)
  outdir.tmp <- paste0(dir.signalp,accession.list) #temporary files for each organism in sigalp
  
  inpath.list <- paste0(dir.genes.aa,accession.list,faa_ex)
  outpath_ex.list <- paste0(dir.extracellular,accession.list,"_extracellular",faa_ex)
  outpath_in.list <- paste0(dir.intracellular,accession.list,"_intracellular",faa_ex)
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    outpath.list <- paste0(outdir.tmp[i],c("_gram_pos","_gram_neg","_arc","_euk"))
    
    cmd.signalp_gram_pos <- sprintf("signalp -fasta %s -org gram+ -format short -prefix %s",
                                    inpath.list[i],
                                    outpath.list[1])
    
    cmd.signalp_gram_neg <- sprintf("signalp -fasta %s -org gram- -format short -prefix %s",
                                    inpath.list[i],
                                    outpath.list[2])
    
    cmd.signalp_gram_arc <- sprintf("signalp -fasta %s -org arch -format short -prefix %s",
                                    inpath.list[i],
                                    outpath.list[3])
    
    cmd.signalp_gram_euk <- sprintf("signalp -fasta %s -org euk -format short -prefix %s",
                                    inpath.list[i],
                                    outpath.list[4])
    
    cmd.signalp_all <- c(cmd.signalp_gram_pos, cmd.signalp_gram_neg, cmd.signalp_gram_arc, cmd.signalp_gram_euk)
    
    result.ex <- data.frame(NULL)
    for (j in 1:4){
      system(cmd.signalp_all[j])
      result.path.tmp <- paste0(outpath.list[j],"_summary.signalp5")
      result.ex <- read.table(result.path.tmp,skip=2, sep = '\t') %>%
        filter(!V2 %in% "OTHER") %>%
        select(-V2) %>%
        group_by(V1) %>%
        pivot_longer(cols=-c(names(.)[1],tail(names(.),1)),names_to = "type",values_to = "prob") %>%
        filter(prob == max(prob)) %>%
        distinct() %>%
        select(-type) %>%
        rename(gene=V1,cleavage=names(.)[2]) %>%
        ungroup() %>%
        rbind(result.ex)
      
    }
    
    #partition into extracellular and intracellular gene lists
    result.ex<-result.ex %>%
      group_by(gene) %>%
      filter(prob == max(prob)) %>%
      distinct() %>%
      mutate(cleavage=as.numeric(substr(cleavage,9,10)))
    
    result.ex$cleavage[is.na(result.ex$cleavage)]<-0
    
    result.in <- read.table(result.path.tmp,skip=2, sep = '\t') %>%
      select(V1) %>%
      filter(!V1 %in% result.ex$gene) %>%
      distinct() %>%
      cbind(cleavage=0)
    
    #overwrite faa files if they exists
    system(sprintf("> %s",outpath_ex.list[i]))
    system(sprintf("> %s",outpath_in.list[i]))
    
    #build command for isolating seq--if cleave true, use sed to trucate sequences 
    cmd.cleave <- c("grep -A 1 \"%s \" %s | sed '2s/^.\\{%s\\}//' >> %s")
    if (cleave == FALSE) {
      result.ex$cleavage <- result.ex$cleavage*0  
    }
    
    for (seq in 1:nrow(result.ex)){
      cmd.isolate.ex <- sprintf(cmd.cleave,
                                result.ex$gene[seq],
                                inpath.list[i],
                                result.ex$cleavage[seq],
                                outpath_ex.list[i])
      
      system(cmd.isolate.ex)
    }
    
    for (seq in 1:nrow(result.in)){
      cmd.isolate.in <- sprintf(cmd.cleave,
                                result.in$V1[seq],
                                inpath.list[i],
                                result.in$cleavage[seq],
                                outpath_in.list[i])
      
      system(cmd.isolate.in)
    }
    
  }
  # 
  # unlink(dir.tmp,recursive = TRUE)
}

#filters fasta files based on gene lists in a target file. gene lists are to have one gene id per line in target file
filter_fasta_with_gene_list <- function(accession.list, dir.gene.list, ex_gene_list, dir.genes="", out_ex="_filtered.faa", ex_gene=".faa", col=1, header.n=0, sep='\t', min.sep=1, max.sep=1, apply.grep=FALSE, grep.filter=""){
  
  if (dir.genes %in% ""){
    dir.genes <- "./data/sequence/genes/aa/"
  }
  
  if (!dir.exists(dir.genes) | (length(list.files(dir.genes)) == 0 )){
    stop("The fasta file directory does not exist or has no files.")
  }
  
  if (!dir.exists(dir.gene.list) | (length(list.files(dir.gene.list)) == 0 )){
    stop("The gene list directory does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.genes,accession.list,ex_gene)
  gene.list <- paste0(dir.gene.list,accession.list,ex_gene_list)
  outpath.list <- paste0(dir.genes,accession.list,out_ex)
  tmp.list <- paste0(dir.gene.list,"tmp_list.txt")
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    
    if (header.n > 0){
      cmd.preprocess <- sprintf("1,%sd;",header.n)
    } else {
      cmd.preprocess <- sprintf("")
    }
    
    if (max.sep == 1){
      cmd.preprocess <- sprintf("sed '%ss/%s/\t/g' %s",
                                cmd.preprocess,
                                sep,
                                gene.list[i])
    } else if (is.infinite(max.sep)) {
      cmd.preprocess <- sprintf("sed '%ss/%s\\{%s,\\}/\t/g' %s",
                                cmd.preprocess,
                                sep,
                                min.sep,
                                gene.list[i])
    } else {
      cmd.preprocess <- sprintf("sed '%ss/%s\\{%s,%s\\}/\t/g' %s",
                                cmd.preprocess,
                                sep,
                                min.sep,
                                max.sep,
                                gene.list[i])
    }
    
    
    if (apply.grep == TRUE){
      cmd.preprocess <- sprintf("%s | grep \"%s\"",
                                cmd.preprocess,
                                grep.filter)
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
  
}

#bowtie2; samtools
#align short-reads to target genes
align_short_reads <- function(accession.list, dir.genes.abundance="", dir.genes.nuc="", dir.fastq.trim="", ex_gene=".fna", align.cpu=30){
  
  if (dir.genes.nuc %in% ""){
    dir.genes.nuc <- "data/sequence/genes/nuc/"
  }
  
  if (dir.fastq.trim %in% ""){
    dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
  }
  
  if (dir.genes.abundance %in% ""){
    dir.genes.abundance <- "./data/sequence/genes/abundance/"
  }
  
  dir.index <- paste0(dir.genes.abundance,'/alignment/')
  dir.results <- paste0(dir.genes.abundance,'/results/')
  
  dir.create(dir.index,recursive = TRUE)
  dir.create(dir.results,recursive = TRUE)
  
  inpath.fastq.list <- paste0(dir.fastq.trim,accession.list)
  inpath.gene.nuc.list <- paste0(dir.genes.nuc,accession.list,ex_gene)
  outpath.results.list <- paste0(dir.results,accession.list)
  outpath.index.list <- paste0(dir.index,accession.list)
  outpath.sam.list <- paste0(dir.index,accession.list,".sam")
  outpath.bam.list <- paste0(dir.index,accession.list,".bam")
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    cmd.bw2.build <- sprintf("bowtie2-build -f %s %s --threads %s",
                             inpath.gene.nuc.list[i],
                             outpath.index.list[i],
                             align.cpu)
    
    system(cmd.bw2.build)
    
    cmd.bw2 <- sprintf("bowtie2 -x %s -1 %s_1.fastq -2 %s_2.fastq -S %s -p %s",
                       outpath.index.list[i],
                       inpath.fastq.list[i],
                       inpath.fastq.list[i],
                       outpath.sam.list[i],
                       align.cpu)
    
    system(cmd.bw2)
    
    
    cmd.bam <- sprintf("samtools view -@ 30 -bS %s | samtools sort -@ 30 -o %s",
                       outpath.sam.list[i],
                       outpath.bam.list[i])
    
    system(cmd.bam)
  }
  
  
}


#featureCounts
#convert bam to abundance files--here, I utilize gff provided by prodigal and convert those gff into gtf using PROKKA's gff to gtf conversion bash commands
alignment_to_abundance <- function(accession.list, dir.genes.abundance="", ex.abund="_abundance.tsv", dir.align="", ex.align=".bam", ex.gtf=".gtf", dir.genes.gtf="", gff2gtf=FALSE, ex.gff=".gff", dir.genes.gff="", feat.cores=30) {
  
  if (dir.genes.abundance %in% ""){
    dir.genes.abundance <- "./data/sequence/genes/abundance/"
  }
  
  if (dir.align %in% ""){
    dir.align <- "./data/sequence/genes/abundance/alignment/"
  }
  
  if (dir.genes.gtf %in% ""){
    dir.genes.gtf <- "./data/sequence/genes/gtf/"
  }
  
  n.accession <- length(accession.list)
  
  if (gff2gtf == TRUE){
    
    if (dir.genes.gff %in% ""){
      dir.genes.gff <- "./data/sequence/genes/gff/"
    }
    
    if (!dir.exists(dir.genes.gtf)){
      dir.create(dir.genes.gtf,recursive = TRUE)
    }
    
    gff.list <- paste0(dir.genes.gff,accession.list,ex.gff)
    gtf.list <- paste0(dir.genes.gtf,accession.list,ex.gtf)
    
    for (i in 1:n.accession){
      
      cmd.gff2gtf <- sprintf("grep -v \"#\" %s | grep \"ID=\" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  awk -v OFS='\t' '{print $1,\"PROKKA\",\"CDS\",$2,$3,\".\",$4,\".\",\"gene_id \" $5}' > %s",
                             gff.list[i],
                             gtf.list[i])
      
      system(cmd.gff2gtf)
      
    }
  }
  
  align.list <- paste0(dir.align,accession.list,ex.align)
  outpath.list <- paste0(dir.genes.abundance,accession.list,ex.abund)
  
  cmd.feature <- sprintf("featureCounts -p -T %s -s 0 -t CDS -a %s -o %s %s",
                         feat.cores,
                         gtf.list[i],
                         outpath.list[i],
                         align.list[i])
  
  system(cmd.feature)
}


#generate null sequences using swiss prot bulk AA frequencies
null_sequence <- function(filepath,seq.len,rep=1,w=c(),alphabet=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")){
  
  if (!dir.exists(dirname(filepath))){ 
    dir.create(dirname(filepath),recursive = TRUE)
  }
  
  if (is.null(w)){
    #weights correspond to AA frequencies reported in swiss prot 2021_03 https://web.expasy.org/docs/relnotes/relstat.html
    w=c(0.0825,0.0553,0.0406,0.0546,0.0138,0.0393,0.0672,0.0707,0.0227,0.0591,0.0965,0.0580,0.0241,0.0386,0.0473,0.0664,0.0535,0.0110,0.0292,0.0686)
    w=w/sum(w) #this deals with precision errors from swiss prot / small presence of Alx, Glx, Xaa
  } 
  
  if (!sum(w)==1){
    stop("The sum of the weight vector needs to equal 1.\n\nIf precision is an issue try: w=w/sum(w)")
  }
  
  if (!length(w)==length(alphabet)){
    stop("The weight vector is not the same length as the alphabet")
  }
  
  system(sprintf("> %s",filepath))
  
  null_sequence.subfun <- function(index,filepath,seq.len,w,alphabet){
    
    header<-paste("\\>sequence_",as.character(seq.len),"_",as.character(index),sep="")
    system(sprintf("echo %s >> %s",header,filepath))
    
    system(sprintf("echo %s >> %s",paste(sample(alphabet,seq.len,replace = TRUE,prob=w),collapse = ""),filepath))
  }
  
  rep.vec=1:rep
  
  for (l in seq.len){
    sapply(rep.vec,null_sequence.subfun,filepath=filepath,seq.len=l,w=w,alphabet=alphabet)
  }
}

#blastp
#query a custom database using blastp
blastp_query <- function(accession.list, db.path, dir.output="", dir.genes="", genes_ex='.faa', blast.options="", blast.cores=30){
  
  warning("Note that this function uses parallelization native to blastp's implementation.\n\nThis is slower than gnu parallel using single threaded blastp on large blast searches.\n\nMy point is that you should use gnu parallel via the shell for large blast searchs...")
  
  if (dir.output %in% ""){
    dir.output <- paste("./data/sequence/blast_results/", basename(db.path), collapse = "")
  }
  
  if (dir.genes %in% ""){
    dir.genes <- "./data/sequence/genes/aa/"
  }
  
  if (!dir.exists(dir.output)){
    dir.create(dir.output,recursive = TRUE)
  }
  
  inpath.list <- paste0(dir.genes,accession.list,genes_ex)
  outpath.list <- paste0(dir.output,accession.list,'.blastp')
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    cmd.blastp <- sprintf("blastp -query %s -out %s -db %s -num_threads %s %s",
                          inpath.list[i],
                          outpath.list[i],
                          db.path,
                          blast.cores,
                          blast.options)
    
    system(cmd.blastp)
  }
  
}

#seqtk sample -- note that this command samples without replacement!!!! 
#randomly sample fasta file
sample_fasta <- function(accession.list, n, gene_ex, dir.gene){
  
  inpath.list <- paste0(dir.gene,accession.list,gene_ex)
  outpath.list <- paste0(dir.gene,accession.list,"_",as.character(n),gene_ex)
  
  n.accession <- length(accession.list)
  for (i in 1:n.accession){
    seed.int <- sample(1:1e4,1)
    cmd.sample <- sprintf("seqtk sample -s %s %s %s > %s",
                          seed.int,
                          inpath.list[i],
                          n,
                          outpath.list[i])
    system(cmd.sample)
  }
}


