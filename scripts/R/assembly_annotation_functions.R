#create standardized directory architecture for metagenome analysis
create_all_subdirectories <- function(override=FALSE){
  
  if (override == TRUE){
    unlink("data/sequence",recursive = TRUE)
  }
  
  dir.fastq.raw <- "data/sequence/fastq/raw/"
  dir.fastq.trim <- "data/sequence/fastq/trimmed/"
  dir.assembly <- "data/sequence/assemblies/"
  dir.genes.aa <- "data/sequence/genes/aa/"
  dir.genes.nuc <- "data/sequence/genes/nuc/"
  dir.genes.gff <- "data/sequence/genes/gff/"
  dir.annotation.hmm <- "data/sequence/annotation/hmm/"
  dir.database <- "data/database/"
  
  if (!dir.exists(dir.fastq.raw)){
    dir.create(dir.fastq.raw,recursive = TRUE)
  }
  
  if (!dir.exists(dir.fastq.trim)){
    dir.create(dir.fastq.trim,recursive = TRUE)
  }
  
  if (!dir.exists(dir.assembly)){
    dir.create(dir.assembly,recursive = TRUE)
  }
  
  if (!dir.exists(dir.genes.aa)){
    dir.create(dir.genes.aa,recursive = TRUE)
  }
  
  if (!dir.exists(dir.genes.nuc)){
    dir.create(dir.genes.nuc,recursive = TRUE)
  }
  
  if (!dir.exists(dir.genes.gff)){
    dir.create(dir.genes.gff,recursive = TRUE)
  }
  
  if (!dir.exists(dir.annotation.hmm)){
    dir.create(dir.annotation.hmm,recursive = TRUE)
  }
  
  if (!dir.exists(dir.database)){
    dir.create(dir.database,recursive = TRUE)
    warning('Files for constructing databases will need to be manually added to this directory if this is where you want to store them...')
  }
}

#create standardized directory architecture for metagenome analysis
create_some_subdirectories <- function(override=FALSE){
  
  if (override == TRUE){
    unlink("data/sequence",recursive = TRUE)
  }
  
  dir.fastq.raw <- "data/sequence/fastq/raw/"
  dir.fastq.trim <- "data/sequence/fastq/trimmed/"
  dir.assembly <- "data/sequence/assemblies/"
  dir.genes.aa <- "data/sequence/genes/aa/"
  dir.genes.nuc <- "data/sequence/genes/nuc/"
  dir.genes.gff <- "data/sequence/genes/gff/"
  
  if (!dir.exists(dir.fastq.raw)){
    dir.create(dir.fastq.raw,recursive = TRUE)
  }
  
  if (!dir.exists(dir.fastq.trim)){
    dir.create(dir.fastq.trim,recursive = TRUE)
  }
  
  if (!dir.exists(dir.assembly)){
    dir.create(dir.assembly,recursive = TRUE)
  }
  
  if (!dir.exists(dir.genes.aa)){
    dir.create(dir.genes.aa,recursive = TRUE)
  }
  
  if (!dir.exists(dir.genes.nuc)){
    dir.create(dir.genes.nuc,recursive = TRUE)
  }
  
  if (!dir.exists(dir.genes.gff)){
    dir.create(dir.genes.gff,recursive = TRUE)
  }
}

#fastq-dump
#download sequence short read sequences from NCBI's SRA
download_SRR <- function(accession.list, numCores=1){
  
  registerDoParallel(numCores)
  
  dir.fastq.raw <- "./data/sequence/fastq/raw/"
  if (!dir.exists(dir.fastq.raw)){
    stop("The output directory for fastq-dump does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  n.accession <- length(accession.list)
  foreach(i=1:n.accession) %dopar% {
    cmd.download <- sprintf("fastq-dump -O %s --split-files %s",dir.fastq.raw,accession.list[i])
    system(cmd.download)
  }
}

#fastp
#trim and quality control short read sequences
trim_reads <- function(accession.list, fastp.cores=16){
  
  #check to make sure user did exceed fastp's core limit 
  if (fastp.cores > 16){
    fastp.cores <- 6
    warning('Fastp caps the number of threads at 16.\n\n fastp.cores has been set to 16.')
  }
  
  dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
  if (!dir.exists(dir.fastq.trim)){
    stop("The output directory for trimmomatic does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  dir.fastq.raw <- "./data/sequence/fastq/raw/"
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
assemble_reads <- function(accession.list, megahit.cores=30, min.contig.length = 1000){
  
  dir.assembly <- "./data/sequence/assemblies/"
  if (!dir.exists(dir.assembly)){
    stop("The output directory for megahit does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
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
predict_ORFs <- function(accession.list,extension=".contigs.fa"){
  
  dir.genes.aa <- "data/sequence/genes/aa/"
  dir.genes.nuc <- "data/sequence/genes/nuc/"
  dir.genes.gff <- "data/sequence/genes/gff/"
  if (!dir.exists(dir.genes.aa) | !dir.exists(dir.genes.nuc) | !dir.exists(dir.genes.gff)){
    stop("The output directory for genes does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  dir.assembly <- "./data/sequence/assemblies/"
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
  #   stop("The output directory for databases does not exist.\n\nConsider executing the create_sequence_directories function first.")
  # }
  
  warning("This generates the database in the same directory as the file.\n\nIf you do not have write permissions, please considering copying 'file' to a new location.")
  
  cmd.makedb <- switch(fun,
                       'hmmpress' = sprintf("hmmpress -f %s",file),
                       'makeblastdb' = sprintf("makeblastdb -dbtype %s -in %s -out %s",type,file,file))
  
  system(cmd.makedb)
  
}

#annotate genes with hmmscan
annotate_hmmscan <- function(accession.list, db, hmm.cores=30, evalue=1e-10){
  
  dir.annotation.hmm <- "./data/sequence/annotation/hmm/"
  if (!dir.exists(dir.genes.aa)){
    stop("The output directory for hmm annotations does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  dir.genes.aa <- "./data/sequence/genes/aa/"
  if (!dir.exists(dir.genes.aa) | (length(list.files(dir.genes.aa)) == 0 )){
    stop("TThe output directory for amino acid fasta files does not exist or has no files.")
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
    invisible(system(cmd.hmmscan)) #I make the terminal output invisible because it slows the code down too much
  }
  
  
}

