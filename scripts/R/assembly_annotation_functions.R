#create standardized directory architecture for metagenome analysis
create_sequence_directories <- function(override=FALSE){
  
  if (override == TRUE){
    unlink("data/sequence",recursive = TRUE)
  }
  
  dir.fastq.raw <- "data/sequence/fastq/raw/"
  dir.fastq.trim <- "data/sequence/fastq/trimmed/"
  dir.assembly.contigs <- "data/sequence/assemblies/contigs/"
  dir.assembly.aa <- "data/sequence/assemblies/genes/aa/"
  dir.assembly.nuc <- "data/sequence/assemblies/genes/nuc/"
  dir.assembly.gff <- "data/sequence/assemblies/genes/gff/"
  
  if (!dir.exists(dir.fastq.raw)){
    dir.create(dir.fastq.raw,recursive = TRUE)
  }
  
  if (!dir.exists(dir.fastq.trim)){
    dir.create(dir.fastq.trim,recursive = TRUE)
  }
  
  if (!dir.exists(dir.assembly.contigs)){
    dir.create(dir.assembly.contigs,recursive = TRUE)
  }
  
  if (!dir.exists(dir.assembly.aa)){
    dir.create(dir.assembly.aa,recursive = TRUE)
  }
  
  if (!dir.exists(dir.assembly.nuc)){
    dir.create(dir.assembly.nuc,recursive = TRUE)
  }
  
  if (!dir.exists(dir.assembly.gff)){
    dir.create(dir.assembly.gff,recursive = TRUE)
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
    cmd.download <- sprintf("fastp -j %s.json -h %s.html -i %s_1.fastq -I %s_2.fastq -o %s_1.fastq -O %s_2.fastq -w %s",
                            outpath.list[i],
                            outpath.list[i],
                            inpath.list[i],
                            inpath.list[i],
                            outpath.list[i],
                            outpath.list[i],
                            fastp.cores)
    system(cmd.download)
  }
}

#megahit
#assembly high quality reads into contigs
assemble_reads <- function(accession.list, megahit.cores=30, min.contig.length = 1000){
  dir.assembly.contigs <- "./data/sequence/assemblies/contigs/"
  if (!dir.exists(dir.assembly.contigs)){
    stop("The output directory for megahit does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
  if (!dir.exists(dir.fastq.trim) | (length(list.files(dir.fastq.trim)) == 0 )){
    stop("The directory containing raw short read sequences either does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.fastq.trim,accession.list)
  outdirs <- paste0(dir.assembly.contigs,accession.list)
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    cmd.download <- sprintf("megahit -1 %s_1.fastq -2 %s_2.fastq --presets meta-sensitive -t %s -o %s --out-prefix %s --min-contig-len %s",
                            inpath.list[i],
                            inpath.list[i],
                            megahit.cores,
                            outdirs[i],
                            accession.list[i],
                            min.contig.length)
    system(cmd.download)
  }
}

#megahit
#assembly high quality reads into contigs
predict_ORFs <- function(accession.list){
  dir.assembly.contigs <- "./data/sequence/assemblies/contigs/"
  if (!dir.exists(dir.assembly.contigs)){
    stop("The output directory for megahit does not exist.\n\nConsider executing the create_sequence_directories function first.")
  }
  
  dir.fastq.trim <- "./data/sequence/fastq/trimmed/"
  if (!dir.exists(dir.fastq.trim) | (length(list.files(dir.fastq.trim)) == 0 )){
    stop("The directory containing raw short read sequences either does not exist or has no files.")
  }
  
  inpath.list <- paste0(dir.fastq.trim,accession.list)
  outdirs <- paste0(dir.assembly.contigs,accession.list)
  
  n.accession <- length(accession.list)
  for(i in 1:n.accession) {
    cmd.download <- sprintf("megahit -1 %s_1.fastq -2 %s_2.fastq --presets meta-sensitive -t %s -o %s --out-prefix %s --min-contig-len %s",
                            inpath.list[i],
                            inpath.list[i],
                            megahit.cores,
                            outdirs[i],
                            accession.list[i],
                            min.contig.length)
    system(cmd.download)
  }
}
