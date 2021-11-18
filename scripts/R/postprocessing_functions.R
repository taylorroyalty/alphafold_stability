###Load transcript data and tabulate counts mapped per read
load_transcript_count_data<-function(directory){
  #sub-function for loading data
  load_data<-function(x){
    tmp_df=read.table(x,header=TRUE,sep='\t')
    tmp_df$dataset=strsplit(tail(strsplit(x,'/')[[1]],n=1),'[.]')[[1]][1]
    
    return(tmp_df)
  }
  
  vec<-list.files(directory,full.names = TRUE)  
  
  df<-map_dfr(vec,load_data) %>%
    group_by(gene) %>%
    summarize(abundance=sum(counts)) %>%
    mutate(transcribed=case_when(abundance==0~FALSE,
                                 TRUE~TRUE)) %>%
    ungroup() %>%
    select(gene,transcribed)
  
  df
  
}


###Load stability data
load_stability_data<-function(directory){
  #sub-function for loading data
  load_data<-function(x){
    tmp_df=read.table(x,header=TRUE,sep='\t')
    tmp_df$dataset=strsplit(tail(strsplit(x,'/')[[1]],n=1),'[.]')[[1]][1]
    
    return(tmp_df)
  }
  
  vec<-list.files(directory,full.names = TRUE)  
  
  df<-map_dfr(vec,load_data) %>%
    select(gene,total_energy,dataset)
  
  df
}

###Load signalp data
load_signalp_predictions<-function(directory){
  #sub-function for loading signalp predictions
  load_data<-function(x){
    tmp_df<-read.table(x,header=TRUE,sep='\t')
    
    return(tmp_df)
  }
  
  vec<-list.files(directory,full.names = TRUE)  
  
  df<-map_dfr(vec,load_data) %>%
    select(ID,Prediction) %>%
    mutate(Prediction=case_when(Prediction=="OTHER"~0,
                                TRUE~1 )) %>%
    group_by(ID) %>% 
    summarize(signal.peptide=sum(Prediction)>0) %>%
    rename("gene"="ID")
  
  
  df
  
  
}

###Rarefy transcript data
rarefy_transcripts<-function(df,counts='min_dataset',reps=3){
  
  #sub-function performing rarefaction
  rarefy<-function(replicate,df){
    df %>%
      pivot_wider(names_from = gene,values_from=abundance,values_fill=0) %>% #convert into wide format to use vegan::rrafey
      column_to_rownames(var="dataset") %>%
      vegan::rrarefy(counts) %>% #vegan::rrarefey
      as.data.frame() %>%
      rownames_to_column(var="dataset") %>%
      pivot_longer(cols=-dataset, 
                   names_to = "gene",
                   values_to = "abundance") %>%
      cbind(replicate=replicate) #append column that indicates replicate
  }
  
  
  #determine if sampling effort is default
  if (counts == 'min_dataset'){
    counts<-df %>%
      group_by(dataset) %>%
      summarize(abundance=sum(abundance)) %>%
      ungroup() %>%
      slice(which.min(abundance)) %>%
      select(abundance) %>%
      as.numeric()
  }
  
  df_o<-map_dfr(1:reps,rarefy,df=df) %>%
    group_by(dataset,gene) %>%
    summarize(abundance.m=mean(abundance),
              abundance.sd=sd(abundance)) %>%
    filter(abundance.m > 0)
  
  df_o
}

###Downsample to 1 meter stability profiles--acounts for redundant depth intervals in chemistry data for different cores
smooth_energy_profile<-function(df) {
  
  smooth_energy_sub<-function(df_s,mean.vec,index.vec,interval){
    max.value<-df_s %>%
      select(index.vec) %>%
      max() %>%
      ceiling()
    
    min.value<-df_s %>%
      select(index.vec) %>%
      min() %>%
      floor()
    
    df_o<-data.frame(NULL)
    for (i in seq(min.value,max.value,by=interval)){
      
      df_o<-df_s %>%
        filter(!!sym(index.vec) >=i & !!sym(index.vec)<i+interval) %>%
        summarize(index.vec=i+interval/2,
                  mean.vec=mean(!!sym(mean.vec))) %>%
        rbind(df_o)
      
    }
    colnames(df_o)<-c(index.vec,mean.vec)
    
    df_o<-df_o %>%
      na.omit()
    
    df_o
  }
  
  df<-df %>%
    select(-event) %>%
    group_by(gene) %>%
    nest() %>%
    mutate(results=map(data,smooth_energy_sub,mean.vec="total_energy",index.vec="depth_m",interval=1)) %>%
    select(-data) %>%
    unnest(results)
  
  df  
  
  
}

###Downsample to 1 meter chemistry profiles--acounts for redundant depth intervals in chemistry data for different cores
smooth_chemistry_profile<-function(df) {
  
  smooth_chemistry_sub<-function(df_s,mean.vec,index.vec,interval){
    max.value<-df_s %>%
      select(index.vec) %>%
      max() %>%
      ceiling()
    
    min.value<-df_s %>%
      select(index.vec) %>%
      min() %>%
      floor()
    
    df_o<-data.frame(NULL)
    for (i in seq(min.value,max.value,by=interval)){
      
      df_o<-df_s %>%
        filter(!!sym(index.vec) >=i & !!sym(index.vec)<i+interval) %>%
        summarize(index.vec=i+interval/2,
                  mean.vec=mean(!!sym(mean.vec))) %>%
        rbind(df_o)
      
    }
    colnames(df_o)<-c(index.vec,mean.vec)
    
    df_o<-df_o %>%
      na.omit()
    
    df_o
  }
  
  df<-df %>%
    select(-event) %>%
    pivot_longer(cols = c(pH,temperature_k,ionic_strength),names_to ="variable",values_to="value") %>%
    group_by(variable) %>%
    nest() %>%
    mutate(results=map(data,smooth_chemistry_sub,mean.vec="value",index.vec="depth_m",interval=1)) %>%
    select(-data) %>%
    unnest(results) %>%
    pivot_wider(names_from=variable,values_from=value)
  
  df  
  
  
}


delta_delta_G <- function(value.vec,factor.vec,it=1000){
  
  split.vec<-split(value.vec,f=factor.vec)
  
  n=length(value.vec)
  n.half<-n/2
  
  
  delta_delta_G<-as.vector(outer(unlist(split.vec[1]),unlist(split.vec[2]),FUN='-'))
  df<-data.frame(DDG=delta_delta_G,it=NaN,type='Observed')
  
  for (i in 1:it){
    indx<-sample(1:n,n.half)
    tmp1.vec<-value.vec[indx]
    tmp2.vec<-value.vec[-indx]
    
    delta_delta_G.tmp<-as.vector(outer(tmp1.vec,tmp2.vec,FUN='-'))
    
    df<-data.frame(DDG=delta_delta_G.tmp,it=it,type="Permutations") %>%
      rbind(df)
  }
  
  df
  
}

###Fit polynomial function to stability profiles.
polynomial_fit<-function(df,power.vec){
  
  polynomial_fit_sub<-function(df_s,power.vec){
    
    #dynamically generate different powers based on power.vec input
    for (i in power.vec){
      df_s<-data.frame(df_s$depth_m^i) %>%
        cbind(df_s)
      colnames(df_s)[1]<-paste("X",as.character(i),sep="")
    }
    
    #select polynomials for regression
    train.data<- df_s %>%
      select(c((paste("X",power.vec,sep="")),"total_energy")) %>%
      as.matrix()
    
    # Set training control
    train_control <- trainControl(method = "repeatedcv",
                                  number = 5,
                                  repeats = 5,
                                  search = "random")
    
    # Train the model
    elastic_net_model <- train(total_energy ~ .,
                               data = train.data,
                               method = "glmnet",
                               preProcess = c("center", "scale"),
                               tuneLength = 25,
                               trControl = train_control)
    
    y_hat_enet <- predict(elastic_net_model, train.data)
    r2 <- cor(as.data.frame(train.data)$total_energy, y_hat_enet)^2
    
    df_o<-coef(elastic_net_model$finalModel, elastic_net_model$bestTune$lambda) %>%
      as.matrix() %>% 
      t() %>% 
      as.data.frame() %>% 
      rename("Intercept"="(Intercept)") %>%
      cbind(data.frame(r2=r2,lambda=elastic_net_model$bestTune$lambda,alpha=elastic_net_model$bestTune$alpha))
    
    rownames(df_o)<-c()
    
    df_o
    
  }
  
  
  df %>%
    group_by(gene) %>%
    nest() %>%
    mutate(poly_fit=map(data,polynomial_fit_sub,power.vec=power.vec)) %>%
    select(-data) %>%
    unnest(poly_fit)
}


###merge list of data.frames based on list of target columns
inner_join_all <- function(df.list,by){
  df_o<-inner_join(df.list[[1]],df.list[[2]],by=by)
  
  for (i in 3:length(df.list)){
    df_o<-inner_join(df_o,df.list[[i]],by=by)
  }
  
  df_o
}

###Evaluate the distribution of gibbs free energy as a function of depth and annotation, I also evaluate transcribed only and those with signal peptide
energy_annotation_analysis <- function(df,depth.list,min.n=5){
  
  f1 <- df %>%
    filter(depth_m %in% depth.list) %>%
    group_by(annotation) %>%
    mutate(n=n()) %>%
    filter(n >= min.n) %>%
    ggplot() +
    geom_boxplot(aes(annotation,total_energy))+
    ylab(expression(paste(Delta,'G',~degree,' (kcal/mol)'))) +
    xlab("dbCAN Annotation") +
    theme_bw() +
    facet_wrap(.~depth_m)
  
  f2 <- df %>%
    filter(depth_m %in% depth.list,
           transcribed == TRUE) %>%
    group_by(annotation) %>%
    mutate(n=n()) %>%
    filter(n >= min.n) %>%
    ggplot() +
    geom_boxplot(aes(annotation,total_energy))+
    ylab(expression(paste(Delta,'G',~degree,' (kcal/mol) (Transcribed)'))) +
    xlab("dbCAN Annotation") +
    theme_bw() +
    facet_wrap(.~depth_m)
  
  f3 <- df %>%
    filter(depth_m %in% depth.list,
           transcribed == TRUE,
           signal.peptide ==TRUE) %>%
    group_by(annotation) %>%
    mutate(n=n()) %>%
    filter(n >= min.n) %>%
    ggplot() +
    geom_boxplot(aes(annotation,total_energy))+
    ylab(expression(paste(Delta,'G',~degree,' (kcal/mol) (Transcribed + Signal Peptide)'))) +
    xlab("dbCAN Annotation") +
    theme_bw() +
    facet_wrap(.~depth_m)
  
  print(f1)
  print(f2)
  print(f3)
  
  
}


#cluster and dimensionally reduce polynomial fit coefficients--here I variance normalize the coefficients
cluster_profiles <- function(df){
  df <- df %>%
    ungroup() %>%
    select(c("genome","annotation","Intercept"),grep("X",colnames(df))) %>%
    distinct()
  
  pca_o <- df %>%
    select("Intercept",grep("X",colnames(df))) %>%
    scale() %>%
    prcomp() 
  
  df_wide <- df %>%
    mutate(pc1=pca_o$x[,1]) %>%
    group_by(genome,annotation) %>%
    summarize(pc1=mean(pc1)) %>%
    ungroup() %>%
    select(genome,annotation,pc1) %>%
    pivot_wider(names_from = annotation,values_from=pc1, values_fill = 0) %>%
    column_to_rownames(var="genome") %>%
    as.matrix()
  
  # dist.matrix <- vegan::vegdist(as.matrix(df)) %>%
  # as.matrix()
  
  labels <- rownames(df_wide) 
  
  # df_wide <- df %>%
  #   ungroup() %>%
  #   select("Intercept",grep("X",colnames(df))) %>%
  #   # group_by(annotation,genome) %>%
  #   scale() %>%
  #   as.matrix()
  
  within.vec=c()
  n.clust<-seq(2,19)
  for (i in n.clust){
    sc.clusters <- kernlab::specc(df_wide,centers=i)
    within.vec<-c(within.vec,sum(kernlab::withinss(sc.clusters)))
  }
  
  indx <- which(within.vec==min(within.vec))
  
  sc.clusters <- kernlab::specc(df_wide %>% as.matrix(),center=)
  
  umap_o <- umap::umap(df_wide,preserve.seed = FALSE)
  
  df_p <- data.frame(f1=umap_o$layout[,1],
                     f2=umap_o$layout[,2],
                     cluster=as.factor(sc.clusters),
                     labels)
  ggplot(df_p,aes(f1,f2,color=cluster))+geom_point()    
}


plot_rank_abundance <- function(accession.list,dir,ex,sep='\t',col=1,header=TRUE,max.rank=20,name.split=" ",parse.index=1){
  
  load_data<-function(x,sep,header,name.split,parse.index=1){
    vec=read.table(x,header=header,sep=sep)[,col]
    dataset=strsplit(basename(x),name.split)[[1]][parse.index]
    
    df=data.frame(vec,dataset)
    
    return(df)
  }
  
  vec <- paste0(dir,accession.list,ex)
  
  df<-map_dfr(vec,load_data,sep=sep,header=header,name.split=name.split,parse.index=parse.index) %>%
    group_by_all() %>%
    summarize(abundance=n()) %>%
    ungroup() %>%
    group_by(dataset) %>%
    mutate(RA=abundance/sum(abundance),
           rank=rank(-RA,ties.method="first")) %>%
    filter(rank<=max.rank)
  
  ggplot(df,aes(rank,RA,fill=dataset))+
    geom_col(position="dodge") +
    geom_text(aes(label=vec),angle=90,hjust=-0.1,
              size=3,position=position_dodge(0.9)) +
    scale_fill_brewer(palette="Dark2") +
    ylim(0,0.2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab('Rank') +
    ylab('Annotation Relative Abundance')

}

