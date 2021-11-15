library(tidyverse)
library(RColorBrewer)
# library(caret)
# library(vegan)
# library(umap)
# library(kernlab)

source('scripts/stability_functions.R')


# dbcan_dir<-'data/raw/annotation/all.dbcan.tsv'
stability_dir<-'data/M0059E/foldx_GH29/stability_results/347-M0059C_0.72.tsv'
anno.abund_dir<-'data/processed/annotation_frequencies/'
# processed_dir<-'data/processed/processed_stability_data.tsv'
# geochem_dir<-'data/processed/chemistry_for_stability_predictions.tsv'
lookup_dir<-'data/M0059E/dbcan_annotation/blast/GH29/SRR7066493_SRR7066492matches.tsv'


df_stab<-read.table(stability_dir,sep='\t',header=TRUE) # %>%
# filter(depth_m == 0.72) %>%
# select(-depth_m,-event,gene,total_energy)
df_lookup<-read.table(lookup_dir,header=TRUE,sep='\t')    


df_stab$dataset<-grep("SRR",unlist(strsplit(df_stab$gene,"_")),value=TRUE)
df_stab$dataset<-gsub("SRR7066492","15 mbsf",df_stab$dataset)
df_stab$dataset<-gsub("SRR7066493","0.25 mbsf",df_stab$dataset)

df_terms<-df_stab %>%
  pivot_longer(cols=c(-dataset,-gene,-residue_number),
               names_to="term",
               values_to="energy_kcal_mol")

ggplot(df_terms,aes(term,energy_kcal_mol,color=dataset))+geom_boxplot()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))



df_chem<-read.table('data/processed/chemistry_for_stability_predictions.tsv',header=TRUE,sep='\t')  %>%
  smooth_chemistry_profile()

tmp<-df_stab %>% 
  filter(dataset %in% "15 mbsf") %>% 
  select(gene,total_energy) %>% 
  rename(gene1=gene,`15 mbsf`=total_energy) %>%
  inner_join(df_lookup,by="gene1")

tmp<-df_stab %>% 
  filter(dataset %in% "0.25 mbsf") %>% 
  select(gene,total_energy) %>% 
  rename(gene2=gene,`0.25 mbsf`=total_energy) %>%
  inner_join(tmp,by="gene2")

df_homolog_paired<- tmp %>%
  mutate(DDG=`15 mbsf`-`0.25 mbsf`) %>%
  mutate(DDG=case_when(DDG<0~"Negative",
                       TRUE~"Postive")) %>%
  pivot_longer(cols=c(`15 mbsf`,`0.25 mbsf`),names_to="depth",values_to="total_energy") %>%
  group_by(gene1,gene2) %>%
  unite(paired,c("gene1","gene2")) 
# mutate(pairing=row_number())
# rowid_to_column("ID")

# df_homology=tmp %>%
#   select(total_energy1,total_energy2,length,pident,bitscore,evalue) %>%
#   mutate(evalue=log10(evalue)) %>%
#   pivot_longer(cols=c(evalue,pident,bitscore),names_to="metric") %>%
#   mutate(DDG=total_energy1-total_energy2)
# 
# ggplot(df_homology,aes(value,DDG))+
#   geom_point()+
#   geom_smooth(method='lm') +
#   facet_wrap(.~metric,scales="free") +
#   theme_bw()+
#   ylab(expression(paste(Delta,Delta,'G',~degree,' (kcal/mol)'))) +
#   xlab('')


anova(lm(total_energy~depth,df_homolog_paired))
kruskal.test(total_energy~depth,df_homolog_paired)

ggplot(df_homolog_paired,aes(depth,total_energy))+
  geom_boxplot() +
  geom_line(aes(group=as.factor(paired),color=DDG),position = position_dodge(0.2)) +
  geom_point(aes(group=as.factor(paired),color=DDG),position = position_dodge(0.2)) +
  scale_color_brewer(expression(paste(Delta,Delta,'G',~degree,' (kcal/mol)')),palette="Set2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab(expression(paste(Delta,'G',~degree,' (kcal/mol)'))) +
  xlab('')





#compare annotation stability variance against annotation coverage/evalue variance



