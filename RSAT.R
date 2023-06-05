# detailed analysis in 11_autoantigen.R

library(seqinr)
library(dplyr)
library(tidyverse)
library(data.table)
library(Biostrings)
library(foreach)
library(doParallel)
library(data.table)
library(grid)
library(gridExtra)
library(effsize)
library(ggpubr)
library(ggplot2)
library(caret)

plot_cumulative_dist_immunogenicity <- function(data, feature = 'LD1', title = NA, remove_legends = F){
  
  data <- data  %>% dplyr::select(Peptide, Immunogenicity, !!sym(feature)) %>% filter(!is.na(!!sym(feature))) %>%
    mutate_if(is.factor, as.character) %>% as.data.frame()  %>% 
    filter(!is.infinite(!!sym(feature)))
  
  d <- data %>% pull(!!sym(feature))
  f <- data %>% pull(Immunogenicity)
  
  if(length(unique(f)) == 0){
    p <- NULL
  }else{
    if(unique(f) %in% c('Positive', 'Negative')){
      data <- data %>%   mutate(Immunogenicity = factor(Immunogenicity, levels = c('Positive', 'Negative')))
      
      npeptides <- data %>% plyr::count(vars = 'Immunogenicity')
      npos <- npeptides$freq[which(npeptides$Immunogenicity == 'Positive')]
      nneg <- npeptides$freq[which(npeptides$Immunogenicity == 'Negative')]
      count_peptides <- paste0('# Pos = ', npos, ' # Neg = ', nneg)
      
    }else{
      npeptides <- data %>% plyr::count(vars = 'Immunogenicity')
      count_peptides <- c()
      for(row in 1:nrow(npeptides)){
        peptide <- paste0(npeptides$Immunogenicity[row], ': ', npeptides$freq[row])
        count_peptides <- c(count_peptides, peptide)
      }; count_peptides <- paste(count_peptides, collapse = ', ')
      
    }
    
    if(length(unique(f)) == 2){
      result <- cohen.d(d ~ f)
      # t_test <- t.test(d~f)
      t_test <- tryCatch( t.test(d~f), error=function(err) NA)
      
      if(!is.na(t_test)){
        if(t_test$p.value < 0.01){
          t_test_p_values <- formatC(t_test$p.value, format = "e", digits = 2)
        }else{
          t_test_p_values <- round(t_test$p.value, 3)
        }
      }else{
        t_test_p_values <- NA
      }
      
      
      p <- data %>% ggplot( aes_string(feature, colour = 'Immunogenicity')) + stat_ecdf(size=2) + theme_bw() +
        scale_color_brewer(palette="Dark2")  +ylab('Cumulative % of peptides') +
        labs(subtitle = paste0('Cohen d = ', round(result$estimate, 3), ', t-test p = ', t_test_p_values,
                               '\n', count_peptides), 
             title = title)  + 
        theme(plot.subtitle = element_text(hjust = 0.5, size = 10.5), plot.title = element_text(hjust = 0.5))
      
    }else{
      
      p <- data %>% ggplot( aes_string(feature, colour = 'Immunogenicity')) + stat_ecdf(size=2) + theme_bw() +
        scale_color_brewer(palette="Dark2")  +ylab('Cumulative % of peptides') +
        labs(subtitle = count_peptides,   title = title)  + 
        theme(plot.subtitle = element_text(hjust = 0.5, size = 10.5), plot.title = element_text(hjust = 0.5))
    }
    
    if(remove_legends == T){
      p <- p + theme(legend.position = "none") 
    }
  }
  return(p)
  
}
plot_violin_immunogenicity <- function(data, feature = 'LD1', title = NA, remove_legends = F){
  
  data <- data  %>% dplyr::select(Peptide, Immunogenicity, !!sym(feature)) %>% filter(!is.na(!!sym(feature))) %>%
    mutate_if(is.factor, as.character) %>% as.data.frame() %>% 
    filter(!is.infinite(!!sym(feature)))
  
  d <- data %>% pull(!!sym(feature))
  f <- data %>% pull(Immunogenicity)
  
  if(length(unique(f)) == 0){
    p <- NULL
  }else{
    if(sort(unique(f)) == c('Negative', 'Positive')){
      data <- data %>%   mutate(Immunogenicity = factor(Immunogenicity, levels = c('Positive', 'Negative')))
      
      npeptides <- data %>% plyr::count(vars = 'Immunogenicity')
      npos <- npeptides$freq[which(npeptides$Immunogenicity == 'Positive')]
      nneg <- npeptides$freq[which(npeptides$Immunogenicity == 'Negative')]
      count_peptides <- paste0('# Pos = ', npos, ' # Neg = ', nneg)
      
    }else{
      npeptides <- data %>% plyr::count(vars = 'Immunogenicity')
      count_peptides <- c()
      for(row in 1:nrow(npeptides)){
        peptide <- paste0(npeptides$Immunogenicity[row], ': ', npeptides$freq[row])
        count_peptides <- c(count_peptides, peptide)
      }; count_peptides <- paste(count_peptides, collapse = ', ')
      
    }
    
    max_y <- max(d, rm.na = T)
    if(length(unique(f)) == 2){
      
      result <- cohen.d(d ~ f)
      t_test <- t.test(d~f)
      if(t_test$p.value < 0.01){
        t_test_p_values <- formatC(t_test$p.value, format = "e", digits = 2)
      }else{
        t_test_p_values <- round(t_test$p.value, 3)
      }
      
      my_comparisons <- list(c('Positive', 'Negative'))
      p <- data %>% ggviolin(x = 'Immunogenicity', y = feature, fill = "Immunogenicity", 
                             palette = c("#006400",  "#FC4E07"),  legend = 'none',
                             add = "boxplot", add.params = list(fill = "white")) +
        stat_compare_means(comparisons = my_comparisons, label = "p.signif") + theme_bw() +
        labs(subtitle = paste0('Cohen d = ', round(result$estimate, 3), ', t-test p = ', t_test_p_values,
                               '\n', count_peptides), 
             title = title)  + 
        theme(plot.subtitle = element_text(hjust = 0.5, size = 10.5), plot.title = element_text(hjust = 0.5)) +
        ylim(NA, max_y*1.25)
      
    }else{
      
      p <- data %>% ggviolin(x = 'Immunogenicity', y = feature, fill = "Immunogenicity", 
                             palette = c("#006400",  "#FC4E07"),  legend = 'none',
                             add = "boxplot", add.params = list(fill = "white")) +
        stat_compare_means(comparisons = my_comparisons, label = "p.signif") + theme_bw() +
        labs(subtitle = count_peptides,   title = title)  + 
        theme(plot.subtitle = element_text(hjust = 0.5, size = 10.5), plot.title = element_text(hjust = 0.5)) +
        ylim(NA, max_y*1.25)
    }
    
    if(remove_legends == T){
      p <- p + theme(legend.position = "none") 
    }
    
  }
  
  return(p)
  
}
plot_boxplot_immunogenicity <- function(data, feature = 'LD1', title = NA, remove_legends = F){
  
  data <- data  %>% dplyr::select(Peptide, Immunogenicity, !!sym(feature)) %>% filter(!is.na(!!sym(feature))) %>%
    mutate_if(is.factor, as.character) %>% as.data.frame()  %>% 
    filter(!is.infinite(!!sym(feature)))
  
  d <- data %>% pull(!!sym(feature))
  f <- data %>% pull(Immunogenicity)
  
  if(length(unique(f)) == 0){
    p <- NULL
  }else{
    if(unique(f) %in% c('Positive', 'Negative')){
      data <- data %>%   mutate(Immunogenicity = factor(Immunogenicity, levels = c('Positive', 'Negative')))
      
      npeptides <- data %>% plyr::count(vars = 'Immunogenicity')
      npos <- npeptides$freq[which(npeptides$Immunogenicity == 'Positive')]
      nneg <- npeptides$freq[which(npeptides$Immunogenicity == 'Negative')]
      count_peptides <- paste0('# Pos = ', npos, ' # Neg = ', nneg)
      
    }else{
      npeptides <- data %>% plyr::count(vars = 'Immunogenicity')
      count_peptides <- c()
      for(row in 1:nrow(npeptides)){
        peptide <- paste0(npeptides$Immunogenicity[row], ': ', npeptides$freq[row])
        count_peptides <- c(count_peptides, peptide)
      }; count_peptides <- paste(count_peptides, collapse = ', ')
      
    }
    
    max_y <- max(d, rm.na = T)
    if(length(unique(f)) == 2){
      result <- cohen.d(d ~ f)
      t_test <- t.test(d~f)
      if(t_test$p.value < 0.01){
        t_test_p_values <- formatC(t_test$p.value, format = "e", digits = 2)
      }else{
        t_test_p_values <- round(t_test$p.value, 3)
      }
      
      my_comparisons <- list(c('Positive', 'Negative'))
      p <- data %>%  ggboxplot(x = 'Immunogenicity', y = feature, fill = "Immunogenicity", 
                               palette = c("#006400",  "#FC4E07"),  legend = 'none') +
        stat_compare_means(comparisons = my_comparisons, label = "p.signif") + theme_bw() +
        labs(subtitle = paste0('Cohen d = ', round(result$estimate, 3), ', t-test p = ', t_test_p_values,
                               '\n', count_peptides), 
             title = title)  + 
        theme(plot.subtitle = element_text(hjust = 0.5, size = 10.5), plot.title = element_text(hjust = 0.5)) +
        ylim(NA, max_y*1.25)
      
    }else{
      
      p <- data %>%  ggboxplot(x = 'Immunogenicity', y = feature, fill = "Immunogenicity", 
                               palette = c("#006400",  "#FC4E07"),  legend = 'none') +
        stat_compare_means(comparisons = my_comparisons, label = "p.signif") + theme_bw() +
        labs(subtitle = count_peptides,   title = title)  + 
        theme(plot.subtitle = element_text(hjust = 0.5, size = 10.5), plot.title = element_text(hjust = 0.5)) +
        ylim(NA, max_y*1.25)
    }
    
    if(remove_legends == T){
      p <- p + theme(legend.position = "none") 
    }
  }
  return(p)
  
}


################ overview of autoantigen dataset ##########################################
autoantigen_data <-  read.csv(file = 'autoantigen_5data.csv', stringsAsFactors=F)
autoantigen_peptides <- autoantigen_data %>% 
  group_by(Peptide) %>%
  dplyr::summarise(MHC_Restriction = paste0(sort(unique(MHC_Restriction)), collapse="|"),
                   wtPeptide = paste0(sort(unique(wtPeptide)), collapse="|"),
                   Antigen_Name = paste0(sort(unique(Antigen_Name)), collapse="|"),
                   Type = paste0(sort(unique(Type)), collapse="|"),
                   data = paste0(sort(unique(data)), collapse="|")) %>%
  filter(!grepl('\\+', Peptide))

autoantigen_data %>% distinct(Peptide, data) %>% plyr::count( vars=c("data")) %>% 
  ggplot(aes(x= reorder(data, -freq), y=freq)) + 
  geom_bar(stat="identity", fill = jcolors('pal7')[1], color = 'black', width = 0.7) + 
  theme_classic() + labs(y = '# of self-peptides', x = 'Database') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
peptide_overlap <- autoantigen_data %>% distinct(Peptide, data) %>%
  plyr::count( vars=c("data", 'Peptide')) %>%
  reshape2::dcast(Peptide ~ data, value.var="freq") %>% rename_all(toupper) %>%
  replace(is.na(.), 0) %>%
  dplyr::rename(CANCER_PEPTIDE_DB = CANCER_PEPTIDE_DATABASE)

library(UpSetR)
upset(peptide_overlap, nsets = ncol(peptide_overlap),
      nintersects = NA, order.by = "freq", line.size = 1.2, point.size = 2.5, text.scale = 1.4,
      # scale.intersections = "log10",
      scale.sets = "log10") 

# histogram of lengths
autoantigen_peptides %>% mutate(length = nchar(as.character(Peptide))) %>%
  plyr::count(vars = c('length') ) %>%
  ggplot(aes(x= length, y=freq)) + geom_bar(stat="identity", position = "dodge2") + 
  # scale_x_continuous(breaks=0:7) +
  theme_classic() + labs(y = '# of autoantigenic peptides', x= 'Length of peptides (aa)') 



########################## compute RSAT ###################################
allsubstr <- function(x, n) unique(substring(x, 1:(nchar(x) - n + 1), n:nchar(x)))

iedbData <- read.csv(file = '1_pathogenic_db.csv') %>% filter(!Antigen_Organism %in% c('Homo sapiens', NA, '')) %>%
  select(Peptide, Immunogenicity, Supertype, Antigen_Organism)
PEPTIDES <- iedbData  %>% pull(Peptide) %>% unique %>% as.data.table() %>% dplyr::rename(Peptide=".")

load('hsProteomeDF.RData')
hsProteomeDF <- str_replace_all(hsProteomeDF, 'U', 'C')
proteome <- AAString(hsProteomeDF)

load(file = 'autoantigen_peptide_vector.RData')
aa1 <- AAString(autoantigen_peptide_vector)


# Run analysis
registerDoParallel(cores=4)
iedb_cognate_epitope = foreach(peptide.i=1:nrow(PEPTIDES),.combine = rbind,.packages = c("dplyr","data.table","magrittr","Biostrings","doParallel")) %dopar% {
  
  dts = foreach(proteome_nrow = 1.,combine = rbind,.packages = c("dplyr","data.table","magrittr","Biostrings")) %do% {
    if(peptide.i == 1){
      print(c(" Processing Peptide Number: ",peptide.i))  
    }
    
    aa2 <- AAString(PEPTIDES[peptide.i,]$Peptide)
    align <- Biostrings::pairwiseAlignment(aa1, aa2, type = "local-global", substitutionMatrix = 'BLOSUM62')
    
    autoantigen_pattern <- AAString(as.character(align@pattern))
    peptid_pattern <- AAString(as.character(align@subject))
    autoantigen_score <- Biostrings::pairwiseAlignment(autoantigen_pattern, autoantigen_pattern, type = "global", substitutionMatrix = 'BLOSUM62', scoreOnly = T)
    peptide_score <- Biostrings::pairwiseAlignment(peptid_pattern, peptid_pattern, type = "global", substitutionMatrix = 'BLOSUM62', scoreOnly = T)
    
    proteome_align <- Biostrings::pairwiseAlignment(proteome, autoantigen_pattern, type = "local-global", substitutionMatrix = 'BLOSUM62')
    proteome_pattern <- AAString(as.character(proteome_align@pattern))
    proteome_score <- Biostrings::pairwiseAlignment(proteome_pattern, proteome_pattern, type = "global", substitutionMatrix = 'BLOSUM62', scoreOnly = T)
    bl_proteoem_iedb_align <- Biostrings::pairwiseAlignment(peptid_pattern, proteome_pattern, type = "global", substitutionMatrix = 'BLOSUM62', scoreOnly = T)
    
    
    DT = data.frame(Peptide =  PEPTIDES[peptide.i,]$Peptide, Peptide_pattern = as.character(align@subject),
                    autoantigen_pattern = as.character(align@pattern), hsproteome_pattern =  as.character(proteome_align@pattern),
                    peptide_score = as.numeric(peptide_score),  autoantigen_score = as.numeric(autoantigen_score),
                    hsproteome_score = as.numeric(proteome_score),  BL_iedb_autoantigen = as.numeric(align@score), 
                    BL_iedb_hsproteome = as.numeric(bl_proteoem_iedb_align)) 
    DT
    
  }
  write.table(dts,
              file = paste0("rsat/PEP_",peptide.i,".txt"),quote=F,row.names = FALSE,sep="\t") # make rsat folder
  
}


setwd('rsat')
temp = list.files(pattern="*.txt")
dt = lapply(temp, read.table, sep="\t", row.names=NULL, colClasses = "character", header = TRUE)
blosum62_autoantigen_iedb <- rbindlist(dt) 
blosum62_autoantigen_iedb <- blosum62_autoantigen_iedb %>% 
  mutate(BL_iedb_autoantigen= as.numeric(BL_iedb_autoantigen),  peptide_score = as.numeric(peptide_score),
         autoantigen_score = as.numeric(autoantigen_score), hsproteome_score = as.numeric(hsproteome_score),
         BL_iedb_hsproteome = as.numeric(BL_iedb_hsproteome)) %>%
  mutate(Matchscore_iedb_autoantigen = BL_iedb_autoantigen/sqrt(peptide_score*autoantigen_score),
         Matchscore_iedb_hsproteome = BL_iedb_hsproteome/sqrt(peptide_score*hsproteome_score)) %>%
  mutate(Ratio_BL_autoantigen_hsproteome = BL_iedb_autoantigen/BL_iedb_hsproteome, 
         Ratio_Matchscore_autoantigen_hsproteome = Matchscore_iedb_autoantigen/Matchscore_iedb_hsproteome) %>%
  left_join(iedbData, by = c('Peptide'))



##################### analysis for similarity with autoantigens ###################################
library(stringdist)
blosum62_autoantigen_iedb <- blosum62_autoantigen_iedb %>%
  filter(!Antigen_Organism %in% c('', 'Homo sapiens') ) %>% 
  mutate(length = nchar(as.character(Peptide))) %>%
  mutate(Matchscore_iedb_autoantigen_binary = ifelse(Matchscore_iedb_autoantigen >= 0.6, 1, 0)) %>%
  mutate(Ratio_BL_iedb_autoantigen_binary = ifelse(Ratio_BL_autoantigen_hsproteome > 1, 1, 0)) %>%
  mutate(Ratio_Matchscore_iedb_autoantigen_binary = ifelse(Ratio_Matchscore_autoantigen_hsproteome > 1, 1, 0)) %>%
  mutate(Ratio_BL_iedb_autoantigen_incomplete = ifelse(Ratio_BL_autoantigen_hsproteome > 1, Ratio_BL_autoantigen_hsproteome, NA)) %>%
  mutate(Ratio_Matchscore_iedb_autoantigen_incomplete = ifelse(Ratio_Matchscore_autoantigen_hsproteome > 1, Ratio_Matchscore_autoantigen_hsproteome, NA)) %>%
  mutate(Ratio_BL_iedb_autoantigen_incomplete_byMatchscore = ifelse(Matchscore_iedb_autoantigen >= 0.6 , Ratio_BL_autoantigen_hsproteome, NA)) %>%
  mutate(Ratio_Matchscore_iedb_autoantigen_incomplete_byMatchscore = ifelse(Matchscore_iedb_autoantigen >= 0.6, Ratio_Matchscore_autoantigen_hsproteome, NA)) %>%
  mutate(hamming_dist = stringdist(Peptide_pattern, autoantigen_pattern, method = 'hamming')) %>% filter(hamming_dist != 0) 
 #filter out 18 peptides overlapping both dataset
9856 -9838 

blosum62_autoantigen_iedb %>% filter(Matchscore_iedb_autoantigen >= 0.6) %>% plyr::count(vars = c('Immunogenicity'))
blosum62_autoantigen_iedb  %>% plyr::count(vars = c('Immunogenicity'))


####### all peptides - FIG 
blosum62_autoantigen_iedb  %>% plot_violin_immunogenicity(feature = 'Matchscore_iedb_autoantigen', 
                                                          title = 'Match score', remove_legends = T) +
  ylim(c(NA, 1.15)) + labs(y = 'Match score\n[pathogenic peptides, autoantigens]') + theme_classic()

# BLOSUM - FIG 
blosum62_autoantigen_iedb  %>% 
  plot_violin_immunogenicity(feature = 'Ratio_BL_iedb_autoantigen_incomplete_byMatchscore', 
                             title = 'Ratio BL score (homologous)', remove_legends =T)+
  ylim(c(NA, 2)) + 
  labs(y = 'Relative similarity to self antigens') + theme_classic()
blosum62_autoantigen_iedb  %>% 
  plot_cumulative_dist_immunogenicity(feature = 'Ratio_BL_iedb_autoantigen_incomplete_byMatchscore', 
                                      title = 'Ratio BL score (homologous)', remove_legends =T) + 
  labs(x = 'Relative similarity to self antigens') + theme_classic()

# Match score 
blosum62_autoantigen_iedb  %>% 
  plot_violin_immunogenicity(feature = 'Ratio_Matchscore_iedb_autoantigen_incomplete_byMatchscore', 
                             title = 'Ratio Match score (homologous)', remove_legends =T) + theme_classic()
blosum62_autoantigen_iedb  %>% 
  plot_cumulative_dist_immunogenicity(feature = 'Ratio_Matchscore_iedb_autoantigen_incomplete_byMatchscore', 
                                      title = 'Ratio Match score (homologous)', remove_legends =T) + theme_classic()


####### output table 
output_table  <- blosum62_autoantigen_iedb %>%
  select(Peptide, Immunogenicity, MHC_Restriction, Supertype, netMHCpan_Rank, NetMHCstabpan_Rank, Antigen_Name, Ratio_BL_iedb_autoantigen_incomplete_byMatchscore)
# write.csv(output_table, file = 'C:/Users/Chloe/Desktop/Oxford-Projects/HashemImmunogenicity/11_final_immunogenicity_model/0_Manuscrip_figures/fig5/RSAT.csv', row.names = FALSE)

####### per length
lengths <- 9:10
length_violin_p <- list()
length_acumul_p <- list()
for(len in 1:length(lengths)){
  length_violin_p[[len]]  <- blosum62_autoantigen_iedb %>%  
    mutate(length = nchar(as.character(Peptide))) %>%
    filter(length %in% lengths[len]) %>%
    plot_violin_immunogenicity(feature = 'Ratio_BL_iedb_autoantigen_incomplete_byMatchscore',  
                               title = lengths[len], remove_legends = T) +
    labs(y = 'Relative similarity to self antigens') + theme_classic()
  
  length_acumul_p[[len]]  <- blosum62_autoantigen_iedb %>%  
    mutate(length = nchar(as.character(Peptide))) %>%
    filter(length %in% lengths[len]) %>%
    plot_cumulative_dist_immunogenicity(feature = 'Ratio_BL_iedb_autoantigen_incomplete_byMatchscore', 
                                        title = lengths[len], remove_legends = T) +
    labs(x = 'Relative similarity to self antigens') + theme_classic()
}
do.call("grid.arrange", c(length_violin_p, ncol=2))
do.call("grid.arrange", c(length_acumul_p, ncol=2))



######## per species effect 
load(file = 'iedb_dissimilar_protome_blosum62.RData')
antigen_organism_col <- blosum62_autoantigen_iedb %>% select(Peptide, Immunogenicity, Antigen_Organism) %>% distinct()
iedb_dissimilar_protome_blosum62 <- iedb_dissimilar_protome_blosum62 %>%
  dplyr::select(-Antigen_Organism) %>%
  left_join(antigen_organism_col, by = c('Peptide', 'Immunogenicity')) %>%
  filter(!Antigen_Organism %in% c('', 'Homo sapiens') )

top_organism <- iedb_dissimilar_protome_blosum62 %>% filter(Antigen_Organism != 'Homo sapiens') %>% 
  filter(Antigen_Organism != "") %>% 
  plyr::count(vars = c('Antigen_Organism', 'Immunogenicity')) %>% 
  filter(freq >= 20) %>% dplyr::select(-freq) %>% plyr::count(vars = 'Antigen_Organism') %>% 
  filter(freq == 2) %>%   pull(Antigen_Organism)

# for incomplete data
org_accum_p <- list()
for(org in 1:length(top_organism)){
  org_accum_p[[org]]  <- blosum62_autoantigen_iedb %>%  
    filter(Antigen_Organism %in% top_organism[org]) %>%
    plot_cumulative_dist_immunogenicity(feature = 'Ratio_BL_iedb_autoantigen_incomplete_byMatchscore',  
                                        title = top_organism[org], remove_legends = T) +
    labs(x = 'Relative similarity to self antigens') + theme_classic()
}
do.call("grid.arrange", c(org_accum_p, ncol=4))
do.call("grid.arrange", c(org_accum_p[c(1, 2, 5, 6, 7)], ncol = 3))


blosum62_autoantigen_iedb %>%  filter(Antigen_Organism %in% top_organism[org]) %>%
  plot_cumulative_dist_immunogenicity(feature = 'Ratio_BL_iedb_autoantigen_incomplete_byMatchscore',  
                                      title = top_organism[org], remove_legends = F) 



