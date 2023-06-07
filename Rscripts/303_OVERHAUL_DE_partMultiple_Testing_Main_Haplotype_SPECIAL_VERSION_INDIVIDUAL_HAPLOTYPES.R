

suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

 
opt = NULL

options(warn = 1)


MT_correction_in_cis_gene = function (option_list)
{
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS))
  cat("\n")
  
  #### Haplotypes ----
  
  
  
  Haplotypes = unlist(strsplit(opt$Haplotypes, split=","))
  
  cat("Haplotypes_\n")
  cat(str(Haplotypes))
  cat("\n")
  
  #### Read previous files ----
  
  setwd(out)
  
  Proxy_file<-as.data.frame(fread(file="Proxy_R2_TABLE_ALL.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("Proxy_file\n")
  cat(str(Proxy_file))
  cat("\n")
  cat(str(unique(Proxy_file$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file$Proxy_VAR)))
  cat("\n")
  
  
  ABSENT_WGS_RNA<-as.data.frame(fread(file="ABSENT_WGS_RNA_Seq.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("ABSENT_WGS_RNA\n")
  cat(str(ABSENT_WGS_RNA))
  cat("\n")
  
  
  #### UPDATE Selected_vars ----
  
  indx.UPDATE<-which(SELECTED_VARS%in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS[-indx.UPDATE]
    
  }else{
    
    SELECTED_VARS_UPDATED = SELECTED_VARS
  }
  
  
  
  cat("SELECTED_VARS_UPDATED_\n")
  cat(str(SELECTED_VARS_UPDATED))
  cat("\n")
  
  
  
  Proxy_file_UPDATED<- Proxy_file[-which(Proxy_file$VAR%in%ABSENT_WGS_RNA$VAR |
                                           Proxy_file$Proxy_VAR%in%ABSENT_WGS_RNA$VAR ),]
  
  
  
  cat("Proxy_file_UPDATED\n")
  cat(str(Proxy_file_UPDATED))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$Proxy_VAR)))
  cat("\n")
  
 ######################### MASTER LOOP #############################
  #### LOOP ----
  
  Condition_DEBUG <- 0
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    SELECTED_VARS_UPDATED_sel<-SELECTED_VARS_UPDATED
    
    # cat("MAIN_VAR---------->\t")
    # cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
    # cat("\n")
    
    
    cat("MAIN_VAR---------->\t")
    cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
    cat("\n")
    
    Proxy_file_UPDATED_sel<-Proxy_file_UPDATED[which(Proxy_file_UPDATED$VAR%in%SELECTED_VARS_UPDATED_sel),]
    
    if(dim(Proxy_file_UPDATED_sel)[1] >0)
    {
      if(Condition_DEBUG == 1)
      {
        cat("Proxy_file_UPDATED_sel\n")
        cat(str(Proxy_file_UPDATED_sel))
        cat("\n")
        cat(str(unique(Proxy_file_UPDATED_sel$VAR)))
        cat("\n")
        cat(str(unique(Proxy_file_UPDATED_sel$Proxy_VAR)))
        cat("\n")
      }
      
      Proxy_array<-unique(Proxy_file_UPDATED_sel$Proxy_VAR)
      
      if(Condition_DEBUG == 1)
      {
        cat("Proxy_array\n")
        cat(str(Proxy_array))
        cat("\n")
      }
      
      Condition_DEBUG <- 1
      for(k in 1:length(Proxy_array))
      {
        Proxy_array_sel<-Proxy_array[k]
        
        
        cat("Proxy_VAR---------->\t")
        cat(sprintf(as.character(Proxy_array_sel)))
        cat("\n")
        
        
        Proxy_df<-Proxy_file_UPDATED_sel[which(Proxy_file_UPDATED_sel$Proxy_VAR == Proxy_array_sel),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Proxy_df\n")
          cat(str(Proxy_df))
          cat("\n")
          cat(str(unique(Proxy_df$VAR)))
          cat("\n")
          cat(str(unique(Proxy_df$Proxy_VAR)))
          cat("\n")
        }
        
        Haplotype_candidate<-paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__")
        
        FLAG<-Haplotype_candidate[which(Haplotype_candidate%in%Haplotypes)]
        
        if(length(FLAG) >0)
        {
          path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
          
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path10)){
            setwd(path10)
            
            filename_rds_haplotype<-paste("DE_LM_Haplotypes_RESULTS_NOMINAL_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'.rds', sep='')
            
            if (file.exists(filename_rds_haplotype)){
              
              Results_Nominal<-readRDS(file=filename_rds_haplotype)
              
              
              Results_Nominal$comparison<-factor(Results_Nominal$comparison,
                                                 levels=c('HOM_REF__HET|HOM_REF','HOM_REF__HOM_REF|HET','HOM_REF__HET|HET','HOM_REF__HOM|HOM_REF','HOM_REF__HOM_REF|HOM','HOM_REF__HOM|HET','HOM_REF__HET|HOM','HOM_REF__HOM|HOM'),
                                                 ordered=T)
              
              if(Condition_DEBUG == 1)
              {
                cat("Results_Nominal\n")
                cat(str(Results_Nominal))
                cat("\n")
                
              }
              
              #### VEP_CSQ ----
              
              VEP_CSQ<-as.data.frame(fread(file=opt$VEP_CSQ,sep=",") , stringsAsFactors=F)
              
              cat("VEP_CSQ_0\n")
              cat(str(VEP_CSQ))
              cat("\n")
              
              VEP_CSQ_sel<-VEP_CSQ[which(VEP_CSQ$VAR%in%c(SELECTED_VARS_UPDATED_sel,Proxy_array_sel)),]
              
              cat("VEP_CSQ_sel_0\n")
              cat(str(VEP_CSQ_sel))
              cat("\n")
              
              rm(VEP_CSQ)
              
              
              #### CIS consequences -----
              
              CIS_LABELS <-c("LOF","MISS","SYN","UTR5","UTR3",
                             "INTRON","SPLICE","UPSTREAM")
              
              
              VEP_CSQ_sel_CIS<-VEP_CSQ_sel[which(VEP_CSQ_sel$VEP_DEF_LABELS%in%CIS_LABELS),]
              
              cat("VEP_CSQ_sel_CIS_0\n")
              cat(str(VEP_CSQ_sel_CIS))
              cat("\n")
              
              if(dim(VEP_CSQ_sel_CIS)[1] >0)
              {
                
                #### SELECTED GENES ----
                
                ENSG_array<-vector()
                
                
                
                if(dim(VEP_CSQ_sel_CIS)[1] >0)
                {
                  
                  ENSG_sel<-unique(VEP_CSQ_sel_CIS$ensembl_gene_id)
                  
                  
                  ENSG_array<-c(ENSG_array,ENSG_sel)
                  
                }
                
                if(Condition_DEBUG == 1)
                {
                  cat("ENSG_array\n")
                  cat(str(ENSG_array))
                  cat("\n")
                }
                
                
                if(length(ENSG_array) >0)
                {
                  Results_Nominal_sel<-droplevels(Results_Nominal[which(Results_Nominal$ensembl_gene_id%in%ENSG_array),])
                  
                  rm(Results_Nominal)
                  
                  if(Condition_DEBUG == 1)
                  { 
                    cat("Results_Nominal_sel\n")
                    cat(str(Results_Nominal_sel))
                    cat("\n")
                    cat(sprintf(as.character(names(summary(as.factor(Results_Nominal_sel$comparison))))))
                    cat("\n")
                    cat(sprintf(as.character(summary(as.factor(Results_Nominal_sel$comparison)))))
                    cat("\n")
                    
                  }
                  
                  if(dim(Results_Nominal_sel)[1] >0)
                  {
                    array_comparisons<-levels(Results_Nominal_sel$comparison)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("array_comparisons\n")
                      cat(str(array_comparisons))
                      cat("\n")
                    }
                    
                    list_comparisons<-list()
                    
                    Condition_DEBUG <- 1
                    
                    for(l in 1:length(array_comparisons))
                    {
                      array_haplotypes_sel<-array_comparisons[l]
                      
                      
                      cat("---------->\t")
                      cat(sprintf(as.character(array_haplotypes_sel)))
                      cat("\n")
                      
                      Results_Nominal_sel_comparison_sel<-Results_Nominal_sel[which(Results_Nominal_sel$comparison == array_haplotypes_sel),]
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("Results_Nominal_sel_comparison_sel\n")
                        cat(str(Results_Nominal_sel_comparison_sel))
                        cat("\n")
                      }
                      
                      indx.select<-c(which(colnames(Results_Nominal_sel_comparison_sel)== "ensembl_gene_id"),which(colnames(Results_Nominal_sel_comparison_sel)== "HGNC"),which(colnames(Results_Nominal_sel_comparison_sel)== "comparison"),
                                     which(colnames(Results_Nominal_sel_comparison_sel)== "coefficient_Haplotypes_specific_CELL_COUNTS"),
                                     which(colnames(Results_Nominal_sel_comparison_sel)== "pvalue_Haplotypes_specific_CELL_COUNTS"),
                                     which(colnames(Results_Nominal_sel_comparison_sel)== "n_breakdown_string"))
                      
                      
                      
                      MT_set<-unique(Results_Nominal_sel_comparison_sel[,indx.select])
                      
                      if(Condition_DEBUG == 1)
                      { 
                        cat("MT_set_0\n")
                        cat(str(MT_set))
                        cat("\n")
                        cat(str(unique(MT_set$ensembl_gene_id)))
                        cat("\n")
                      }
                      
                      MT_set$ajusted.pvalue_Haplotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Haplotypes_specific_CELL_COUNTS, method = "BH")
                      MT_set$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Haplotypes_specific_CELL_COUNTS),2)
                      
                      if(Condition_DEBUG == 1)
                      { 
                        cat("MT_set_1\n")
                        cat(str(MT_set))
                        cat("\n")
                      }
                      
                      list_comparisons[[l]]<-MT_set
                      
                    }#l in 1:length(array_comparisons)
                    
                    if(length(list_comparisons) >0)
                    {
                      
                      Results_corrected_all_comparisons = as.data.frame(data.table::rbindlist(list_comparisons, fill=T), stringsAsFactors=F)
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("Results_uncorrected_all_comparisons_0\n")
                        cat(str(Results_corrected_all_comparisons))
                        cat("\n")
                        cat(str(unique(Results_corrected_all_comparisons$ensembl_gene_id)))
                        cat("\n")
                        
                      }
                      
                      DE_CIS_gene<-Results_corrected_all_comparisons[which(Results_corrected_all_comparisons$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
                      
                      cat("DE_CIS_gene_1\n")
                      cat(str(DE_CIS_gene))
                      cat("\n")
                      cat(str(unique(DE_CIS_gene$ensembl_gene_id)))
                      cat("\n")
                      
                      
                      # quit(status = 1)
                      path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                      
                      
                      setwd(path10)
                      if(dim(Results_corrected_all_comparisons)[1] >0)
                      {
                        saveRDS(file="DE_CIS_gene.rds", Results_corrected_all_comparisons)  
                      }
                    }#length(list_comparisons) >0
                  }#dim(Results_Nominal_sel)[1] >0
                }#length(ENSG_array) >0
              }else{
                
                cat("Hello_world\n")
              }#dim(VEP_CSQ_sel_CIS)[1] >0
              
              
            }else{
              
            }#file.exists(filename_rds_haplotype)
            
            
            
          } else {
            
            
          }#file.exists(path10))
        }#length(FLAG) >0
      }#k in 1:length(Proxy_array)
      # ############################################################################################################
      # quit(status = 1)
      
    }#dim(Proxy_file_UPDATED_sel)[1] >0
  }# length(SELECTED_VARS_UPDATED) >0
}


MT_correction_in_block_Plus_PCHi_C = function (option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS))
  cat("\n")
  
  #### Haplotypes ----
  
  
  
  Haplotypes = unlist(strsplit(opt$Haplotypes, split=","))
  
  cat("Haplotypes_\n")
  cat(str(Haplotypes))
  cat("\n")
  
  #### Read previous files ----
  
  setwd(out)
  
  Proxy_file<-as.data.frame(fread(file="Proxy_R2_TABLE_ALL.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("Proxy_file\n")
  cat(str(Proxy_file))
  cat("\n")
  cat(str(unique(Proxy_file$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file$Proxy_VAR)))
  cat("\n")
  
  
  ABSENT_WGS_RNA<-as.data.frame(fread(file="ABSENT_WGS_RNA_Seq.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("ABSENT_WGS_RNA\n")
  cat(str(ABSENT_WGS_RNA))
  cat("\n")
  
  
  #### UPDATE Selected_vars ----
  
  indx.UPDATE<-which(SELECTED_VARS%in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS[-indx.UPDATE]
    
  }else{
    
    SELECTED_VARS_UPDATED = SELECTED_VARS
  }
  
  
  
  cat("SELECTED_VARS_UPDATED_\n")
  cat(str(SELECTED_VARS_UPDATED))
  cat("\n")
  
  
  
  Proxy_file_UPDATED<- Proxy_file[-which(Proxy_file$VAR%in%ABSENT_WGS_RNA$VAR |
                                           Proxy_file$Proxy_VAR%in%ABSENT_WGS_RNA$VAR ),]
  
  
  
  cat("Proxy_file_UPDATED\n")
  cat(str(Proxy_file_UPDATED))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$Proxy_VAR)))
  cat("\n")
  
  ######################### MASTER LOOP #############################
  #### LOOP ----
  
  Condition_DEBUG <- 0
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    SELECTED_VARS_UPDATED_sel<-SELECTED_VARS_UPDATED
    
    # cat("MAIN_VAR---------->\t")
    # cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
    # cat("\n")
    
    
    cat("MAIN_VAR---------->\t")
    cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
    cat("\n")
    
    Proxy_file_UPDATED_sel<-Proxy_file_UPDATED[which(Proxy_file_UPDATED$VAR%in%SELECTED_VARS_UPDATED_sel),]
    
    if(dim(Proxy_file_UPDATED_sel)[1] >0)
    {
      if(Condition_DEBUG == 1)
      {
        cat("Proxy_file_UPDATED_sel\n")
        cat(str(Proxy_file_UPDATED_sel))
        cat("\n")
        cat(str(unique(Proxy_file_UPDATED_sel$VAR)))
        cat("\n")
        cat(str(unique(Proxy_file_UPDATED_sel$Proxy_VAR)))
        cat("\n")
      }
      
      Proxy_array<-unique(Proxy_file_UPDATED_sel$Proxy_VAR)
      
      if(Condition_DEBUG == 1)
      {
        cat("Proxy_array\n")
        cat(str(Proxy_array))
        cat("\n")
      }
      
      Condition_DEBUG <- 0
      for(k in 1:length(Proxy_array))
      {
        Proxy_array_sel<-Proxy_array[k]
        
        
        cat("Proxy_VAR---------->\t")
        cat(sprintf(as.character(Proxy_array_sel)))
        cat("\n")
        
        
        Proxy_df<-Proxy_file_UPDATED_sel[which(Proxy_file_UPDATED_sel$Proxy_VAR == Proxy_array_sel),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Proxy_df\n")
          cat(str(Proxy_df))
          cat("\n")
          cat(str(unique(Proxy_df$VAR)))
          cat("\n")
          cat(str(unique(Proxy_df$Proxy_VAR)))
          cat("\n")
        }
        
        Haplotype_candidate<-paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__")
        
        FLAG<-Haplotype_candidate[which(Haplotype_candidate%in%Haplotypes)]
        
        if(length(FLAG) >0)
        {
          path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
          
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path10)){
            setwd(path10)
            
            filename_rds_haplotype<-paste("DE_LM_Haplotypes_RESULTS_NOMINAL_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'.rds', sep='')
            
            if (file.exists(filename_rds_haplotype)){
              
              Results_Nominal<-readRDS(file=filename_rds_haplotype)
              
              
              Results_Nominal$comparison<-factor(Results_Nominal$comparison,
                                                 levels=c('HOM_REF__HET|HOM_REF','HOM_REF__HOM_REF|HET','HOM_REF__HET|HET','HOM_REF__HOM|HOM_REF','HOM_REF__HOM_REF|HOM','HOM_REF__HOM|HET','HOM_REF__HET|HOM','HOM_REF__HOM|HOM'),
                                                 ordered=T)
              
              if(Condition_DEBUG == 1)
              {
                cat("Results_Nominal\n")
                cat(str(Results_Nominal))
                cat("\n")
                
              }
              
              #### Read ALL_dB file ----
              
              ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
              
              if(Condition_DEBUG == 1)
              {
                cat("ALL_dB_0\n")
                cat(str(ALL_dB))
                cat("\n")
              }
              
              
              indx.dep<-c(which(colnames(ALL_dB) == "maf_origin"))
              
              ALL_dB_subset<-unique(ALL_dB[,-indx.dep])
              
              # cat("ALL_dB_subset\n")
              # cat(str(ALL_dB_subset))
              # cat("\n")
              
              rm(ALL_dB)
              
              
              ALL_dB_subset$Allelic_Series_ID<-paste(ALL_dB_subset$phenotype,ALL_dB_subset$block_no,sep='__')
              
              # cat("ALL_dB_subset_1\n")
              # cat(str(ALL_dB_subset))
              # cat("\n")
              
              
              indx.int<-c(which(colnames(ALL_dB_subset) == "VAR"),which(colnames(ALL_dB_subset) == "Allelic_Series_ID"))
              
              ALL_FINAL<-unique(ALL_dB_subset[,indx.int])
              
              if(Condition_DEBUG == 1)
              {
                cat("ALL_FINAL\n")
                cat(str(ALL_FINAL))
                cat("\n")
              }
              
              
              ALL_FINAL_sel<-ALL_FINAL[which(ALL_FINAL$VAR%in%c(SELECTED_VARS_UPDATED,Proxy_array_sel)),]
              
              rm(ALL_FINAL)
              
              if(Condition_DEBUG == 1)
              {
                cat("ALL_FINAL_sel\n")
                cat(str(ALL_FINAL_sel))
                cat("\n")
              }
              
              if(dim(ALL_FINAL_sel)[1] >0)
              {
                
                PCHiC_info<-as.data.frame(fread(file=opt$PCHiC_info, sep =",", header =T), stringsAsFactors=F)
                
                if(Condition_DEBUG == 1)
                {
                  cat("PCHiC_info:\n")
                  cat(str(PCHiC_info))
                  cat("\n")
                }
                
                PCHiC_info_sel<-PCHiC_info[which(PCHiC_info$VAR%in%c(SELECTED_VARS_UPDATED,Proxy_array_sel)),]
                
                rm(PCHiC_info)
                
                if(Condition_DEBUG == 1)
                {
                  cat("PCHiC_info_sel:\n")
                  cat(str(PCHiC_info_sel))
                  cat("\n")
                }
                
                #### GENES_PER_BLOCKS ----
                
                GENES_PER_BLOCKS = as.data.frame(fread(file=opt$GENES_PER_BLOCKS, sep="\t", stringsAsFactors = F, header = T))
                
                # cat("GENES_PER_BLOCKS\n")
                # cat(str(GENES_PER_BLOCKS))
                # cat("\n")
                
                GENES_PER_BLOCKS$phenotype<-gsub("^.+__","",GENES_PER_BLOCKS$BLOCK, perl=T)
                GENES_PER_BLOCKS$block_no<-gsub("__.+$","",GENES_PER_BLOCKS$BLOCK, perl=T)
                
                GENES_PER_BLOCKS$Allelic_Series_ID<-paste(GENES_PER_BLOCKS$phenotype,GENES_PER_BLOCKS$block_no,sep="__")
                
                GENES_PER_BLOCKS<-GENES_PER_BLOCKS[,-c(which(colnames(GENES_PER_BLOCKS) == "phenotype"),
                                                       which(colnames(GENES_PER_BLOCKS) == "block_no"))]
                
                
                # cat("GENES_PER_BLOCKS\n")
                # cat(str(GENES_PER_BLOCKS))
                # cat("\n")
                # 
                
                AS_Gene_source<-merge(ALL_FINAL_sel,
                                      GENES_PER_BLOCKS,
                                      by="Allelic_Series_ID")
                
                rm(GENES_PER_BLOCKS)
                
                if(Condition_DEBUG == 1)
                {
                  cat("AS_Gene_source\n")
                  cat(str(AS_Gene_source))
                  cat("\n")
                }
                
                
                #### SELECTED GENES ----
                
                ENSG_array<-vector()
                
                
                AS_Gene_source_sel<-AS_Gene_source[which(AS_Gene_source$VAR%in%c(SELECTED_VARS_UPDATED,Proxy_array_sel)),]
                
                if(Condition_DEBUG == 1)
                {
                  cat("AS_Gene_source_sel\n")
                  cat(str(AS_Gene_source_sel))
                  cat("\n")
                }
                
                if(dim(AS_Gene_source_sel)[1] >0)
                {
                  
                  ENSG_sel<-unique(AS_Gene_source_sel$ensembl_gene_id)
                  
                  
                  ENSG_array<-c(ENSG_array,ENSG_sel)
                  
                }
                
                if(dim(PCHiC_info_sel)[1] >0)
                {
                  if(Condition_DEBUG == 1)
                  {
                    cat("PCHiC_info_sel\n")
                    cat(str(PCHiC_info_sel))
                    cat("\n")
                  }
                  
                  ENSG_sel<-unique(PCHiC_info_sel$ensembl_gene_id)
                  
                  
                  ENSG_array<-c(ENSG_array,ENSG_sel)
                }
                
                # if(SELECTED_VARS_UPDATED == "chr9_135920196_C_T")
                # {
                #   ENSG_array<-unique(c(ENSG_array,"ENSG00000165702"))
                # }
                
                if(Condition_DEBUG == 1)
                {
                  cat("ENSG_array\n")
                  cat(str(ENSG_array))
                  cat("\n")
                }
                
                if(length(ENSG_array) >0)
                {
                  Results_Nominal_sel<-droplevels(Results_Nominal[which(Results_Nominal$ensembl_gene_id%in%ENSG_array),])
                  
                  rm(Results_Nominal)
                  
                  if(Condition_DEBUG == 1)
                  { 
                    cat("Results_Nominal_sel\n")
                    cat(str(Results_Nominal_sel))
                    cat("\n")
                    cat(sprintf(as.character(names(summary(as.factor(Results_Nominal_sel$comparison))))))
                    cat("\n")
                    cat(sprintf(as.character(summary(as.factor(Results_Nominal_sel$comparison)))))
                    cat("\n")
                    
                  }
                  
                  if(dim(Results_Nominal_sel)[1] >0)
                  {
                    array_comparisons<-levels(Results_Nominal_sel$comparison)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("array_comparisons\n")
                      cat(str(array_comparisons))
                      cat("\n")
                    }
                    
                    list_comparisons<-list()
                    
                    Condition_DEBUG <- 1
                    
                    for(l in 1:length(array_comparisons))
                    {
                      array_haplotypes_sel<-array_comparisons[l]
                      
                      
                      cat("---------->\t")
                      cat(sprintf(as.character(array_haplotypes_sel)))
                      cat("\n")
                      
                      Results_Nominal_sel_comparison_sel<-Results_Nominal_sel[which(Results_Nominal_sel$comparison == array_haplotypes_sel),]
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("Results_Nominal_sel_comparison_sel\n")
                        cat(str(Results_Nominal_sel_comparison_sel))
                        cat("\n")
                      }
                      
                      indx.select<-c(which(colnames(Results_Nominal_sel_comparison_sel)== "ensembl_gene_id"),which(colnames(Results_Nominal_sel_comparison_sel)== "HGNC"),which(colnames(Results_Nominal_sel_comparison_sel)== "comparison"),
                                     which(colnames(Results_Nominal_sel_comparison_sel)== "coefficient_Haplotypes_specific_CELL_COUNTS"),
                                     which(colnames(Results_Nominal_sel_comparison_sel)== "pvalue_Haplotypes_specific_CELL_COUNTS"),
                                     which(colnames(Results_Nominal_sel_comparison_sel)== "n_breakdown_string"))
                      
                      
                      
                      MT_set<-unique(Results_Nominal_sel_comparison_sel[,indx.select])
                      
                      if(Condition_DEBUG == 1)
                      { 
                        cat("MT_set_0\n")
                        cat(str(MT_set))
                        cat("\n")
                        cat(str(unique(MT_set$ensembl_gene_id)))
                        cat("\n")
                      }
                      
                      MT_set$ajusted.pvalue_Haplotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Haplotypes_specific_CELL_COUNTS, method = "BH")
                      MT_set$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Haplotypes_specific_CELL_COUNTS),2)
                      
                      if(Condition_DEBUG == 1)
                      { 
                        cat("MT_set_1\n")
                        cat(str(MT_set))
                        cat("\n")
                      }
                      
                      list_comparisons[[l]]<-MT_set
                      
                    }#l in 1:length(array_comparisons)
                    
                    if(length(list_comparisons) >0)
                    {
                      
                      Results_corrected_all_comparisons = as.data.frame(data.table::rbindlist(list_comparisons, fill=T), stringsAsFactors=F)
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("Results_uncorrected_all_comparisons_0\n")
                        cat(str(Results_corrected_all_comparisons))
                        cat("\n")
                        cat(str(unique(Results_corrected_all_comparisons$ensembl_gene_id)))
                        cat("\n")
                        
                      }
                      
                      DE_Block<-Results_corrected_all_comparisons[which(Results_corrected_all_comparisons$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
                      
                      cat("DE_Block_1\n")
                      cat(str(DE_Block))
                      cat("\n")
                      cat(str(unique(DE_Block$ensembl_gene_id)))
                      cat("\n")
                      
                      
                      # quit(status = 1)
                      path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                      
                      
                      setwd(path10)
                      if(dim(Results_corrected_all_comparisons)[1] >0)
                      {
                        saveRDS(file="DE_Block.rds", Results_corrected_all_comparisons)  
                      }
                    }#length(list_comparisons) >0
                  }#dim(Results_Nominal_sel)[1] >0
                }#length(ENSG_array) >0
              }#dim(ALL_FINAL_sel)[1] >0
              
              
              
            }else{
              
            }#file.exists(filename_rds_haplotype)
            
            
            
          } else {
            
            
          }#file.exists(path10))
        }#length(FLAG) >0
      }#k in 1:length(Proxy_array)
      # ############################################################################################################
      # quit(status = 1)
      
    }#dim(Proxy_file_UPDATED_sel)[1] >0
  }# length(SELECTED_VARS_UPDATED) >0

}



MT_correction_genome_wide = function (option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS))
  cat("\n")
  
  #### Haplotypes ----
  
  
  
  Haplotypes = unlist(strsplit(opt$Haplotypes, split=","))
  
  cat("Haplotypes_\n")
  cat(str(Haplotypes))
  cat("\n")
  
  #### Read previous files ----
  
  setwd(out)
  
  Proxy_file<-as.data.frame(fread(file="Proxy_R2_TABLE_ALL.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("Proxy_file\n")
  cat(str(Proxy_file))
  cat("\n")
  cat(str(unique(Proxy_file$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file$Proxy_VAR)))
  cat("\n")
  
  
  ABSENT_WGS_RNA<-as.data.frame(fread(file="ABSENT_WGS_RNA_Seq.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("ABSENT_WGS_RNA\n")
  cat(str(ABSENT_WGS_RNA))
  cat("\n")
  
  
  #### UPDATE Selected_vars ----
  
  indx.UPDATE<-which(SELECTED_VARS%in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS[-indx.UPDATE]
    
  }else{
    
    SELECTED_VARS_UPDATED = SELECTED_VARS
  }
  
  
  
  cat("SELECTED_VARS_UPDATED_\n")
  cat(str(SELECTED_VARS_UPDATED))
  cat("\n")
  
  
  
  Proxy_file_UPDATED<- Proxy_file[-which(Proxy_file$VAR%in%ABSENT_WGS_RNA$VAR |
                                           Proxy_file$Proxy_VAR%in%ABSENT_WGS_RNA$VAR ),]
  
  
  
  cat("Proxy_file_UPDATED\n")
  cat(str(Proxy_file_UPDATED))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$Proxy_VAR)))
  cat("\n")
  
  ######################### MASTER LOOP #############################
  #### LOOP ----
  
  Condition_DEBUG <- 0
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    SELECTED_VARS_UPDATED_sel<-SELECTED_VARS_UPDATED
    
    # cat("MAIN_VAR---------->\t")
    # cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
    # cat("\n")
    
    
    cat("MAIN_VAR---------->\t")
    cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
    cat("\n")
    
    Proxy_file_UPDATED_sel<-Proxy_file_UPDATED[which(Proxy_file_UPDATED$VAR%in%SELECTED_VARS_UPDATED_sel),]
    
    if(dim(Proxy_file_UPDATED_sel)[1] >0)
    {
      if(Condition_DEBUG == 1)
      {
        cat("Proxy_file_UPDATED_sel\n")
        cat(str(Proxy_file_UPDATED_sel))
        cat("\n")
        cat(str(unique(Proxy_file_UPDATED_sel$VAR)))
        cat("\n")
        cat(str(unique(Proxy_file_UPDATED_sel$Proxy_VAR)))
        cat("\n")
      }
      
      Proxy_array<-unique(Proxy_file_UPDATED_sel$Proxy_VAR)
      
      if(Condition_DEBUG == 1)
      {
        cat("Proxy_array\n")
        cat(str(Proxy_array))
        cat("\n")
      }
      
      Condition_DEBUG <- 0
      for(k in 1:length(Proxy_array))
      {
        Proxy_array_sel<-Proxy_array[k]
        
        
        cat("Proxy_VAR---------->\t")
        cat(sprintf(as.character(Proxy_array_sel)))
        cat("\n")
        
        
        Proxy_df<-Proxy_file_UPDATED_sel[which(Proxy_file_UPDATED_sel$Proxy_VAR == Proxy_array_sel),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Proxy_df\n")
          cat(str(Proxy_df))
          cat("\n")
          cat(str(unique(Proxy_df$VAR)))
          cat("\n")
          cat(str(unique(Proxy_df$Proxy_VAR)))
          cat("\n")
        }
        
        Haplotype_candidate<-paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__")
        
        FLAG<-Haplotype_candidate[which(Haplotype_candidate%in%Haplotypes)]
        
        if(length(FLAG) >0)
        {
          path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
          
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path10)){
            setwd(path10)
            
            filename_rds_haplotype<-paste("DE_LM_Haplotypes_RESULTS_NOMINAL_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'.rds', sep='')
            
            if (file.exists(filename_rds_haplotype)){
              
              Results_Nominal<-readRDS(file=filename_rds_haplotype)
              
              
              Results_Nominal$comparison<-factor(Results_Nominal$comparison,
                                                 levels=c('HOM_REF__HET|HOM_REF','HOM_REF__HOM_REF|HET','HOM_REF__HET|HET','HOM_REF__HOM|HOM_REF','HOM_REF__HOM_REF|HOM','HOM_REF__HOM|HET','HOM_REF__HET|HOM','HOM_REF__HOM|HOM'),
                                                 ordered=T)
              
              Results_Nominal$comparison<-droplevels(Results_Nominal$comparison)
              if(Condition_DEBUG == 1)
              {
                cat("Results_Nominal\n")
                cat(str(Results_Nominal))
                cat("\n")
                
              }
              
              array_comparisons<-levels(Results_Nominal$comparison)
              
              if(Condition_DEBUG == 1)
              {
                cat("array_comparisons\n")
                cat(str(array_comparisons))
                cat("\n")
              }
              
              list_comparisons<-list()
              
              Condition_DEBUG <- 1
              
              for(l in 1:length(array_comparisons))
              {
                array_haplotypes_sel<-array_comparisons[l]
                
                
                cat("---------->\t")
                cat(sprintf(as.character(array_haplotypes_sel)))
                cat("\n")
                
                Results_Nominal_comparison_sel<-Results_Nominal[which(Results_Nominal$comparison == array_haplotypes_sel),]
                
                if(Condition_DEBUG == 1)
                {
                  cat("Results_Nominal_comparison_sel\n")
                  cat(str(Results_Nominal_comparison_sel))
                  cat("\n")
                }
                
                indx.select<-c(which(colnames(Results_Nominal_comparison_sel)== "ensembl_gene_id"),which(colnames(Results_Nominal_comparison_sel)== "HGNC"),which(colnames(Results_Nominal_comparison_sel)== "comparison"),
                               which(colnames(Results_Nominal_comparison_sel)== "coefficient_Haplotypes_specific_CELL_COUNTS"),
                               which(colnames(Results_Nominal_comparison_sel)== "pvalue_Haplotypes_specific_CELL_COUNTS"),
                               which(colnames(Results_Nominal_comparison_sel)== "n_breakdown_string"))
                
                
                
                MT_set<-unique(Results_Nominal_comparison_sel[,indx.select])
                
                if(Condition_DEBUG == 1)
                { 
                  cat("MT_set_0\n")
                  cat(str(MT_set))
                  cat("\n")
                  cat(str(unique(MT_set$ensembl_gene_id)))
                  cat("\n")
                }
                
                MT_set$ajusted.pvalue_Haplotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Haplotypes_specific_CELL_COUNTS, method = "BH")
                MT_set$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Haplotypes_specific_CELL_COUNTS),2)
                
                if(Condition_DEBUG == 1)
                { 
                  cat("MT_set_1\n")
                  cat(str(MT_set))
                  cat("\n")
                }
                
                list_comparisons[[l]]<-MT_set
                
              }#l in 1:length(array_comparisons)
              
              if(length(list_comparisons) >0)
              {
                
                Results_corrected_all_comparisons = as.data.frame(data.table::rbindlist(list_comparisons, fill=T), stringsAsFactors=F)
                
                
                if(Condition_DEBUG == 1)
                {
                  cat("Results_uncorrected_all_comparisons_0\n")
                  cat(str(Results_corrected_all_comparisons))
                  cat("\n")
                  cat(str(unique(Results_corrected_all_comparisons$ensembl_gene_id)))
                  cat("\n")
                  
                }
                
                DE_genome_wide<-Results_corrected_all_comparisons[which(Results_corrected_all_comparisons$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
                
                cat("DE_genome_wide_1\n")
                cat(str(DE_genome_wide))
                cat("\n")
                cat(str(unique(DE_genome_wide$ensembl_gene_id)))
                cat("\n")
                
                
                # quit(status = 1)
                path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                
                
                setwd(path10)
                if(dim(Results_corrected_all_comparisons)[1] >0)
                {
                  saveRDS(file="DE_genome_wide.rds", Results_corrected_all_comparisons)  
                }
              }#length(list_comparisons) >0
            }else{
              
            }#file.exists(filename_rds_haplotype)
            
            
            
          } else {
            
            
          }#file.exists(path10))
          
        }#length(FLAG) >0
      }#k in 1:length(Proxy_array)
      # ############################################################################################################
      # quit(status = 1)
      
    }#dim(Proxy_file_UPDATED_sel)[1] >0
  }# length(SELECTED_VARS_UPDATED) >0
}



printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----




main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--SELECTED_VARS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VEP_CSQ"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC_info"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Haplotypes"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GENES_PER_BLOCKS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)

  MT_correction_in_cis_gene(opt)
  MT_correction_in_block_Plus_PCHi_C(opt)
  MT_correction_genome_wide(opt)
  
}




###########################################################################

system.time( main() )
