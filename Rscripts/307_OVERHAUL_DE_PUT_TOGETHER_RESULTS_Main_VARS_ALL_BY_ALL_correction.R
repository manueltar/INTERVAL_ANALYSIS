

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


collate_results_function = function(option_list)
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
  
  
  
  SELECTED_VARS_UPDATED  = unlist(strsplit(opt$SELECTED_VARS_UPDATED , split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS_UPDATED ))
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
  
  indx.UPDATE<-which(SELECTED_VARS_UPDATED %in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS_UPDATED [-indx.UPDATE]
    
  }else{
    
    SELECTED_VARS_UPDATED = SELECTED_VARS_UPDATED 
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
  
  Condition_DEBUG <- 0
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    
    ##### LOOP TO READ ALL VARIABLES -----
    list_GW<-list()
    list_Block<-list()
    list_CIS<-list()
    
    
    
    for(i in 1:length(SELECTED_VARS_UPDATED))
    {
      
      SELECTED_VARS_sel<-SELECTED_VARS_UPDATED[i]
      
      cat("--------------->\t")
      cat(sprintf(as.character(SELECTED_VARS_sel)))
      cat("\n")
      
      
      ##### path6 ---
      
      path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
      
      # cat("path6\n")
      # cat(sprintf(as.character(path6)))
      # cat("\n")
      
      
      setwd(path6)
      
      filename="DE_genome_wide.rds"
      
      
      
      if (file.exists(filename)) {
        DE_genome_wide = readRDS(filename)
        
        # quit(status = 1)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_genome_wide_0\n")
          cat(str(DE_genome_wide))
          cat("\n")
        }
        
        DE_genome_wide_subset<-unique(DE_genome_wide[,-c(which(colnames(DE_genome_wide) == "ajusted.pvalue_Genotypes_specific_CELL_COUNTS"),
                                                         which(colnames(DE_genome_wide) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))])
        
        DE_genome_wide_subset$VAR<-SELECTED_VARS_sel
        if(Condition_DEBUG == 1)
        {
          cat("DE_genome_wide_subset_0\n")
          cat(str(DE_genome_wide_subset))
          cat("\n")
        }
        
        list_GW[[i]]<-DE_genome_wide_subset
        
        # ##########################################
        # quit(status = 1)
        
      }#file.exists(DE_genome_wide.rds)
      
      ### DTU block
      
      setwd(path6)
      
      filename="DE_Block.rds"
      
      
      
      if (file.exists(filename)) {
        DE_BLOCK = readRDS(filename)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_BLOCK_0\n")
          cat(str(DE_BLOCK))
          cat("\n")
        }
        
        DE_BLOCK_subset<-unique(DE_BLOCK[,-c(which(colnames(DE_BLOCK) == "ajusted.pvalue_Genotypes_specific_CELL_COUNTS"),
                                             which(colnames(DE_BLOCK) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))])
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_BLOCK_0\n")
          cat(str(DE_BLOCK))
          cat("\n")
        }
        
        DE_BLOCK_subset<-unique(DE_BLOCK[,-c(which(colnames(DE_BLOCK) == "ajusted.pvalue_Genotypes_specific_CELL_COUNTS"),
                                             which(colnames(DE_BLOCK) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))])
        
        DE_BLOCK_subset$VAR<-SELECTED_VARS_sel
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_BLOCK_subset_0\n")
          cat(str(DE_BLOCK_subset))
          cat("\n")
        }
        
        list_Block[[i]]<-DE_BLOCK_subset
        
      }
      
      # cat("DE_BLOCK_FINAL\n")
      
      
      ### CIS gene
      setwd(path6)
      
      filename="DE_CIS_gene.rds"
      
      
      
      if (file.exists(filename)) {
        DE_CIS_gene = readRDS(filename)
        
        # cat("DE_CIS_gene\n")
        # cat(str(DE_CIS_gene))
        # cat("\n")
        
        DE_CIS_gene_subset<-unique(DE_CIS_gene[,-c(which(colnames(DE_CIS_gene) == "ajusted.pvalue_Genotypes_specific_CELL_COUNTS"),
                                                   which(colnames(DE_CIS_gene) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))])
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_CIS_gene_0\n")
          cat(str(DE_CIS_gene))
          cat("\n")
        }
        
        DE_CIS_gene_subset<-unique(DE_CIS_gene[,-c(which(colnames(DE_CIS_gene) == "ajusted.pvalue_Genotypes_specific_CELL_COUNTS"),
                                                   which(colnames(DE_CIS_gene) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))])
        
        DE_CIS_gene_subset$VAR<-SELECTED_VARS_sel
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_CIS_gene_subset_0\n")
          cat(str(DE_CIS_gene_subset))
          cat("\n")
        }
        
        list_CIS[[i]]<-DE_CIS_gene_subset
        
        
        
        #  quit(status = 1)
        
        
      }
      
      
    }# i in 1:length(SELECTED_VARS_UPDATED)
    
    
    Condition_DEBUG <- 1
    
    
    
    # list_DEF[[3]]<-Results_CIS_subset 
    
    
    
    
    Results_DEF<-data.frame()
    
    if(length(list_CIS) >0)
    {
      
      Results_CIS_nominal = as.data.frame(data.table::rbindlist(list_CIS, fill=T), stringsAsFactors=F)
      
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_CIS_nominal_0\n")
        cat(str(Results_CIS_nominal))
        cat("\n")
        cat(str(unique(Results_CIS_nominal$VAR)))
        cat("\n")
      }
      
      
      indx.select<-c(which(colnames(Results_CIS_nominal)== "VAR"),
                     which(colnames(Results_CIS_nominal)== "ensembl_gene_id"),which(colnames(Results_CIS_nominal)== "HGNC"),
                     which(colnames(Results_CIS_nominal)== "coefficient_Genotypes_specific_CELL_COUNTS"),
                     which(colnames(Results_CIS_nominal)== "pvalue_Genotypes_specific_CELL_COUNTS"),
                     which(colnames(Results_CIS_nominal)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_CIS_nominal[,indx.select])
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_0\n")
        cat(str(MT_set))
        cat("\n")
        cat(str(unique(MT_set$ensembl_gene_id)))
        cat("\n")
      }
      
      MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Genotypes_specific_CELL_COUNTS, method = "BH")
      MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS),2)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_1\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      
      MT_set.dt<-data.table(MT_set, key="VAR")
      
      MT_set_collapsed<-as.data.frame(MT_set.dt[,.N,by=key(MT_set.dt)], stringsAsFactors=F)
      
      colnames(MT_set_collapsed)<-c("VAR","n_tested_CIS")
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_collapsed_0\n")
        cat(str(MT_set_collapsed))
        cat("\n")
      }
      
      
      
      DE_CIS_gene<-MT_set[which(MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_CIS_gene_1\n")
        cat(str(DE_CIS_gene))
        cat("\n")
      }
      
      
      DE_CIS_gene.dt<-data.table(DE_CIS_gene, key="VAR")
      
      DE_CIS_gene_collapsed<-as.data.frame(DE_CIS_gene.dt[,.N,by=key(DE_CIS_gene.dt)], stringsAsFactors=F)
      
      colnames(DE_CIS_gene_collapsed)<-c("VAR","n_SIG_CIS")
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_CIS_gene_collapsed_0\n")
        cat(str(DE_CIS_gene_collapsed))
        cat("\n")
      }
      
      MT_set<-merge(MT_set,
                    MT_set_collapsed,
                    by="VAR",
                    all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_2\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      MT_set<-merge(MT_set,
                    DE_CIS_gene_collapsed,
                    by="VAR",
                    all=T)
      
      MT_set$n_SIG_CIS[is.na(MT_set$n_SIG_CIS)]<-0
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_3\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      MT_set$Perc_SIG_CIS<-round(100*(MT_set$n_SIG_CIS/MT_set$n_tested_CIS),1)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_4\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      indx.int<-c(which(colnames(MT_set) == "VAR"),which(colnames(MT_set) == "ensembl_gene_id"),which(colnames(MT_set) == "HGNC"),
                  which(colnames(MT_set) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"),which(colnames(MT_set) == "coefficient_Genotypes_specific_CELL_COUNTS"),
                  which(colnames(MT_set) == "n_breakdown_string"),
                  which(colnames(MT_set) == "n_SIG_CIS"),which(colnames(MT_set) == "n_tested_CIS"),which(colnames(MT_set) == "Perc_SIG_CIS"))
      
      MT_set_subset<-unique(MT_set[,indx.int])
      
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS")]<-"CIS_gene_minuslogpvalue"
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "coefficient_Genotypes_specific_CELL_COUNTS")]<-"CIS_Beta"
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "n_breakdown_string")]<-"CIS_n_breakdown_string"
      
      
      Results_DEF<-rbind(MT_set_subset,Results_DEF)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_0\n")
        cat(str(Results_DEF))
        cat("\n")
      }
      
      # ##################################################################
      # quit(status = 1)
      
    }# length(list_CIS) >0
    
    
    if(length(list_Block) >0)
    {
      
      Results_Block_nominal = as.data.frame(data.table::rbindlist(list_Block, fill=T), stringsAsFactors=F)
      
      
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_Block_nominal_0\n")
        cat(str(Results_Block_nominal))
        cat("\n")
        cat(str(unique(Results_Block_nominal$VAR)))
        cat("\n")
      }
      
      
      indx.select<-c(which(colnames(Results_Block_nominal)== "VAR"),
                     which(colnames(Results_Block_nominal)== "ensembl_gene_id"),which(colnames(Results_Block_nominal)== "HGNC"),
                     which(colnames(Results_Block_nominal)== "coefficient_Genotypes_specific_CELL_COUNTS"),
                     which(colnames(Results_Block_nominal)== "pvalue_Genotypes_specific_CELL_COUNTS"),
                     which(colnames(Results_Block_nominal)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_Block_nominal[,indx.select])
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_0\n")
        cat(str(MT_set))
        cat("\n")
        cat(str(unique(MT_set$VAR)))
        cat("\n")
        cat(str(unique(MT_set$ensembl_gene_id)))
        cat("\n")
      }
      
      MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Genotypes_specific_CELL_COUNTS, method = "BH")
      MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS),2)
      
      MT_set.dt<-data.table(MT_set, key="VAR")
      
      
      MT_set_collapsed<-as.data.frame(MT_set.dt[,.N,by=key(MT_set.dt)], stringsAsFactors=F)
      
      colnames(MT_set_collapsed)<-c("VAR","n_tested_Block")
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_collapsed_0\n")
        cat(str(MT_set_collapsed))
        cat("\n")
      }
      
      
      
      DE_Block<-MT_set[which(MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_Block_1\n")
        cat(str(DE_Block))
        cat("\n")
      }
      
      
      DE_Block.dt<-data.table(DE_Block, key="VAR")
      
      DE_Block_collapsed<-as.data.frame(DE_Block.dt[,.N,by=key(DE_Block.dt)], stringsAsFactors=F)
      
      colnames(DE_Block_collapsed)<-c("VAR","n_SIG_Block")
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_Block_collapsed_0\n")
        cat(str(DE_Block_collapsed))
        cat("\n")
      }
      
      MT_set<-merge(MT_set,
                    MT_set_collapsed,
                    by="VAR",
                    all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_2\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      MT_set<-merge(MT_set,
                    DE_Block_collapsed,
                    by="VAR",
                    all=T)
      
      MT_set$n_SIG_Block[is.na(MT_set$n_SIG_Block)]<-0
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_3\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      MT_set$Perc_SIG_Block<-round(100*(MT_set$n_SIG_Block/MT_set$n_tested_Block),1)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_4\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      indx.int<-c(which(colnames(MT_set) == "VAR"),which(colnames(MT_set) == "ensembl_gene_id"),which(colnames(MT_set) == "HGNC"),
                  which(colnames(MT_set) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"),which(colnames(MT_set) == "coefficient_Genotypes_specific_CELL_COUNTS"),
                  which(colnames(MT_set) == "n_breakdown_string"),
                  which(colnames(MT_set) == "n_SIG_Block"),which(colnames(MT_set) == "n_tested_Block"),which(colnames(MT_set) == "Perc_SIG_Block"))
      
      MT_set_subset<-unique(MT_set[,indx.int])
      
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS")]<-"Block_PCHiC_minuslogpvalue"
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "coefficient_Genotypes_specific_CELL_COUNTS")]<-"Block_PCHiC_Beta"
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "n_breakdown_string")]<-"Block_PCHiC_n_breakdown_string"
      
      
      
      Results_DEF<-merge(Results_DEF,
                         MT_set_subset,
                         by=c("VAR","ensembl_gene_id","HGNC"),
                         all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_Block\n")
        cat(str(Results_DEF))
        cat("\n")
      }
      
      
      
    }# length(list_Block) >0
    
    
    
    if(length(list_GW) >0)
    {
      
      Results_GW_nominal = as.data.frame(data.table::rbindlist(list_GW, fill=T), stringsAsFactors=F)
      
      
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_GW_nominal_0\n")
        cat(str(Results_GW_nominal))
        cat("\n")
        cat(str(unique(Results_GW_nominal$VAR)))
        cat("\n")
      }
      
      
      indx.select<-c(which(colnames(Results_GW_nominal)== "VAR"),
                     which(colnames(Results_GW_nominal)== "ensembl_gene_id"),which(colnames(Results_GW_nominal)== "HGNC"),
                     which(colnames(Results_GW_nominal)== "coefficient_Genotypes_specific_CELL_COUNTS"),
                     which(colnames(Results_GW_nominal)== "pvalue_Genotypes_specific_CELL_COUNTS"),
                     which(colnames(Results_GW_nominal)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_GW_nominal[,indx.select])
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_0\n")
        cat(str(MT_set))
        cat("\n")
        cat(str(unique(MT_set$VAR)))
        cat("\n")
        cat(str(unique(MT_set$ensembl_gene_id)))
        cat("\n")
      }
      
      MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Genotypes_specific_CELL_COUNTS, method = "BH")
      MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS),2)
      
      MT_set.dt<-data.table(MT_set, key="VAR")
      
      
      MT_set_collapsed<-as.data.frame(MT_set.dt[,.N,by=key(MT_set.dt)], stringsAsFactors=F)
      
      colnames(MT_set_collapsed)<-c("VAR","n_tested_GW")
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_collapsed_0\n")
        cat(str(MT_set_collapsed))
        cat("\n")
      }
      
      
      
      DE_GW<-MT_set[which(MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_GW_1\n")
        cat(str(DE_GW))
        cat("\n")
      }
      
      
      DE_GW.dt<-data.table(DE_GW, key="VAR")
      
      DE_GW_collapsed<-as.data.frame(DE_GW.dt[,.N,by=key(DE_GW.dt)], stringsAsFactors=F)
      
      colnames(DE_GW_collapsed)<-c("VAR","n_SIG_GW")
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_GW_collapsed_0\n")
        cat(str(DE_GW_collapsed))
        cat("\n")
      }
      
      MT_set<-merge(MT_set,
                    MT_set_collapsed,
                    by="VAR",
                    all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_2\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      MT_set<-merge(MT_set,
                    DE_GW_collapsed,
                    by="VAR",
                    all=T)
      
      MT_set$n_SIG_GW[is.na(MT_set$n_SIG_GW)]<-0
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_3\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      MT_set$Perc_SIG_GW<-round(100*(MT_set$n_SIG_GW/MT_set$n_tested_GW),1)
      
      if(Condition_DEBUG == 1)
      {
        cat("MT_set_4\n")
        cat(str(MT_set))
        cat("\n")
      }
      
      indx.int<-c(which(colnames(MT_set) == "VAR"),which(colnames(MT_set) == "ensembl_gene_id"),which(colnames(MT_set) == "HGNC"),
                  which(colnames(MT_set) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"),which(colnames(MT_set) == "coefficient_Genotypes_specific_CELL_COUNTS"),
                  which(colnames(MT_set) == "n_breakdown_string"),
                  which(colnames(MT_set) == "n_SIG_GW"),which(colnames(MT_set) == "n_tested_GW"),which(colnames(MT_set) == "Perc_SIG_GW"))
      
      MT_set_subset<-unique(MT_set[,indx.int])
      
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS")]<-"Genome_wide_minuslogpvalue"
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "coefficient_Genotypes_specific_CELL_COUNTS")]<-"Genome_wide_Beta"
      colnames(MT_set_subset)[which(colnames(MT_set_subset) == "n_breakdown_string")]<-"Genome_wide_n_breakdown_string"
      
      
      Results_DEF<-merge(Results_DEF,
                         MT_set_subset,
                         by=c("VAR","ensembl_gene_id","HGNC"),
                         all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_GW\n")
        cat(str(Results_DEF))
        cat("\n")
      }
      
      # quit(status = 1)
      
    }# length(list_GW) >0
    
    if(dim(Results_DEF)[1] >0)
    {
      
      cat("Results_DEF_BORN\n")
      cat(str(Results_DEF))
      cat("\n")
      cat(str(unique(Results_DEF$VAR)))
      cat("\n")
      
      SIG_CIS<-Results_DEF[which(Results_DEF$CIS_gene_minuslogpvalue >= 1.3),]
      
      cat("SIG_CIS_0\n")
      cat(str(SIG_CIS))
      cat("\n")
      cat(str(unique(SIG_CIS$VAR)))
      cat("\n")
      
      SIG_CIS.dt<-data.table(SIG_CIS, key=c("VAR"))
      
      SIG_CIS_collapsed<-data.frame(SIG_CIS.dt[,.(strings_ENSG_CIS=paste(ensembl_gene_id, collapse = ";"),
                                                  strings_HGNC_CIS=paste(HGNC, collapse = ";"),
                                                  strings_CIS_gene_minuslogpvalue=paste(CIS_gene_minuslogpvalue, collapse = ";"),
                                                  n_SIG_CIS=n_SIG_CIS,
                                                  n_tested_CIS=n_tested_CIS,
                                                  Perc_SIG_CIS=Perc_SIG_CIS),
                                               by=key(SIG_CIS.dt)], stringsAsFactors=F)
      
      
      
      
      cat("SIG_CIS_collapsed_0\n")
      cat(str(SIG_CIS_collapsed))
      cat("\n")
      cat(str(unique(SIG_CIS_collapsed$VAR)))
      cat("\n")
      
      SIG_BLOCK<-Results_DEF[which(Results_DEF$Block_PCHiC_minuslogpvalue >= 1.3),]
      
      cat("SIG_BLOCK_0\n")
      cat(str(SIG_BLOCK))
      cat("\n")
      cat(str(unique(SIG_BLOCK$VAR)))
      cat("\n")
      
      SIG_BLOCK.dt<-data.table(SIG_BLOCK, key=c("VAR"))
      
      SIG_BLOCK_collapsed<-data.frame(SIG_BLOCK.dt[,.(strings_ENSG_Block=paste(ensembl_gene_id, collapse = ";"),
                                                      strings_HGNC_Block=paste(HGNC, collapse = ";"),
                                                      strings_Block_minuslogpvalue=paste(Block_PCHiC_minuslogpvalue, collapse = ";"),
                                                      n_SIG_Block=n_SIG_Block,
                                                      n_tested_Block=n_tested_Block,
                                                      Perc_SIG_Block=Perc_SIG_Block),
                                                   by=key(SIG_BLOCK.dt)], stringsAsFactors=F)
      
      
      
      cat("SIG_BLOCK_collapsed_0\n")
      cat(str(SIG_BLOCK_collapsed))
      cat("\n")
      cat(str(unique(SIG_BLOCK_collapsed$VAR)))
      cat("\n")
      
      SIG_GW<-Results_DEF[which(Results_DEF$Genome_wide_minuslogpvalue >= 1.3),]
      
      cat("SIG_GW_0\n")
      cat(str(SIG_GW))
      cat("\n")
      cat(str(unique(SIG_GW$VAR)))
      cat("\n")
      
      SIG_GW.dt<-data.table(SIG_GW, key=c("VAR"))
      
      SIG_GW_collapsed<-data.frame(SIG_GW.dt[,.(strings_ENSG_GW=paste(ensembl_gene_id, collapse = ";"),
                                                strings_HGNC_GW=paste(HGNC, collapse = ";"),
                                                strings_Genome_wide_minuslogpvalue=paste(Genome_wide_minuslogpvalue, collapse = ";"),
                                                n_SIG_GW=n_SIG_GW,
                                                n_tested_GW=n_tested_GW,
                                                Perc_SIG_GW=Perc_SIG_GW),
                                             by=key(SIG_GW.dt)], stringsAsFactors=F)
      
      
      cat("SIG_GW_collapsed_0\n")
      cat(str(SIG_GW_collapsed))
      cat("\n")
      cat(str(unique(SIG_GW_collapsed$VAR)))
      cat("\n")
      
      DEF_SIG<-unique(merge(SIG_CIS_collapsed,
                            SIG_BLOCK_collapsed,
                            by=c("VAR"),
                            all=T))
      
      #"strings_ENSG","strings_HGNC"
      
      DEF_SIG<-unique(merge(DEF_SIG,
                            SIG_GW_collapsed,
                            by=c("VAR"),
                            all=T))
      
      # ,"strings_ENSG","strings_HGNC"
      
      cat("DEF_SIG_0\n")
      cat(str(DEF_SIG))
      cat("\n")
      cat(str(unique(DEF_SIG$VAR)))
      cat("\n")
      
      path6<-paste(out,'FINAL_RESULTS','/', sep='')
      
      # cat("path6\n")
      # cat(sprintf(as.character(path6)))
      # cat("\n")
      
      
      if (file.exists(path6)){
        
        
        
        
      } else {
        dir.create(file.path(path6))
        
      }
      
      setwd(path6)
      
      cat("Results_DEF_FINAL_PRE_Z\n")
      cat(str(Results_DEF))
      cat("\n")
      cat(str(unique(Results_DEF$VAR)))
      cat("\n")
      
      ### Z-score the effect sizes ----
      
      Z_score_effect_sizes_BLOCK<-c(Results_DEF$CIS_Beta[!is.na(Results_DEF$CIS_Beta)],
                              Results_DEF$Block_PCHiC_Beta[!is.na(Results_DEF$Block_PCHiC_Beta)])
      
      cat("Z_score_effect_sizes_BLOCK\n")
      cat(str(Z_score_effect_sizes_BLOCK))
      cat("\n")
      
      mean_Z_score_effect_sizes_BLOCK<-mean(Z_score_effect_sizes_BLOCK)
      
      cat("mean_Z_score_effect_sizes_BLOCK\n")
      cat(str(mean_Z_score_effect_sizes_BLOCK))
      cat("\n")
      
      sd_Z_score_effect_sizes_BLOCK<-sd(Z_score_effect_sizes_BLOCK)
      
      cat("sd_Z_score_effect_sizes_BLOCK\n")
      cat(str(sd_Z_score_effect_sizes_BLOCK))
      cat("\n")
      
      Results_DEF$CIS_Beta_Z_score<-(Results_DEF$CIS_Beta-mean_Z_score_effect_sizes_BLOCK)/sd_Z_score_effect_sizes_BLOCK
      Results_DEF$Block_PCHiC_Beta_Z_score<-(Results_DEF$Block_PCHiC_Beta-mean_Z_score_effect_sizes_BLOCK)/sd_Z_score_effect_sizes_BLOCK
      
      Z_score_effect_sizes_GW<-c(Results_DEF$Genome_wide_Beta[!is.na(Results_DEF$Genome_wide_Beta)])
      
      cat("Z_score_effect_sizes_GW\n")
      cat(str(Z_score_effect_sizes_GW))
      cat("\n")
      
      mean_Z_score_effect_sizes_GW<-mean(Z_score_effect_sizes_GW)
      
      cat("mean_Z_score_effect_sizes_GW\n")
      cat(str(mean_Z_score_effect_sizes_GW))
      cat("\n")
      
      sd_Z_score_effect_sizes_GW<-sd(Z_score_effect_sizes_GW)
      
      cat("sd_Z_score_effect_sizes_GW\n")
      cat(str(sd_Z_score_effect_sizes_GW))
      cat("\n")
      
      Results_DEF$Genome_wide_Beta_Z_score<-(Results_DEF$Genome_wide_Beta-mean_Z_score_effect_sizes_GW)/sd_Z_score_effect_sizes_GW
      
      
      cat("Results_DEF_FINAL_POST_Z\n")
      cat(str(Results_DEF))
      cat("\n")
      cat(str(unique(Results_DEF$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Results_DEF$CIS_Beta_Z_score)))))
      cat("\n")
      cat(sprintf(as.character(summary(Results_DEF$CIS_Beta_Z_score))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Results_DEF$Block_PCHiC_Beta_Z_score)))))
      cat("\n")
      cat(sprintf(as.character(summary(Results_DEF$Block_PCHiC_Beta_Z_score))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Results_DEF$Genome_wide_Beta_Z_score)))))
      cat("\n")
      cat(sprintf(as.character(summary(Results_DEF$Genome_wide_Beta_Z_score))))
      cat("\n")
      # ###############################
      # quit(status = 1)
      # 
      
      write.table(file="ALL_BY_ALL_DE_LM_FPKM_results_Main_VARS.tsv", Results_DEF, sep="\t", row.names = F,quote = F)
      write.table(file="ALL_BY_ALL_DE_LM_FPKM_results_SIG_collapsed_Main_VARS.tsv", DEF_SIG, sep="\t", row.names = F,quote = F)
      
    
    }# dim(Results_DEF)[1] >0
    
  }#length(SELECTED_VARS_UPDATED) >0
  
  
  
  # ##################################
  # quit(status = 1)
  
  
 
  
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
    make_option(c("--SELECTED_VARS_UPDATED"), type="character", default=NULL,
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
  
 
  collate_results_function(opt)
 
}


###########################################################################

system.time( main() )
