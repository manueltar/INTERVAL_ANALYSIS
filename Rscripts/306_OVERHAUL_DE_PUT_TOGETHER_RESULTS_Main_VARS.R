

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
    
    for(i in 1:length(SELECTED_VARS_UPDATED ))
    {
      
      SELECTED_VARS_sel<-SELECTED_VARS_UPDATED [i]
      
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
      
      n_SIG_GW<-0
      n_tested_GW<-0
      Perc_SIG_GW<-NA
      
      if (file.exists(filename)) {
        DE_genome_wide = readRDS(filename)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_genome_wide_0\n")
          cat(str(DE_genome_wide))
          cat("\n")
        }
        
        DE_genome_wide_SIG<-DE_genome_wide[which(DE_genome_wide$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_genome_wide_SIG_0\n")
          cat(str(DE_genome_wide_SIG))
          cat("\n")
        }
        
        n_tested_GW<-length(unique(DE_genome_wide$ensembl_gene_id))
        n_SIG_GW<-length(unique(DE_genome_wide_SIG$ensembl_gene_id))
        Perc_SIG_GW<-round(100*(n_SIG_GW/n_tested_GW),2)
        
        DE_genome_wide_SIG$n_tested_GW<-n_tested_GW
        DE_genome_wide_SIG$n_SIG_GW<-n_SIG_GW
        DE_genome_wide_SIG$Perc_SIG_GW<-Perc_SIG_GW
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_genome_wide_SIG_1\n")
          cat(str(DE_genome_wide_SIG))
          cat("\n")
        }
        
        if(dim(DE_genome_wide_SIG)[1] >0)
        {
          DE_genome_wide_SIG$VAR<-SELECTED_VARS_sel
          
          if(Condition_DEBUG == 1)
          {
            cat("DE_genome_wide_SIG_2\n")
            cat(str(DE_genome_wide_SIG))
            cat("\n")
          }
          
          list_GW[[i]]<-DE_genome_wide_SIG
          
          # list_Block<-list()
          # list_CIS<-list()
          
        }#dim(DE_genome_wide_SIG)[1] >0
      }#file.exists(DE_genome_wide.rds)
      
      ### DTU block
      
      setwd(path6)
      
      filename="DE_Block.rds"
      
      n_SIG_Block<-0
      n_tested_Block<-0
      Perc_SIG_Block<-NA
      
      if (file.exists(filename)) {
        DE_BLOCK = readRDS(filename)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_BLOCK_INTIAL\n")
          cat(str(DE_BLOCK))
          cat("\n")
        }
        
        DE_BLOCK_SIG<-DE_BLOCK[which(DE_BLOCK$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_BLOCK_SIG_0\n")
          cat(str(DE_BLOCK_SIG))
          cat("\n")
        }
        
        n_tested_Block<-length(unique(DE_BLOCK$ensembl_gene_id))
        n_SIG_Block<-length(unique(DE_BLOCK_SIG$ensembl_gene_id))
        Perc_SIG_Block<-round(100*(n_SIG_Block/n_tested_Block),2)
        
        DE_BLOCK$n_tested_Block<-n_tested_Block
        DE_BLOCK$n_SIG_Block<-n_SIG_Block
        DE_BLOCK$Perc_SIG_Block<-Perc_SIG_Block
        
        
        # cat(str(DE_BLOCK))
        # cat("\n")
        
        if(dim(DE_BLOCK)[1] >0)
        {
          DE_BLOCK$VAR<-SELECTED_VARS_sel
          
          if(Condition_DEBUG == 1)
          {
            cat("DE_BLOCK_2\n")
            cat(str(DE_BLOCK))
            cat("\n")
          }
                    
          list_Block[[i]]<-DE_BLOCK
          
          
          # list_CIS<-list()
          
        }
        
        
      }
      
      # cat("DE_BLOCK_FINAL\n")
      
      
      ### CIS gene
      setwd(path6)
      
      filename="DE_CIS_gene.rds"
      
      n_SIG_CIS<-0
      n_tested_CIS<-0
      Perc_SIG_CIS<-NA
      
      if (file.exists(filename)) {
        DE_CIS_gene = readRDS(filename)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_CIS_gene\n")
          cat(str(DE_CIS_gene))
          cat("\n")
        }
        
        DE_CIS_gene_SIG<-DE_CIS_gene[which(DE_CIS_gene$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_CIS_gene_SIG_0\n")
          cat(str(DE_CIS_gene_SIG))
          cat("\n")
        }
        
        n_tested_CIS<-length(unique(DE_CIS_gene$ensembl_gene_id))
        n_SIG_CIS<-length(unique(DE_CIS_gene_SIG$ensembl_gene_id))
        Perc_SIG_CIS<-round(100*(n_SIG_CIS/n_tested_CIS),2)
        
        DE_CIS_gene$n_tested_CIS<-n_tested_CIS
        DE_CIS_gene$n_SIG_CIS<-n_SIG_CIS
        DE_CIS_gene$Perc_SIG_CIS<-Perc_SIG_CIS
        
        
        # cat(str(DE_CIS_gene))
        # cat("\n")
        
        if(dim(DE_CIS_gene)[1] >0)
        {
          DE_CIS_gene$VAR<-SELECTED_VARS_sel
          
          if(Condition_DEBUG == 1)
          {
            cat("DE_CIS_gene_2\n")
            cat(str(DE_CIS_gene))
            cat("\n")
          }
                    
          list_CIS[[i]]<-DE_CIS_gene
          
          
          # list_CIS<-list()
          
        }
        
        #  quit(status = 1)
        
        
      }
    }# i in 1:length(SELECTED_VARS_UPDATED )
    
    
    Condition_DEBUG <- 1
    
    Results_DEF<-data.frame()
    
    # list_DEF[[3]]<-Results_CIS_subset 
    
    if(length(list_CIS) >0)
    {
      
      Results_CIS = as.data.frame(data.table::rbindlist(list_CIS, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_CIS_0\n")
        cat(str(Results_CIS))
        cat("\n")
        cat(str(unique(Results_CIS$VAR)))
        cat("\n")
      }
      
      
      Results_CIS_SIG<-Results_CIS[which(Results_CIS$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_CIS_SIG_0\n")
        cat(str(Results_CIS_SIG))
        cat("\n")
        cat(str(unique(Results_CIS_SIG$VAR)))
        cat("\n")
      }
      
      
      
      indx.int<-c(which(colnames(Results_CIS) == "VAR"),which(colnames(Results_CIS) == "ensembl_gene_id"),which(colnames(Results_CIS) == "HGNC"),
                  which(colnames(Results_CIS) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"),which(colnames(Results_CIS) == "coefficient_Genotypes_specific_CELL_COUNTS"),
                  which(colnames(Results_CIS) == "n_SIG_CIS"),which(colnames(Results_CIS) == "n_tested_CIS"),which(colnames(Results_CIS) == "Perc_SIG_CIS"))
      
      Results_CIS_subset<-unique(Results_CIS[,indx.int])
      
      colnames(Results_CIS_subset)[which(colnames(Results_CIS_subset) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS")]<-"CIS_gene_minuslogpvalue"
      colnames(Results_CIS_subset)[which(colnames(Results_CIS_subset) == "coefficient_Genotypes_specific_CELL_COUNTS")]<-"CIS_Beta"
      
      Results_DEF<-rbind(Results_CIS_subset,Results_DEF)
      
      
      
      
    }# length(list_CIS) >0
    
    if(length(list_Block) >0)
    {
      
      Results_Block = as.data.frame(data.table::rbindlist(list_Block, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_Block_0\n")
        cat(str(Results_Block))
        cat("\n")
        cat(str(unique(Results_Block$VAR)))
        cat("\n")
      }
      
      
      Results_Block_SIG<-Results_Block[which(Results_Block$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_Block_SIG_0\n")
        cat(str(Results_Block_SIG))
        cat("\n")
        cat(str(unique(Results_Block_SIG$VAR)))
        cat("\n")
      }
      
      
      
      
      indx.int<-c(which(colnames(Results_Block) == "VAR"),which(colnames(Results_Block) == "ensembl_gene_id"),which(colnames(Results_Block) == "HGNC"),
                  which(colnames(Results_Block) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"),which(colnames(Results_Block) == "coefficient_Genotypes_specific_CELL_COUNTS"),
                  which(colnames(Results_Block) == "n_SIG_Block"),which(colnames(Results_Block) == "n_tested_Block"),which(colnames(Results_Block) == "Perc_SIG_Block"))
      
      Results_Block_subset<-unique(Results_Block[,indx.int])
      
      colnames(Results_Block_subset)[which(colnames(Results_Block_subset) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS")]<-"Block_PCHiC_minuslogpvalue"
      colnames(Results_Block_subset)[which(colnames(Results_Block_subset) == "coefficient_Genotypes_specific_CELL_COUNTS")]<-"Block_PCHiC_Beta"
      
      
      Results_DEF<-merge(Results_DEF,
                         Results_Block_subset,
                         by=c("VAR","ensembl_gene_id","HGNC"),
                         all=T)
      
    }# length(list_Block) >0
    
    if(length(list_GW) >0)
    {
      
      Results_GW = as.data.frame(data.table::rbindlist(list_GW, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_GW_0\n")
        cat(str(Results_GW))
        cat("\n")
        cat(str(unique(Results_GW$VAR)))
        cat("\n")
      }
      
      
      
      indx.int<-c(which(colnames(Results_GW) == "VAR"),which(colnames(Results_GW) == "ensembl_gene_id"),which(colnames(Results_GW) == "HGNC"),
                  which(colnames(Results_GW) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"),which(colnames(Results_GW) == "coefficient_Genotypes_specific_CELL_COUNTS"),
                  which(colnames(Results_GW) == "n_SIG_GW"),which(colnames(Results_GW) == "n_tested_GW"),which(colnames(Results_GW) == "Perc_SIG_GW"))
      
      Results_GW_subset<-unique(Results_GW[,indx.int])
      
      colnames(Results_GW_subset)[which(colnames(Results_GW_subset) == "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS")]<-"Genome_wide_minuslogpvalue"
      colnames(Results_GW_subset)[which(colnames(Results_GW_subset) == "coefficient_Genotypes_specific_CELL_COUNTS")]<-"Genome_wide_Beta"
      
      
      Results_GW_subset_SIG<-Results_GW_subset[which(Results_GW_subset$Genome_wide_minuslogpvalue >= 1.3),]
      
      Results_DEF<-merge(Results_DEF,
                         Results_GW_subset_SIG,
                         by=c("VAR","ensembl_gene_id","HGNC"),
                         all=T)
      
      
      
      
      
    }# length(list_GW) >0
    
    if(dim(Results_DEF)[1] >0)
    {
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_BORN\n")
        cat(str(Results_DEF))
        cat("\n")
        cat(str(unique(Results_DEF$VAR)))
        cat("\n")
      }
      
      # ##################################################################################
      # quit(status = 1)
      
      
      Results_DEF.dt<-data.table(Results_DEF, key=c("VAR","ensembl_gene_id","HGNC"))
      
      Results_DEF_max_GW.per.gene<-data.frame(Results_DEF.dt[,.SD[which.max(Genome_wide_minuslogpvalue)],by=key(Results_DEF.dt)], stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_max_GW.per.gene_0\n")
        cat(str(Results_DEF_max_GW.per.gene))
        cat("\n")
        cat(str(unique(Results_DEF_max_GW.per.gene$VAR)))
        cat("\n")
      }
      
      SIG_GW<-Results_DEF_max_GW.per.gene[which(Results_DEF_max_GW.per.gene$Genome_wide_minuslogpvalue >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("SIG_GW_0\n")
        cat(str(SIG_GW))
        cat("\n")
        cat(str(unique(SIG_GW$VAR)))
        cat("\n")
      }
      
      SIG_GW.dt<-data.table(SIG_GW, key=c("VAR"))
      
      SIG_GW_collapsed<-data.frame(SIG_GW.dt[,.(strings_ENSG_GW=paste(ensembl_gene_id, collapse = ";"),
                                                strings_HGNC_GW=paste(HGNC, collapse = ";"),
                                                strings_Genome_wide_minuslogpvalue=paste(Genome_wide_minuslogpvalue, collapse = ";"),
                                                n_SIG_GW=n_SIG_GW,
                                                n_tested_GW=n_tested_GW,
                                                Perc_SIG_GW=Perc_SIG_GW),
                                             by=key(SIG_GW.dt)], stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("SIG_GW_collapsed_0\n")
        cat(str(SIG_GW_collapsed))
        cat("\n")
        cat(str(unique(SIG_GW_collapsed$VAR)))
        cat("\n")
      }
      
      # quit(status = 1)
      
      
      Results_DEF_max_Block.per.gene<-data.frame(Results_DEF.dt[,.SD[which.max(Block_PCHiC_minuslogpvalue)],by=key(Results_DEF.dt)], stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_max_Block.per.gene_0\n")
        cat(str(Results_DEF_max_Block.per.gene))
        cat("\n")
        cat(str(unique(Results_DEF_max_Block.per.gene$VAR)))
        cat("\n")
      }
      
      SIG_BLOCK<-Results_DEF_max_Block.per.gene[which(Results_DEF_max_Block.per.gene$Block_PCHiC_minuslogpvalue >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("SIG_BLOCK_0\n")
        cat(str(SIG_BLOCK))
        cat("\n")
        cat(str(unique(SIG_BLOCK$VAR)))
        cat("\n")
      }
      
      SIG_BLOCK.dt<-data.table(SIG_BLOCK, key=c("VAR"))
      
      SIG_BLOCK_collapsed<-data.frame(SIG_BLOCK.dt[,.(strings_ENSG_Block=paste(ensembl_gene_id, collapse = ";"),
                                                      strings_HGNC_Block=paste(HGNC, collapse = ";"),
                                                      strings_Block_minuslogpvalue=paste(Block_PCHiC_minuslogpvalue, collapse = ";"),
                                                      n_SIG_Block=n_SIG_Block,
                                                      n_tested_Block=n_tested_Block,
                                                      Perc_SIG_Block=Perc_SIG_Block),
                                                   by=key(SIG_BLOCK.dt)], stringsAsFactors=F)
      
      
      if(Condition_DEBUG == 1)
      {
        cat("SIG_BLOCK_collapsed_0\n")
        cat(str(SIG_BLOCK_collapsed))
        cat("\n")
        cat(str(unique(SIG_BLOCK_collapsed$VAR)))
        cat("\n")
      }
      
      # quit(status = 1)
      
      Results_DEF_max_CIS.per.gene<-data.frame(Results_DEF.dt[,.SD[which.max(CIS_gene_minuslogpvalue)],by=key(Results_DEF.dt)], stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_max_CIS.per.gene_0\n")
        cat(str(Results_DEF_max_CIS.per.gene))
        cat("\n")
        cat(str(unique(Results_DEF_max_CIS.per.gene$VAR)))
        cat("\n")
      }
      
      
      SIG_CIS<-Results_DEF_max_CIS.per.gene[which(Results_DEF_max_CIS.per.gene$CIS_gene_minuslogpvalue >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("SIG_CIS_0\n")
        cat(str(SIG_CIS))
        cat("\n")
        cat(str(unique(SIG_CIS$VAR)))
        cat("\n")
      }
      
      SIG_CIS.dt<-data.table(SIG_CIS, key=c("VAR"))
      
      SIG_CIS_collapsed<-data.frame(SIG_CIS.dt[,.(strings_ENSG_CIS=paste(ensembl_gene_id, collapse = ";"),
                                                  strings_HGNC_CIS=paste(HGNC, collapse = ";"),
                                                  strings_CIS_gene_minuslogpvalue=paste(CIS_gene_minuslogpvalue, collapse = ";"),
                                                  n_SIG_CIS=n_SIG_CIS,
                                                  n_tested_CIS=n_tested_CIS,
                                                  Perc_SIG_CIS=Perc_SIG_CIS),
                                               by=key(SIG_CIS.dt)], stringsAsFactors=F)
      
      
      
      if(Condition_DEBUG == 1)
      {
        cat("SIG_CIS_collapsed_0\n")
        cat(str(SIG_CIS_collapsed))
        cat("\n")
        cat(str(unique(SIG_CIS_collapsed$VAR)))
        cat("\n")
      }
      
      
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
      
      path6<-paste(out,'FINAL_RESULTS','/', sep='')
      
      # cat("path6\n")
      # cat(sprintf(as.character(path6)))
      # cat("\n")
      
      
      if (file.exists(path6)){
        
        
        
        
      } else {
        dir.create(file.path(path6))
        
      }
      
      setwd(path6)
      
      
      write.table(file="DE_LM_FPKM_results_Main_VARS.tsv", Results_DEF, sep="\t", row.names = F,quote = F)
      write.table(file="DE_LM_FPKM_results_SIG_collapsed_Main_VARS.tsv", DEF_SIG, sep="\t", row.names = F,quote = F)
      
      
    }# length(list_DEF) >0
  }#length(SELECTED_VARS_UPDATED) >0
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
    make_option(c("--SELECTED_VARS_UPDATED "), type="character", default=NULL,
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
