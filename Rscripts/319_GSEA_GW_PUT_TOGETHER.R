

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
  suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GO.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("org.Hs.eg.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("Homo.sapiens", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("gwascat", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("rtracklayer", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("clusterProfiler", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
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
  
  
  
  SELECTED_VARS_UPDATED  = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
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
    list_GW<-list()
    
    for(i in 1:length(SELECTED_VARS_UPDATED))
    {
      
      SELECTED_VARS_sel<-SELECTED_VARS_UPDATED[i]
      
      cat("------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(SELECTED_VARS_sel)))
      cat("\n")
      
      #### Read files for Residuals violin plots ----
      
      path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
      
      
      setwd(path6)
      
      if(file.exists("DE_RESULTS_FOR_GSEA.rds"))
      {
        GSEA_FILE<-readRDS(file="DE_RESULTS_FOR_GSEA.rds")
        
        if(Condition_DEBUG == 1)
        {
          cat("GSEA_FILE_0\n")
          cat(str(GSEA_FILE))
          cat("\n")
          cat(str(unique(GSEA_FILE$VAR)))
          cat("\n")
          cat(str(unique(GSEA_FILE$ensembl_gene_id)))
          cat("\n")
        }
        
        if(file.exists("GSEA_result.rds"))
        {
          RESULT_GSEA<-readRDS(file="GSEA_result.rds")
          
          if(Condition_DEBUG == 1)
          {
            cat("RESULT_GSEA_0\n")
            cat(str(RESULT_GSEA))
            cat("\n")
            
          }
          
          RESULT_GSEA_SIG<-RESULT_GSEA@result
          
          # cat("RESULT_GSEA_SIG_0\n")
          # cat(str(RESULT_GSEA_SIG))
          # cat("\n")
          
          if(dim(RESULT_GSEA_SIG)[1] >0)
          {
            
            array_core_enrichment<-unique(RESULT_GSEA_SIG[,c(which(colnames(RESULT_GSEA_SIG) == "ID"),
                                                             which(colnames(RESULT_GSEA_SIG) == "core_enrichment"))])
            
            # cat("array_core_enrichment\n")
            # cat(str(array_core_enrichment))
            # cat("\n")
            
            array_core_enrichment.dt<-data.table(array_core_enrichment, key="ID")
            
            array_core_enrichment_unlisted<-as.data.frame(array_core_enrichment.dt[,.(ENTREZID_id=unlist(strsplit(core_enrichment, split='/'))),
                                                                                   by=key(array_core_enrichment.dt)], stringsAsFactors=F)
            
            
            # cat("array_core_enrichment_unlisted_0\n")
            # cat(str(array_core_enrichment_unlisted))
            # cat("\n")
            
            
            array_core_enrichment_unlisted$ensembl_gene_id= mapIds(org.Hs.eg.db,
                                                                   keys=array_core_enrichment_unlisted$ENTREZID_id, 
                                                                   column="ENSEMBL",
                                                                   keytype="ENTREZID",
                                                                   multiVals="first")
            
            
            # cat("array_core_enrichment_unlisted_1\n")
            # cat(str(array_core_enrichment_unlisted))
            # cat("\n")
            
            array_core_enrichment_unlisted.dt<-data.table(array_core_enrichment_unlisted, key="ID")
            
            array_core_enrichment_collapsed<-as.data.frame(array_core_enrichment_unlisted.dt[,.(string_ENSG_core=paste(ensembl_gene_id, collapse=";")),
                                                                                             by=key(array_core_enrichment_unlisted.dt)], stringsAsFactors=F)
            
            
            # cat("array_core_enrichment_collapsed_0\n")
            # cat(str(array_core_enrichment_collapsed))
            # cat("\n")
            
            RESULT_GSEA_SIG<-merge(RESULT_GSEA_SIG,
                                   array_core_enrichment_collapsed,
                                   by="ID",
                                   all.x=T)
            
            # cat("RESULT_GSEA_SIG_2\n")
            # cat(str(RESULT_GSEA_SIG))
            # cat("\n")
            
            
            RESULT_GSEA_SIG$VAR<-SELECTED_VARS_sel
            
            # cat("RESULT_GSEA_SIG_3\n")
            # cat(str(RESULT_GSEA_SIG))
            # cat("\n")
            
            # quit(status = 1)
            
            
            list_GW[[i]]<-RESULT_GSEA_SIG
            
            # list_Block<-list()
            # list_CIS<-list()
            
          }#dim(RESULT_GSEA_SIG)[1] >0
        }#file.exists("GSEA_result.rds")
      }#file.exists("DE_RESULTS_FOR_GSEA.rds")
      
      
    }#i in 1:length(SELECTED_VARS)
    
    Condition_DEBUG <- 1
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
      
      
      
      
      path7<-paste(out,'FINAL_RESULTS','/', sep='')
      
      # cat("path7\n")
      # cat(sprintf(as.character(path7)))
      # cat("\n")
      
      
      if (file.exists(path7)){
        
        
        
        
      } else {
        dir.create(file.path(path7))
        
      }
      
      setwd(path7)
      
      write.table(file="GSEA_GW_RESULTS.tsv", Results_GW, sep="\t", row.names = F,quote = F)
      
    }# length(list_GW) >0
    
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
    make_option(c("--SELECTED_VARS"), type="character", default=NULL,
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
