

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



GSEA = function (option_list)
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
  library("enrichplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  Condition_DEBUG <- 1
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    SELECTED_VARS_sel<-SELECTED_VARS_UPDATED
    
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
      
      geneList = as.numeric(GSEA_FILE[,which(colnames(GSEA_FILE)== "FC")])
      
      # geneList = GSEA_FILE[,c(which(colnames(GSEA_FILE)== "ENTREZID_id"),
      #                           which(colnames(GSEA_FILE)== "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))]
      
      cat("geneList_0\n")
      cat(str(geneList))
      cat("\n")
      
      # ## feature 2: named vector
      names(geneList) = as.character(GSEA_FILE[,which(colnames(GSEA_FILE)== "ENTREZID_id")])
      
      # ## feature 3: decreasing orde
      geneList = sort(geneList, decreasing = TRUE)
      
      cat("geneList_1\n")
      cat(str(geneList))
      cat("\n")
      
      GSEA<-gseGO(
        geneList,
        OrgDb='org.Hs.eg.db',
        ont = "ALL",
        keyType = "ENTREZID",
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea"
      )
      
      #eps = 1e-10,
      
      cat("GSEA_0\n")
      cat(str(GSEA))
      cat("\n")
      
      setwd(path6)
      
      saveRDS(GSEA,file="GSEA_result.rds")
      
    }#file.exists("DE_RESULTS_FOR_GSEA.rds")
  }#length(SELECTED_VARS_UPDATED) >0
}

GSEA_plot = function (option_list)
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
  library("enrichplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  library("DOSE", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
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
    SELECTED_VARS_sel<-SELECTED_VARS_UPDATED
    
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
      
      geneList = as.numeric(GSEA_FILE[,which(colnames(GSEA_FILE)== "FC")])
      
      # geneList = GSEA_FILE[,c(which(colnames(GSEA_FILE)== "ENTREZID_id"),
      #                           which(colnames(GSEA_FILE)== "ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"))]
      
      cat("geneList_0\n")
      cat(str(geneList))
      cat("\n")
      
      # ## feature 2: named vector
      names(geneList) = as.character(GSEA_FILE[,which(colnames(GSEA_FILE)== "ENTREZID_id")])
      
      # ## feature 3: decreasing orde
      geneList = sort(geneList, decreasing = TRUE)
      
      cat("geneList_1\n")
      cat(str(geneList))
      cat("\n")
    
      setwd(path6)
      
      GSEA_result<-readRDS(GSEA,file="GSEA_result.rds")
      
      if(Condition_DEBUG == 1)
      {
        cat("GSEA_result_0\n")
        cat(str(GSEA_result))
        cat("\n")
      }
    
      
      array_GO_ID<-unique(GSEA_result@result$ID)  
      
      if(Condition_DEBUG == 1)
      {
        cat("array_GO_ID_0\n")
        cat(str(array_GO_ID))
        cat("\n")
      }
    
      for(i in 1:length(array_GO_ID))
      {
        array_GO_ID_sel<-array_GO_ID[i]
        
        cat("---->\t")
        cat(sprintf(as.character(array_GO_ID_sel)))
        cat("\n")
        
        GSEA_result_sel<-GSEA_result[which(GSEA_result@result$ID == array_GO_ID_sel),]
        
        if(Condition_DEBUG == 1)
        {
          cat("GSEA_result_sel_0\n")
          cat(str(GSEA_result_sel))
          cat("\n")
        }
        

        Description_sel<-gsub("\\s+","_",unique(GSEA_result_sel$Description))
        Description_sel<-gsub(",","_",Description_sel)
        Description_sel<-gsub("-","_",Description_sel)
        Description_sel<-gsub("\\/","_",Description_sel)
        
        
        cat("---->\t")
        cat(sprintf(as.character(Description_sel)))
        cat("\t")
        
        minus_logpval_adjusted_GO_sel<-unique(round(-1*log10(GSEA_result_sel$p.adjust),2))
        
        cat("---->\t")
        cat(sprintf(as.character(minus_logpval_adjusted_GO_sel)))
        cat("\t")
          
        p3 <- gseaplot(GSEA_result, geneSetID = i, title = GSEA_result$Description[i])
        
        path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
        
        # cat("path7\n")
        # cat(sprintf(as.character(path7)))
        # cat("\n")
        
        
        if (file.exists(path7)){
          
          
          
          
        } else {
          dir.create(file.path(path7))
          
        }
        
        
        
        
        path8<-paste(path7,SELECTED_VARS_sel,'/', sep='')
        
        if (file.exists(path8)){
          
          
          
          
        } else {
          dir.create(file.path(path8))
          
        }
        
        path9<-paste(path8,'GSEA_Leading_Edge','/', sep='')
        
        if (file.exists(path9)){
          
          
          
          
        } else {
          dir.create(file.path(path9))
          
        }
        
        setwd(path9)
        
        array_GO_ID_sel<-gsub(":","_",array_GO_ID_sel)
        
        svgname<-paste("Leading_Edge_",minus_logpval_adjusted_GO_sel,"_",Description_sel,"_",array_GO_ID_sel,"_",SELECTED_VARS_sel,".svg",sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= p3,
                 device="svg",
                 height=10, width=12)
        }
        
        if(Condition_DEBUG == 1)
        {
          cat("PRINTING#1\n")
          
        }
        
        # quit(status = 1)
      }# i in 1:length(array_GO_ID)
      
      
      
     
      
      # ######################### HERE HERE
      # quit(status = 1)
      
    }#file.exists("DE_RESULTS_FOR_GSEA.rds")
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
    make_option(c("--ALL_dB"), type="character", default=NULL,
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
  
  
  # GSEA(opt)
  GSEA_plot(opt)
  
}

# echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
# echo "--PCHiC_info $PCHiC_info \\" >> $output
# echo "--ALL_dB $ALL_dB \\" >> $output


###########################################################################

system.time( main() )
