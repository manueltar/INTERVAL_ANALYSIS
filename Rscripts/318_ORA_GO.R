

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



enrichGO = function (option_list)
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
  # library("enrichplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
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
  
  #### READ Threshold_logpval ----
  
  Threshold_logpval = opt$Threshold_logpval
  
  cat("Threshold_logpval_\n")
  cat(sprintf(as.character(Threshold_logpval)))
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
      
      #### from https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/ ----
      
      
      # we want the log2 fold change 
      original_gene_list <- GSEA_FILE$Genome_wide_minuslogpvalue
      
      # name the vector
      names(original_gene_list) <- GSEA_FILE$ensembl_gene_id
      
      # omit any NA values 
      geneList<-na.omit(original_gene_list)
      
      # sort the list in decreasing order (required for clusterProfiler)
      geneList = sort(geneList, decreasing = TRUE)
      
      cat("geneList\n")
      cat(str(geneList))
      cat("\n")
      
      # # Exctract significant results ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3
      # 
      # GSEA_FILE_SIG<-GSEA_FILE[which(GSEA_FILE$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      # 
      # cat("GSEA_FILE_SIG_0\n")
      # cat(str(GSEA_FILE_SIG))
      # cat("\n")
      
      gene <- geneList[which(geneList >= Threshold_logpval)]
      
      if(length(gene) >0)
      {
        cat("gene_0\n")
        cat(str(gene))
        cat("\n")
        
        gene <- names(gene)
        
        cat("gene_1\n")
        cat(str(gene))
        cat("\n")
        
        
        
        
        
        cat("EGO_START\n")
        
        ego <- clusterProfiler::enrichGO(gene          = gene,
                                         universe      = names(geneList),
                                         OrgDb         = org.Hs.eg.db,
                                         ont           = "ALL",
                                         pAdjustMethod = "BH",
                                         minGSSize = 10,
                                         maxGSSize = 500,
                                         pvalueCutoff  = 0.01,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)
        
        
        
        #eps = 1e-10,
        
        cat("ego\n")
        cat(str(ego))
        cat("\n")
        
        FLAG_NULL<-sum(is.null(ego))
        
        
        cat("FLAG_NULL\n")
        cat(str(FLAG_NULL))
        cat("\n")
        
        # quit(status = 1)
        
        if(FLAG_NULL == 0)
        {
          setwd(path6)
          
          saveRDS(ego,file=paste("enricher_object_GO_",Threshold_logpval,".rds",sep=''))
          
          
        }#FLAG_NULL == 0
        
      }#length(gene) >0
      
      
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
    make_option(c("--Threshold_logpval"), type="numeric", default=NULL,
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
  
  
  enrichGO(opt)
    
  
}

# echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
# echo "--PCHiC_info $PCHiC_info \\" >> $output
# echo "--ALL_dB $ALL_dB \\" >> $output


###########################################################################

system.time( main() )
