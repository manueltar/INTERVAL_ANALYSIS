

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


enricher = function (option_list)
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
  
  TERM2GENE<-data.frame()
  
  
  ### EPO genes----
  
  EPO_HGNC<-c("EPO","EPOR","GRB2","HRAS","JAK2","PLCG1","PTPN6","SHC1","SOS1","STAT5A","STAT5B")
  EPO_TERM<-rep("EPO_signaling_pathway", length(EPO_HGNC))
  
  EPO.df<-as.data.frame(cbind(EPO_HGNC,EPO_TERM), stringsAsFactors=F)
  
  colnames(EPO.df)<-c("SYMBOL","TERM")
  
  cat("EPO.df_0\n")
  cat(str(EPO.df))
  cat("\n")
  
  
  EPO.df$ensembl_gene_id= mapIds(org.Hs.eg.db,
                                 keys=EPO.df$SYMBOL,
                                 column="ENSEMBL",
                                 keytype="SYMBOL",
                                 multiVals="first")
  
  
  cat("EPO.df_1\n")
  cat(str(EPO.df))
  cat("\n")
  
  ### THPO genes----
  
  THPO_HGNC<- c("CSNK2A1","FOS","GRB2","HRAS","JAK2","JUN","MAP2K1","MAPK3","MPL","PIK3CA","PIK3CG","PIK3R1","PLCG1","PRKCA","PRKCB","RAF1","RASA1","SHC1","SOS1","STAT1","STAT3","STAT5A","STAT5B","THPO")
  THPO_TERM<-rep("THPO_signaling_pathway", length(THPO_HGNC))
  
  THPO.df<-as.data.frame(cbind(THPO_HGNC,THPO_TERM), stringsAsFactors=F)
  
  colnames(THPO.df)<-c("SYMBOL","TERM")
  
  cat("THPO.df_0\n")
  cat(str(THPO.df))
  cat("\n")
  
  
  THPO.df$ensembl_gene_id= mapIds(org.Hs.eg.db,
                                  keys=THPO.df$SYMBOL,
                                  column="ENSEMBL",
                                  keytype="SYMBOL",
                                  multiVals="first")
  
  
  cat("THPO.df_1\n")
  cat(str(THPO.df))
  cat("\n")
  
  ### HAMP genes----
  
  HAMP_HGNC<- c("BMP6","HAMP","HFE","HJV","ID1","SMAD7","TMPRSS6")
  HAMP_TERM<-rep("HAMP_signaling_pathway", length(HAMP_HGNC))
  
  HAMP.df<-as.data.frame(cbind(HAMP_HGNC,HAMP_TERM), stringsAsFactors=F)
  
  colnames(HAMP.df)<-c("SYMBOL","TERM")
  
  cat("HAMP.df_0\n")
  cat(str(HAMP.df))
  cat("\n")
  
  
  HAMP.df$ensembl_gene_id= mapIds(org.Hs.eg.db,
                                  keys=HAMP.df$SYMBOL,
                                  column="ENSEMBL",
                                  keytype="SYMBOL",
                                  multiVals="first")
  
  
  cat("HAMP.df_1\n")
  cat(str(HAMP.df))
  cat("\n")
  
  
  
  ### CEBPE genes----
  
  CEBPE_HGNC<- c("MAPK8","PPP1CB","PPP1R3G","XBP1","ATF6","PPP1CA","EIF2AK2","ERN1","ATF4","DDIT3","MBTPS1","EIF2AK3","PPP1R12A","MBTPS2","PPP1R15A","PPP1R3D","PPP1CC","PPP1R2","PPP1R8","PPP1R1A","PPP1R7","PPP1R3A","PPP1R3F","PPP1R3B","PPP1R14C","PPP1R1C","PPP1R14B","PPP1R14A","PPP1R13B","PPP1R16B","PPP1R10","PPP1R12C","PPP1R3E","PPP1R14D","PPP1R9A","PPP1R3C","MAPK9","MAPK10","PPP1R12B","CEBPB","CEBPA","CEBPD","CEBPZ","CEBPE","EIF2A","PPP1R1B","CEBPG","PPP1R11","PPP1R9B","PPP1R15B","PPP1R16A","BCL2","HSPA5","EIF2AK1")
  
  CEBPE_TERM<-rep("CEBPE_signaling_pathway", length(CEBPE_HGNC))
  
  CEBPE.df<-as.data.frame(cbind(CEBPE_HGNC,CEBPE_TERM), stringsAsFactors=F)
  
  colnames(CEBPE.df)<-c("SYMBOL","TERM")
  
  cat("CEBPE.df_0\n")
  cat(str(CEBPE.df))
  cat("\n")
  
  
  CEBPE.df$ensembl_gene_id= mapIds(org.Hs.eg.db,
                                   keys=CEBPE.df$SYMBOL,
                                   column="ENSEMBL",
                                   keytype="SYMBOL",
                                   multiVals="first")
  
  
  cat("CEBPE.df_1\n")
  cat(str(CEBPE.df))
  cat("\n")
  
  ####  merge ALL ----
  
  TERM2GENE<-rbind(EPO.df,THPO.df,HAMP.df,CEBPE.df)
  
  cat("TERM2GENE_1\n")
  cat(str(TERM2GENE))
  cat("\n")
  
  TERM2GENE$ENTREZID = mapIds(org.Hs.eg.db,
                              keys=TERM2GENE$ensembl_gene_id, 
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
  
  cat("TERM2GENE_1\n")
  cat(str(TERM2GENE))
  cat("\n")
  
  
  TERM2GENE_subset<-TERM2GENE[,c(which(colnames(TERM2GENE) == "TERM"),
                                 which(colnames(TERM2GENE) == "ensembl_gene_id"))]
  
  
  cat("TERM2GENE_subset_1\n")
  cat(str(TERM2GENE_subset))
  cat("\n")
  
  colnames(TERM2GENE_subset)<-c("TERM","ensembl_gene_id")
  
  cat("TERM2GENE_subset_0\n")
  cat(str(TERM2GENE_subset))
  cat("\n")
  
  TERM2GENE_subset<-TERM2GENE_subset[!is.na(TERM2GENE_subset$ensembl_gene_id),]
  
  cat("TERM2GENE_subset_1\n")
  cat(str(TERM2GENE_subset))
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
        
        
        #### enricher ----
        
        cat("EGO_START\n")
        
        enricher_object <- clusterProfiler::enricher(gene          = gene,
                                                     universe      = names(geneList),
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = "BH",
                                                     minGSSize = 10,
                                                     maxGSSize = 500,
                                                     qvalueCutoff = 0.2,
                                                     TERM2GENE_subset,
                                                     TERM2NAME = NA)
        
        
        
        #eps = 1e-10,
        
        cat("enricher_object\n")
        cat(str(enricher_object))
        cat("\n")
        
        # quit(status = 1)
        
        FLAG_NULL<-sum(is.null(enricher_object))
        
        cat("FLAG_NULL\n")
        cat(str(FLAG_NULL))
        cat("\n")
        
        # quit(status = 1)
        
        if(FLAG_NULL == 0)
        {
          setwd(path6)
          
          saveRDS(enricher_object,file=paste("enricher_object_",Threshold_logpval,".rds",sep=''))
          
        }#FLAG_NULL == 0
        
      }#length(gene) >0
      
    
    }#file.exists("DE_RESULTS_FOR_GSEA.rds")
  }#length(SELECTED_VARS_UPDATED) >0
  
  #### from https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/ ----
  
  
 
  
  
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
  
 
  enricher(opt)
    
}

# echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
# echo "--PCHiC_info $PCHiC_info \\" >> $output
# echo "--ALL_dB $ALL_dB \\" >> $output


###########################################################################

system.time( main() )
