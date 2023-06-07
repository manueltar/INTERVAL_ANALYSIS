

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
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gtools", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
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
suppressMessages(library("liftOver",lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

# options(warn = 1)


graphical_function = function(option_list)
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
 
  
  #### metabolon_VARS ----
  
  
  
  metabolon_VARS  = unlist(strsplit(opt$metabolon_VARS , split=","))
  
  cat("metabolon_VARS_\n")
  cat(str(metabolon_VARS ))
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
  
  
  #### UPDATE metabolon_VARS ----
  
  indx.UPDATE<-which(metabolon_VARS %in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    metabolon_VARS = metabolon_VARS [-indx.UPDATE]
    
  }else{
    
    metabolon_VARS = metabolon_VARS 
  }
  
  
  
  cat("metabolon_VARS_\n")
  cat(str(metabolon_VARS))
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
  
  if(length(metabolon_VARS) >0)
  {
    #### Read Kousik matrix of metabolite abundance ----
    
    metabolon_1<-as.data.frame(fread(file=opt$metabolon_1, sep=",", header=T) , stringsAsFactors=F)
    
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_1\n")
      cat(str(metabolon_1))
      cat("\n")
    }
    
    metabolon_1.m<-melt(metabolon_1, id.vars=c("locus","alleles","s","GT","Hard_GT","VarID","rsid","Genotype"),value.name = "residuals_LM",)
    colnames(metabolon_1.m)[which(colnames(metabolon_1.m) == "s")]<-"EGAN_ID"
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_1.m\n")
      cat(str(metabolon_1.m))
      cat("\n")
    }
    
    metabolon_2<-as.data.frame(fread(file=opt$metabolon_2, sep=",", header=T) , stringsAsFactors=F)
    
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_2\n")
      cat(str(metabolon_2))
      cat("\n")
    }
    
    metabolon_2.m<-melt(metabolon_2, id.vars=c("locus","alleles","s","GT","Hard_GT","VarID","rsid","Genotype"),value.name = "residuals_LM")
    colnames(metabolon_2.m)[which(colnames(metabolon_2.m) == "s")]<-"EGAN_ID"
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_2.m\n")
      cat(str(metabolon_2.m))
      cat("\n")
    }
    
    metabolon_DEF<-rbind(metabolon_1.m,metabolon_2.m)
   
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_DEF_0\n")
      cat(str(metabolon_DEF))
      cat("\n")
    }
    
    metabolon_DEF$VAR<-gsub(":","_",metabolon_DEF$VarID)
    metabolon_DEF$variable<-gsub('phenotype.',"",metabolon_DEF$variable)
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_DEF_1\n")
      cat(str(metabolon_DEF))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(metabolon_DEF$variable))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(metabolon_DEF$variable)))))
      cat("\n")
    }
    
    metabolon_DEF$chr<-gsub("_.+$","",metabolon_DEF$VAR)
    metabolon_DEF$pos38<-gsub("^chr[^_]+_","",metabolon_DEF$VAR)
    metabolon_DEF$pos38<-gsub("_.+$","",metabolon_DEF$pos38)
    metabolon_DEF$ref<-gsub("^chr[^_]+_[^_]+_","",metabolon_DEF$VAR)
    metabolon_DEF$ref<-gsub("_.+$","",metabolon_DEF$ref)
    metabolon_DEF$alt<-gsub("^chr[^_]+_[^_]+_[^_]+_","",metabolon_DEF$VAR)
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_DEF_2\n")
      cat(str(metabolon_DEF))
      cat("\n")
    }
    
    # cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_Prob_subset$phenotype_DEF))))))
    # cat("\n")
    # cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_Prob_subset$phenotype_DEF)))))
    # cat("\n")
    
    
    
    ###### LiftOver 38 -> 37 ----
    
    
    gr_VARS <- GRanges(
      seqnames = as.character(gsub("chr","",metabolon_DEF$chr)),
      ranges=IRanges(
        start=as.numeric(metabolon_DEF$pos38),
        end=as.numeric(metabolon_DEF$pos38),
        name=metabolon_DEF$VAR))
    
    # cat("gr_VARS\n")
    # str(gr_VARS)
    # cat("\n")
    
    VAR_df<-unique(data.frame(chr=as.character(paste('chr',seqnames(gr_VARS), sep='')),
                       pos38=start(gr_VARS),
                       ref=metabolon_DEF$ref,
                       alt=metabolon_DEF$alt,
                       VAR=metabolon_DEF$VAR,
                       stringsAsFactors = F))
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_df_\n")
      str(VAR_df)
      cat("\n")
    }
    
    
    path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch = import.chain(path)
    
    seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
    gr_VARS37 = liftOver(gr_VARS, ch)
    gr_VARS37 = unlist(gr_VARS37)
    genome(gr_VARS37) = "hg19"
    
    if(length(gr_VARS37) >0)
    {
      
      chr_37<-as.character(seqnames(gr_VARS37))
      names_37<-as.character(names(gr_VARS37))
      
      ref_VAR38<-gsub("^chr[^_]+_[0-9]+_","",names_37)
      ref_VAR38<-gsub("_.+$","",ref_VAR38)
      
      
      # cat("ref_VAR38\n")
      # cat(sprintf(as.character(ref_VAR38)))
      # cat("\n")
      
      alt_VAR38<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_37)
      # alt_VAR38<-gsub("_.+$","",alt_VAR38)
      
      
      # cat("alt_VAR38\n")
      # cat(sprintf(as.character(alt_VAR38)))
      # cat("\n")
      
      
      
      
      VAR_37_df<-data.frame(chr=as.character(seqnames(gr_VARS37)),
                            pos_37=start(gr_VARS37),
                            ref=ref_VAR38,
                            alt=alt_VAR38,
                            VAR=names(gr_VARS37),
                            stringsAsFactors = F)
      
      VAR_37_df$VAR_37<-paste(VAR_37_df$chr,VAR_37_df$pos_37,VAR_37_df$ref,VAR_37_df$alt,sep='_')
      
      if(Condition_DEBUG == 1)
      {
        cat("VAR_37_df_1\n")
        str(VAR_37_df)
        cat("\n")
      }
      
      
      VAR_metabolon_DEF_df<-unique(merge(VAR_df,
                                                             VAR_37_df,
                                                             by=c("chr","ref","alt","VAR"),
                                                             all=T))
      
      VAR_metabolon_DEF_df$VAR_37[is.na(VAR_metabolon_DEF_df$VAR_37)]<-"ABSENT"
      
      if(Condition_DEBUG == 1)
      {
        cat("VAR_metabolon_DEF_df_2\n")
        str(VAR_metabolon_DEF_df)
        cat("\n")
      }
      
      check.ABSENT<-VAR_metabolon_DEF_df[which(VAR_metabolon_DEF_df$VAR_37 == "ABSENT"),]
      #
      if(Condition_DEBUG == 1)
      {
        cat("check.ABSENT\n")
        str(check.ABSENT) #7
        cat("\n")
      }
      
      
      
      metabolon_DEF<-unique(merge(metabolon_DEF,
                                  VAR_metabolon_DEF_df[which(VAR_metabolon_DEF_df$VAR_37 != "ABSENT"),],
                                  by=c("VAR","chr","pos38","ref","alt")))
      if(Condition_DEBUG == 1)
      {
        cat("df_2\n")
        str(df)
        cat("\n")
      }
      
    }else{
      
      
      stop("NO_LIFT_OVER\n")
      
    }# length(gr_VARS37) >0
    
    
    
    # metabolon_DEF$chr.pos_ref_alt<-gsub("^chr","",metabolon_DEF$VAR_38)
    # metabolon_DEF$chr.pos_ref_alt<-sub("_",":",metabolon_DEF$chr.pos_ref_alt)
    
    if(Condition_DEBUG == 1)
    {
      cat("metabolon_DEF_3\n")
      cat(str(metabolon_DEF))
      cat("\n")
      cat(str(unique(metabolon_DEF$VAR)))
      cat("\n")
    }
    
    
   
    
    
    ##### LOOP TO READ ALL VARIABLES -----
   
    
    for(i in 1:length(metabolon_VARS))
    {
      
      metabolon_VARS_sel<-metabolon_VARS[i]
      
      cat("------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(metabolon_VARS_sel)))
      cat("\n")
      
      path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
      
      # cat("path7\n")
      # cat(sprintf(as.character(path7)))
      # cat("\n")
      
      
      if (file.exists(path7)){
        
        
        
        
      } else {
        dir.create(file.path(path7))
        
      }
      
      
      
      
      path8<-paste(path7,metabolon_VARS_sel,'/', sep='')
      
      if (file.exists(path8)){
        
        
        
        
      } else {
        dir.create(file.path(path8))
        
      }
      
      metabolon_DEF_sel<-metabolon_DEF[which(metabolon_DEF$VAR_37 == metabolon_VARS_sel),]
      
      if(Condition_DEBUG == 1)
      {
        cat("metabolon_DEF_sel_1\n")
        cat(str(metabolon_DEF_sel))
        cat("\n")
      }
      
      metabolon_array<-unique(metabolon_DEF_sel$variable)
     
     
      ############################# VIOLIN PLOTS ---------------------------------------------------------------------------------------------------------------------------------------
      if(Condition_DEBUG == 1)
      {
        cat("metabolon_array_FOR_VIOLIN_PLOTS\n")
        cat(str(metabolon_array))
        cat("\n")
      }
      
     
      
      if(length(metabolon_array) >0)
      {
        Condition_DEBUG <- 1
        
        list_graphs<-list()
        for(z in 1:length(metabolon_array))
        {
          
          metabolon_array_sel<-metabolon_array[z]
          
          cat("---------------->\t")
          cat(sprintf(as.character(metabolon_array_sel)))
          cat("\n")
          
          metabolon_DEF_sel_metabolon_sel<-metabolon_DEF_sel[which(metabolon_DEF_sel$variable%in%metabolon_array_sel),]
          
          if(Condition_DEBUG == 1)
          {
            cat("metabolon_DEF_sel_metabolon_sel\n")
            cat(str(metabolon_DEF_sel_metabolon_sel))
            cat("\n")
            
          }
          
         
          
          A<-round(summary(metabolon_DEF_sel_metabolon_sel$residuals_LM[!is.na(metabolon_DEF_sel_metabolon_sel$residuals_LM)]),2)
          # A2<-round(summary(metabolon_DEF_sel_metabolon_sel$FPKM[!is.na(metabolon_DEF_sel_metabolon_sel$FPKM)]),2)
          
          metabolon_DEF_sel_metabolon_sel$residuals_LM_no_negative<-metabolon_DEF_sel_metabolon_sel$residuals_LM+abs(A[1])
          
          A3<-round(summary(metabolon_DEF_sel_metabolon_sel$residuals_LM_no_negative[!is.na(metabolon_DEF_sel_metabolon_sel$residuals_LM_no_negative)]),0)
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("metabolon_DEF_sel_metabolon_sel_2\n")
            cat(str(metabolon_DEF_sel_metabolon_sel))
            cat("\n")
            
            cat("Summary_residuals_LM:\t")
            cat(sprintf(as.character(names(A))))
            cat("\n")
            cat(sprintf(as.character(A)))
            cat("\n")
            
            # cat("Summary_FPKM:\t")
            # cat(sprintf(as.character(names(A2))))
            # cat("\n")
            # cat(sprintf(as.character(A2)))
            # cat("\n")
            
            
            cat("Summary_residuals_LM_no_negative:\t")
            cat(sprintf(as.character(names(A3))))
            cat("\n")
            cat(sprintf(as.character(A3)))
            cat("\n")
            
            
          }
          
          ########### Selects only HET genotypes ----------------------
          
          metabolon_DEF_sel_metabolon_sel$Genotype_2<-"NA"
          
          metabolon_DEF_sel_metabolon_sel$Genotype_2[which(metabolon_DEF_sel_metabolon_sel$GT%in%c("0/0"))]<-"HOM_REF"
          
          metabolon_DEF_sel_metabolon_sel$Genotype_2[which(metabolon_DEF_sel_metabolon_sel$GT%in%c("0/1","1/0"))]<-"HET"
          
          metabolon_DEF_sel_metabolon_sel$Genotype_2[which(metabolon_DEF_sel_metabolon_sel$GT%in%c("1/1"))]<-"HOM"
          
          metabolon_DEF_sel_metabolon_sel$Genotype_2<-factor(metabolon_DEF_sel_metabolon_sel$Genotype_2,
                                                             levels=c("HOM_REF","HET","HOM"),
                                                             ordered=T)
          if(Condition_DEBUG == 1)
          {
            cat("metabolon_DEF_sel_metabolon_sel_Genotype_2\n")
            cat(str(metabolon_DEF_sel_metabolon_sel))
            cat("\n")
            cat(sprintf(as.character(names(summary(metabolon_DEF_sel_metabolon_sel$Genotype_2)))))
            cat("\n")
            cat(sprintf(as.character(summary(metabolon_DEF_sel_metabolon_sel$Genotype_2))))
            cat("\n")
          }
          
          
          metabolon_DEF_sel_metabolon_sel_HET<-droplevels(metabolon_DEF_sel_metabolon_sel[which(metabolon_DEF_sel_metabolon_sel$Genotype_2%in%c("HOM_REF","HET")),])
          
          if(Condition_DEBUG == 1)
          {
            cat("metabolon_DEF_sel_metabolon_sel_HET_Genotype_0\n")
            cat(str(metabolon_DEF_sel_metabolon_sel_HET))
            cat("\n")
            cat(sprintf(as.character(names(summary(metabolon_DEF_sel_metabolon_sel_HET$Genotype_2)))))
            cat("\n")
            cat(sprintf(as.character(summary(metabolon_DEF_sel_metabolon_sel_HET$Genotype_2))))
            cat("\n")
          }
          
          
          ##### Violin visual FPKM_adjusted reduced model ----
          
       
          A<-round(summary(metabolon_DEF_sel_metabolon_sel_HET$residuals_LM_no_negative[!is.na(metabolon_DEF_sel_metabolon_sel_HET$residuals_LM_no_negative)]),2)
         
          step<-abs(A[6]-A[1])/10
          
          if(Condition_DEBUG == 1)
          {
            cat("summary_metabolon_DEF_sel_metabolon_sel_HET$residuals_LM_no_negative\n")
            cat(sprintf(as.character(names(A))))
            cat("\n")
            cat(sprintf(as.character(A)))
            cat("\n")
          }
          
          if(step == 0)
          {
            
            step<-1
          }
          
          breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
          labels.Rank<-as.character(round(breaks.Rank,3))
          
          levels_variable<-unique(metabolon_DEF_sel_metabolon_sel_HET$variable)
          
          if(Condition_DEBUG == 1)
          {
            cat("labels.Rank:\t")
            cat(sprintf(as.character(labels.Rank)))
            cat("\n")
            
            cat("levels_variable:\t")
            cat(sprintf(as.character(levels_variable)))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          
          
          graph_adjusted_metabolon<-ggplot(data=metabolon_DEF_sel_metabolon_sel_HET,aes(x=Genotype_2, y=residuals_LM_no_negative, fill=Genotype_2)) +
            geom_violin()+
            stat_summary(fun = median, fun.min = median, fun.max = median,
                         geom = "crossbar", width = 0.5)+
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_blank(),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_discrete(name=NULL, drop=F)+
            scale_y_continuous(name="Residuals model metabolites adjusted for covariates",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
            scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
            theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=12))+
            ggeasy::easy_center_title()
          
          # list_graphs[[z]]<-graph_adjusted_metabolon
        
          ##### Print per metabolite ----
          
          path9<-paste(path8,metabolon_array_sel,'/', sep='')
          
          if (file.exists(path9)){
            
            
            
            
          } else {
            dir.create(file.path(path9))
            
          }
          
          setwd(path9)
          # setwd(out)
          
          if(Condition_DEBUG == 1)
          {
            cat("path9:\t")
            cat(sprintf(as.character(path9)))
            cat("\n")
          }
          
          
          
          
          graph_DEF<-plot_grid(graph_adjusted_metabolon,
                               nrow = 1,
                               ncol = 1)
          
          
          
          
          
          
          
          svgname<-paste("METABOLITE_EVALUATION_",metabolon_array_sel,".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= graph_DEF,
                   device="svg",
                   height=10, width=12)
          }
          
          # ##############################################################
          # quit(status = 1)
        }#z in 1:length(metabolon_array
      }#length(metabolon_array) >0
    }# i in 1:length(metabolon_VARS)
  }#length(metabolon_VARS) >0
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
    make_option(c("--metabolon_VARS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--metabolon_1"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--metabolon_2"), type="character", default=NULL, 
                metavar="filename", 
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
  
 
  graphical_function(opt)
 
 
  
  
  
}


###########################################################################

system.time( main() )
