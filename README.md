# bash 309_1_FIND_SNPS_in_WGS.sh \<output_dir\> 4000 1 normal

- Finds the samples in interval WGS (EGAN IDs) that caryy the SNP in HET and HOM

# bash 309_2_Add_specific_covariates.sh \<output_dir\> 4000 1 normal

- Adds to each sample-genotype the set of covariates to be included:

  CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                        "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                        "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                        "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                        "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")

  specific covariates: see Dependencies/covariate_phenotype_correspondence.tsv

# nohup bash 309_6_Loop_haplotypes_SPECIAL_EDITION_TOO_MANY_HAPLOTYPES.sh \<output_dir\> 4000 1 normal  > my_SPECIAL_Haplotypes_only.out 2>&1 &

  -For a subet of variants (chr18_42041131_T_G,chr6_34947254_A_G,chr4_1008212_C_T,chr2_24091099_C_T,chr5_35476470_G_T,chr1_158613314_G_A)
  -Haplotype per haplotype DE model for all the haplotypes of a variant
  -Haplotype per haplotype DTU model for all the haplotypes of a variant
  -Haplotype per haplotype Multiple testing correction for DE and DTU
  
# bash 309_3_Loop_haplotypes.sh \<output_dir\> 4000 1 long

  - Run 309_6 before and this one after to get all the haplotypes in the put together part
  - Done for most of the variants except the ones that have many haplotypes (long queue cannot cope with them)
  - DE model for all the haplotypes of a variant
  - DTU model for all the haplotypes of a variant
  - Multiple testing correction for DE and DTU
  - Put together results DE + ALL by ALL correction
  - Put together results DTU
  - DE graphical representation
  - DTU graphical representation

# bash 309_4_Loop_main_VARS.sh \<output_dir\> 4000 1 normal

  - For all the variants and their proxys individually (no haplotypes)
  - DE model for all the haplotypes of a variant
  - DTU model for all the haplotypes of a variant
       1. Filter genes and transcripts:
         - if(n_transcripts_per_gene > 1) # genes to be considered for DTU have to have more than one transcript
         -  Summary_table<-as.data.frame(Transposed_Isoform_Expression_sel.m_ADAPTED_HET.dt[, .(median=round(as.numeric(summary(TPM)[3]),3)),
                                                                              by=key(Transposed_Isoform_Expression_sel.m_ADAPTED_HET.dt)],stringsAsFactors=F)   # Calculate the median TPM expression per genotype (Hom_ref vs HET)
         -  Summary_table_TOTAL<-as.data.frame(Summary_table.dt[, .(TOTAL_GENE_EXP_median=sum(median)),
                                                              by=key(Summary_table.dt)],stringsAsFactors=F) # Calculate total gene expression as the summatory of the medians per genotype
         
         -  Summary_table$Transcript_Ratio<-(Summary_table$median/Summary_table$TOTAL_GENE_EXP_median) # calculate Transcript Ratio as the median transcript ratio per genotype divided by the total gene expression by gentoype
         - Summary_table.dt<-data.table(Summary_table,
                                       key=c("transcript_id")) \ change the key to transcript
         - Summary_table_FILTERED<-as.data.frame(Summary_table.dt[,.SD[which.max(Transcript_Ratio)],
                                                                 by=key(Summary_table.dt)],stringsAsFactors=F) # get the maximum transcript ratio by transcript (it will happen in any of the two genotypes Hom_ref or HET)
         -  Summary_table_FILTERED<-Summary_table_FILTERED[which(Summary_table_FILTERED$Transcript_Ratio >= 0.1),] # filter transcripts that make 10% or more of the Transcript ratio in at least one of the two genotypes
         -  saveRDS(file=paste("DTU_PASS_Transcripts_",SELECTED_VARS_UPDATED_sel,'.rds', sep=''), TRANSCRIPTS_table_FILTERED_FINAL) # save the filtered genes and transcripts

      2. DTU model      
          - if(n_transcripts_per_gene > 1) # only for genes that have two or more transcripts after the transcript filtering
          - Impute the 0 to 0.65 of the lowest expression:
              - Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO<-Transposed_Isoform_Expression_sel.m_REDUCED[which(Transposed_Isoform_Expression_sel.m_REDUCED$TPM >0),]
              - Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO,
                                                                     key=c("transcript_id"))
              - Zero_imputation_values<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO.dt[,.(min_TPM=min(TPM)),by=key(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO.dt)], stringsAsFactors=F) # min expression value for transcript irrespective of the genotype
              - Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$TPM<-0.65*Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$min_TPM # All the values with 0 expression are imputed to 0.65 of the minimum
              - Transposed_Isoform_Expression_sel.m_REDUCED<-rbind(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO,
                                                     Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO) # Rejoin the dataframes
                                                     
                                                     
           -  Obtain the reference transcript for the model, the transcript with the highest mean expression across the genoytpes
              -  Transposed_Isoform_Expression_sel.m_REDUCED.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED,
                                                             key=c("transcript_id"))
              -  Reference_value<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED.dt[,.(mean_TPM=mean(TPM)),by=key(Transposed_Isoform_Expression_sel.m_REDUCED.dt)], stringsAsFactors=F)
              -   Reference_value_MAX<-Reference_value[which(Reference_value$mean_TPM == max_value),]
              -   Reference_transcript<-Reference_value_MAX$transcript_id[1] # here we get the transcript_id of the transcript with the highest mean expression across genotypes for a given gene. This is done without having corrected the expression by the covariates

            - Calculate transcript ratios
              - Transposed_Isoform_Expression_sel.m_REDUCED.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED, key=c("sample_id"))
              - Summary_table_GENE_EXP<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED.dt[, .(sum_GENE_EXP=sum(TPM)), by=key(Transposed_Isoform_Expression_sel.m_REDUCED.dt)],stringsAsFactors=F)
              - Transposed_Isoform_Expression_sel.m_REDUCED<-merge(Transposed_Isoform_Expression_sel.m_REDUCED, Summary_table_GENE_EXP,by="sample_id")
              - Transposed_Isoform_Expression_sel.m_REDUCED.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED,key=c("sample_id","transcript_id"))
              - Ratio_df<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED.dt[, .(TPM=TPM,sum_GENE_EXP=sum_GENE_EXP, Ratio=TPM/sum_GENE_EXP),by=key(Transposed_Isoform_Expression_sel.m_REDUCED.dt)],stringsAsFactors=F)

            - Scale the transcript ratios to the value of the reference transcript for the model
               - Ratio_df_reference<-droplevels(Ratio_df[which(Ratio_df$transcript_id == Reference_transcript),c(which(colnames(Ratio_df) == "sample_id"), which(colnames(Ratio_df) == "transcript_id"), which(colnames(Ratio_df) == "Ratio"))])
               - colnames(Ratio_df_reference)[which(colnames(Ratio_df_reference) == "Ratio")]<-"Reference_ratio_value"
               - Ratio_df<-merge(Ratio_df,Ratio_df_reference,by="sample_id",all.x=T)
               - Ratio_df$scaled_ratio<-Ratio_df$Ratio/Ratio_df$Reference_ratio_value # Here we scaled all the ratio values to the ratio value of the reference transcript
             - Calculate the log
               - Ratio_df$logscaled_ratio<-log(Ratio_df$scaled_ratio)
             - Full Linear regression model per transcript
                - for(k in 1:length(ENST_array))
                - Ratio_df_ENST_sel<-merge(Ratio_df_ENST_sel,INTERVAL_covariates_and_PEER_factors_sel,by="sample_id")
                - ACCEPTED_genotypes<-c("HOM_REF","HET")
                - Ratio_df_ENST_sel_HET<-Ratio_df_ENST_sel[which(Ratio_df_ENST_sel$Genotype%in%ACCEPTED_genotypes),]
                - Ratio_df_ENST_sel_HET<-droplevels(Ratio_df_ENST_sel_HET)
                - unselected_columns<-c("sample_id","transcript_id","TPM","sum_GENE_EXP","Ratio","Reference_ratio_value","scaled_ratio")
                - model<-lm(logscaled_ratio ~ ., data=Ratio_df_ENST_sel_HET_LM) # This is the full model, logscaled_ratio ratio against all the covariates (constitutive and specific) + Genotype
                - results<-A$coefficients
                - pvalue_Genotypes_1<-as.numeric(results[which(row.names(results) == "Genotype.L"),4])
                - coefficient_Genotypes_1<-as.numeric(results[which(row.names(results) == "Genotype.L"),1])
                - A.df<-as.data.frame(cbind(ENST_array_sel,pvalue_Genotypes_1, coefficient_Genotypes_1, RUN_OUT_OF_NAMES_DEF))
                - colnames(A.df)<-c("transcript_id","pvalue_Genotypes_specific_CELL_COUNTS","coefficient_Genotypes_specific_CELL_COUNTS","n_breakdown_string")
                - A.df$pvalue_Genotypes_specific_CELL_COUNTS<-as.numeric(A.df$pvalue_Genotypes_specific_CELL_COUNTS)
                - A.df$coefficient_Genotypes_specific_CELL_COUNTS<-as.numeric(A.df$coefficient_Genotypes_specific_CELL_COUNTS)
           - Reduced linear model to plot residuals + intercept
                - unselected_columns<-c("sample_id","transcript_id","sum_GENE_EXP","TPM","Reference_ratio_value","scaled_ratio","Genotype","logscaled_ratio") # exclude the genotype from the analysis
                - model<-lm(Ratio ~ ., data=Ratio_df_ENST_sel_HET_LM) # reduced model with the ratio not the logscaled ratio
                -  residual_results=residuals(model)
                -  residual_names<-names(residual_results)
                -  colnames(residual_results)<-"residuals_full_model" # Obtain residuals
                -  intercept<-summary_model.m$value[which(summary_model.m$Parameters == "Estimate" & summary_model.m$Terms == "(Intercept)")] # Obtain intercept
                -   residual_df<-as.data.frame(cbind(residual_names,residual_results),stringsAsFactors=F)
                -   colnames(residual_df)<-c("sample_id","residuals")
                -   residual_df$residuals<-as.numeric(residual_df$residuals)
                -   residual_df$residuals<-as.numeric(intercept) + residual_df$residuals # Add intercept to residuals                                                            
  - Multiple testing correction for DE
  - Multiple testing correction for DTU:
             - CIS consequences
                - Per variant
                - CIS_LABELS <-c("LOF","MISS","SYN","UTR5","UTR3","INTRON","SPLICE","UPSTREAM")
                - VEP_CSQ_sel_CIS<-VEP_CSQ_sel[which(VEP_CSQ_sel$VEP_DEF_LABELS%in%CIS_LABELS),] # Only variants that have CIS consequences in the genes
                - Results_Nominal_sel<-Results_Nominal[which(Results_Nominal$ensembl_gene_id%in%ENSG_array),] # Select nominal pvalues for the genes in which the variant has CIS consequences
                - MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Genotypes_specific_CELL_COUNTS, method = "BH")
                - MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS),2) # Adjust the pvalues to the tests (different transcripts of the CIS genes)
             - in_block_Plus_PCHi_C
                -  AS_Gene_source<-merge(ALL_FINAL_sel,GENES_PER_BLOCKS,by="Allelic_Series_ID") # select all the genes in the GWAS block
                -  PCHiC_info_sel<-PCHiC_info[which(PCHiC_info$VAR%in%SELECTED_VARS_UPDATED),] # Select all the genes linked to the variant by PCHiC
                -  ENSG_sel<-unique(AS_Gene_source_sel$ensembl_gene_id)
                -  ENSG_array<-c(ENSG_array,ENSG_sel)
                -  ENSG_sel<-unique(PCHiC_info_sel$ensembl_gene_id)
                -  ENSG_array<-c(ENSG_array,ENSG_sel) # Add both subsets of genes to the space of genes to be corrected
                -  Results_Nominal_sel<-Results_Nominal[which(Results_Nominal$ensembl_gene_id%in%ENSG_array),] # select nominal pvalues
                -  MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Genotypes_specific_CELL_COUNTS, method = "BH")
                -  MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS),2) # Correct the pvalues
              -  correction_genome_wide
                -  MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS<-p.adjust(MT_set$pvalue_Genotypes_specific_CELL_COUNTS, method = "BH")
                -  MT_set$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes_specific_CELL_COUNTS),2) # Select all the genes with DTU model for a variant and correct (all genesx all transcripts)
             
  - Put together results DE + ALL by ALL correction
  - Put together results DTU
  - DE graphical representation
  - DTU graphical representation

# bash 309_5_GSEA_analysis.sh \<output_dir\> 4000 1 normal

  - Order all genes (no DE) by Fold Chnage between genotypes
  - Perform GSEA analysis using GO terms. Specifications:

    	minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05

  - Perform ORA analysis
  - Put together results
  - Graphical representation

# bash 304_Figure4_pannels_v3.sh \<output_dir\> 4000 1 normal

  - Manual curation results into the table
  - Enrichments for Manual curation vs MPRA

# bash 314_explore_science_adh_7699.sh \<output_dir\> 4000 1 normal

  - Compare out set of variants with the results from:

  Morris, John A., Christina Caragine, Zharko Daniloski, Júlia Domingo, Timothy Barry, Lu Lu, Kyrie Davis, et al. 2023. “Discovery of Target Genes and Pathways at GWAS Loci by Pooled Single-Cell CRISPR Screens.” Science, May, eadh7699.
