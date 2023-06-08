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
      
          
          
                                                              
                                                              
  - Multiple testing correction for DE and DTU
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
