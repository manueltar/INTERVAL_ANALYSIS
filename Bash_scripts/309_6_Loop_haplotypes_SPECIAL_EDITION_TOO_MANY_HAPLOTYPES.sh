#!/bin/bash>
  
  
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/310_ME.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output

SELECTED_VARS=$(echo 'chr18_42041131_T_G,chr6_34947254_A_G,chr4_1008212_C_T,chr2_24091099_C_T,chr5_35476470_G_T,chr1_158613314_G_A')


MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4



#rm -rf $MASTER_ROUTE
 
####################################################################################################### LOOP HAPLOTYPES ###################################################################################################################################################################################
####################################################################################################### LOOP HAPLOTYPES ####################################
###############################################################################################################################################
####################################################################################################### LOOP HAPLOTYPES ###################################################################################################################################################################################
####################################################################################################### LOOP HAPLOTYPES ###################################################################################################################################################################################
####################################################################################################### LOOP HAPLOTYPES ###################################################################################################################################################################################
####################################################################################################### LOOP HAPLOTYPES ###################################################################################################################################################################################

queue=$(echo "long")

a=($(echo "$SELECTED_VARS" | tr "," '\n'))
 
declare -a arr

for i  in "${a[@]}"
    do
        VAR_sel=${i}
        echo "$VAR_sel"

        VAR_ROUTE=$(echo "$MASTER_ROUTE""/""$VAR_sel""/")


	Haplotypes_array_1=$(echo 'chr18_42041131_T_G__chr18_41921592_C_A,chr18_42041131_T_G__chr18_42025923_C_G,chr18_42041131_T_G__chr18_42205085_C_T,chr18_42041131_T_G__chr18_41932704_A_G,chr18_42041131_T_G__chr18_42088811_A_T,chr18_42041131_T_G__chr18_42209961_T_G,chr18_42041131_T_G__chr18_41953230_C_T,chr18_42041131_T_G__chr18_42104777_C_T')

	Haplotypes_array_2=$(echo 'chr6_34947254_A_G__chr6_34583858_T_C,chr6_34947254_A_G__chr6_34683940_C_T,chr6_34947254_A_G__chr6_34798985_A_T,chr6_34947254_A_G__chr6_35359978_G_A,chr6_34947254_A_G__chr6_34640870_C_T,chr6_34947254_A_G__chr6_34690380_G_A,chr6_34947254_A_G__chr6_34816070_A_T,chr6_34947254_A_G__chr6_35378040_A_G,chr6_34947254_A_G__chr6_34654245_C_T,chr6_34947254_A_G__chr6_34750671_C_T,chr6_34947254_A_G__chr6_34919656_C_T,chr6_34947254_A_G__chr6_35407202_G_A,chr6_34947254_A_G__chr6_34676324_G_A,chr6_34947254_A_G__chr6_34763012_G_A,chr6_34947254_A_G__chr6_35182249_C_T,chr6_34947254_A_G__chr6_35413408_C_G,chr6_34947254_A_G__chr6_34676324_G_T,chr6_34947254_A_G__chr6_34798985_A_C,chr6_34947254_A_G__chr6_35284939_A_G,chr6_34947254_A_G__chr6_35413408_C_T')

	Haplotypes_array_3=$(echo 'chr4_1008212_C_T__chr4_1000383_G_A,chr4_1008212_C_T__chr4_1007322_C_G,chr4_1008212_C_T__chr4_1014325_C_T,chr4_1008212_C_T__chr4_1018499_G_A,chr4_1008212_C_T__chr4_698456_C_A,chr4_1008212_C_T__chr4_802027_CATAA_C,chr4_1008212_C_T__chr4_802144_C_T,chr4_1008212_C_T__chr4_995868_C_T,chr4_1008212_C_T__chr4_1000383_G_T,chr4_1008212_C_T__chr4_1009703_G_T,chr4_1008212_C_T__chr4_1015204_C_T,chr4_1008212_C_T__chr4_1019011_C_T,chr4_1008212_C_T__chr4_698456_C_G,chr4_1008212_C_T__chr4_802027_C_CATAA,chr4_1008212_C_T__chr4_917017_C_T,chr4_1008212_C_T__chr4_998777_C_T,chr4_1008212_C_T__chr4_1000799_C_T,chr4_1008212_C_T__chr4_1010030_T_C,chr4_1008212_C_T__chr4_1015789_T_C,chr4_1008212_C_T__chr4_1099797_T_A,chr4_1008212_C_T__chr4_749620_T_G,chr4_1008212_C_T__chr4_802027_C_CATAAATAA,chr4_1008212_C_T__chr4_993388_G_A')

	Haplotypes_array_4=$(echo 'chr2_24091099_C_T__chr2_24109099_C_A,chr2_24091099_C_T__chr2_24109099_C_T,chr2_24091099_C_T__chr2_24134679_T_C,chr2_24091099_C_T__chr2_24245713_C_G,chr2_24091099_C_T__chr2_24613564_C_T,chr2_24091099_C_T__chr2_24664837_G_A')

	Haplotypes_array_5=$(echo 'chr5_35476470_G_T__chr5_35303868_G_T,chr5_35476470_G_T__chr5_35396231_C_T,chr5_35476470_G_T__chr5_35591234_C_T,chr5_35476470_G_T__chr5_35742361_A_G,chr5_35476470_G_T__chr5_35748576_G_A,chr5_35476470_G_T__chr5_35788555_G_A')

	Haplotypes_array_6=$(echo 'chr1_158613314_G_A__chr1_158384566_G_T,chr1_158613314_G_A__chr1_158427208_A_C,chr1_158613314_G_A__chr1_158503018_T_C,chr1_158613314_G_A__chr1_158521490_A_G,chr1_158613314_G_A__chr1_158803254_G_A')

      if [ $VAR_sel == "chr18_42041131_T_G" ]; then

	  haplotypes_DEF=$Haplotypes_array_1
      else
	  if [ $VAR_sel == "chr6_34947254_A_G" ]; then
	      
              haplotypes_DEF=$Haplotypes_array_2
	  else

	      if [ $VAR_sel == "chr4_1008212_C_T" ]; then
		  
		  haplotypes_DEF=$Haplotypes_array_3
	      else
		  if [ $VAR_sel == "chr2_24091099_C_T" ]; then
		      
		      haplotypes_DEF=$Haplotypes_array_4
		  else
		       if [ $VAR_sel == "chr5_35476470_G_T" ]; then

			   haplotypes_DEF=$Haplotypes_array_5
                       else


			   haplotypes_DEF=$Haplotypes_array_6
		       fi
		      
		  fi  
		  
	      fi  
	  fi

      fi

      
      b=($(echo "$haplotypes_DEF" | tr "," '\n'))

      for j  in "${b[@]}"
      do
        Haplotype_sel=${j}

	Haplotypes=$Haplotype_sel
      ### DE LM Haplotypes ###################

      Rscript_DE_LM_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/300_OVERHAUL_DE_partLM_Haplotypes_SPECIAL_VERSION_INDIVIDUAL_HAPLOTYPES.R


      type=$(echo "DE_LM_Haplotypes""_""$VAR_sel""_""$Haplotypes")
      outfile_DE_LM_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_DE_LM_Haplotypes
      echo -n "" > $outfile_DE_LM_Haplotypes
      name_DE_LM_Haplotypes=$(echo "$type""_job")

      GENES_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_GENES_table.txt")
      INTERVAL_GENE_EXP=$(echo "/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/processed_v97/manuel/covariates_peer/ge_matrix_residuals_scaled_final.csv")


      DE_mem=$(expr $mem \* 2)
      DE_pc=$(expr $pc \* 2)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"



      echo "bsub -G team151 -o $outfile_DE_LM_Haplotypes -M $DE_mem -J $name_DE_LM_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_DE_LM_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_table $GENES_table \\" >> $output
      echo "--Haplotypes $Haplotypes \\" >> $output
      echo "--INTERVAL_GENE_EXP $INTERVAL_GENE_EXP \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output
      
###############       MT_correction_DE_Haplotypes

      Rscript_MT_correction_DE_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/303_OVERHAUL_DE_partMultiple_Testing_Main_Haplotype_SPECIAL_VERSION_INDIVIDUAL_HAPLOTYPES.R

      type=$(echo "MT_correction_DE_Haplotypes""_""$VAR_sel""_""$Haplotypes")
      outfile_MT_correction_DE_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_MT_correction_DE_Haplotypes
      echo -n "" > $outfile_MT_correction_DE_Haplotypes
      name_MT_correction_DE_Haplotypes=$(echo "$type""_job")


      ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
      GENES_PER_BLOCKS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/Allelic_Series_candidates_Allelic_Series_Generation.tsv")
      PCHiC_info=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
      VEP_CSQ=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/VEP_CSQ.csv")

      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


      echo "bsub -G team151 -o $outfile_MT_correction_DE_Haplotypes -M $DE_mem -w\"done($name_DE_LM_Haplotypes)\" -J $name_MT_correction_DE_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#      echo "bsub -G team151 -o $outfile_MT_correction_DE_Haplotypes -M $DE_mem -J $name_MT_correction_DE_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_MT_correction_DE_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
      echo "--PCHiC_info $PCHiC_info \\" >> $output
      echo "--Haplotypes $Haplotypes \\" >> $output
      echo "--VEP_CSQ $VEP_CSQ \\" >> $output
      echo "--ALL_dB $ALL_dB \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output


	

 ##############     DTU LM Haplotypes ###################

      Rscript_DTU_logRatio_LM_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/301_OVERHAUL_DTU_part_LogRatio_LM_Haplotypes_SPECIAL_VERSION_INDIVIDUAL_HAPLOTYPES.R


      type=$(echo "DTU_logRatio_LM_Haplotypes""_""$VAR_sel""_""$Haplotypes")
      outfile_DTU_logRatio_LM_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_DTU_logRatio_LM_Haplotypes
      echo -n "" > $outfile_DTU_logRatio_LM_Haplotypes
      name_DTU_logRatio_LM_Haplotypes=$(echo "$type""_job")


      TRANSCRIPTS_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_Transcripts_table.txt")
      INTERVAL_isoform_EXP=$(echo "/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/salmon/tpm_matrix/post_gene_qc/txi_transcript_log2_tmm_tpm_smpls_rmvd_swapped_post_geneqc_01tpm_20perc.csv")
      Transposed_Isoform_Expression=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/Transposed_Isoform_Expression_df.rds")



      DTU_mem=$(expr $mem \* 2)
      DTU_pc=$(expr $pc \* 2)


      echo "$mem""->""$DTU_mem"
      echo "$pc""->""$DTU_pc"


 
      echo "bsub -G hematopoiesis -o $outfile_DTU_logRatio_LM_Haplotypes -M $DTU_mem -J $name_DTU_logRatio_LM_Haplotypes -R\"select[model==Intel_Platinum]\" -R\" select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_DTU_logRatio_LM_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--Haplotypes $Haplotypes \\" >> $output
      echo "--TRANSCRIPTS_table $TRANSCRIPTS_table \\" >> $output
      echo "--INTERVAL_isoform_EXP $INTERVAL_isoform_EXP \\" >> $output
      echo "--Transposed_Isoform_Expression $Transposed_Isoform_Expression \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

     #   ### MT_correction_DTU_Haplotypes
 
      Rscript_MT_correction_DTU_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/305_OVERHAUL_DTU_partMultiple_Testing_Haplotype_SPECIAL_VERSION_INDIVIDUAL_HAPLOTYPES.R

      type=$(echo "MT_correction_DTU_Haplotypes""_""$VAR_sel""_""$Haplotypes")
      outfile_MT_correction_DTU_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_MT_correction_DTU_Haplotypes
      echo -n "" > $outfile_MT_correction_DTU_Haplotypes
        name_MT_correction_DTU_Haplotypes=$(echo "$type""_job")


      ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
      GENES_PER_BLOCKS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/Allelic_Series_candidates_Allelic_Series_Generation.tsv")
      PCHiC_info=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
      VEP_CSQ=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/VEP_CSQ.csv")

      DTU_mem=$(expr $mem \* 1)
      DTU_pc=$(expr $pc \* 1)


      echo "$mem""->""$DTU_mem"
      echo "$pc""->""$DTU_pc"



      echo "bsub -G team151 -o $outfile_MT_correction_DTU_Haplotypes -M $DTU_mem -w\"done($name_MT_correction_DE_Haplotypes)\" -J $name_MT_correction_DTU_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
#      echo "bsub -G team151 -o $outfile_MT_correction_DTU_Haplotypes -M $DTU_mem -w\"done($name_DTU_logRatio_LM_Haplotypes) && done($name_MT_correction_DE_Haplotypes)\" -J $name_MT_correction_DTU_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_MT_correction_DTU_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
      echo "--Haplotypes $Haplotypes \\" >> $output
      echo "--PCHiC_info $PCHiC_info \\" >> $output
      echo "--VEP_CSQ $VEP_CSQ \\" >> $output
      echo "--ALL_dB $ALL_dB \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

#      bash $output
#      exit

      done
done
 
bash $output




