#!/bin/bash>
 
  
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript
 
output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/309_5.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output

SELECTED_VARS=$(echo 'chr10_9032555_C_G,chr1_101744646_C_G,chr1_158613314_G_A,chr1_158638260_G_A,chr4_1008212_C_T,chr14_74218569_C_T,chr12_110331967_G_A,chr2_169313518_C_T,chr12_112092462_G_A,chr12_22610292_G_C,chr12_53309769_A_G,chr9_135857646_C_T,chr7_50444152_G_T,chr14_23587046_A_C,chr14_51132622_T_C,chr6_7143859_C_A,chr14_74630170_A_G,chr1_202129205_G_A,chr8_130641322_C_T,chr15_65174494_A_G,chr15_65174438_C_A,chr16_24761046_G_A,chr16_333719_G_C,chr16_67250992_C_T,chr16_67690688_C_T,chr16_85595360_G_C,chr16_86016328_C_T,chr16_89094897_C_T,chr17_16949211_C_A,chr17_27197056_G_T,chr17_27778073_C_T,chr17_38764524_T_A,chr17_47856909_A_G,chr17_56339594_A_C,chr13_28604007_T_C,chr17_58602131_G_A,chr17_7106378_G_A,chr18_42041131_T_G,chr18_60880701_T_C,chr18_60920854_C_T,chr18_67856078_G_A,chr19_10676941_G_A,chr8_41589736_T_G,chr6_34947254_A_G,chr1_91606142_A_G,chr9_114663385_T_C,chr19_11210157_C_T,chr19_3203962_C_T,chr1_93482787_G_A,chr19_35776481_C_T,chr20_25409287_A_G,chr20_37544151_C_T,chr20_55990370_A_T,chr2_144084356_C_T,chr2_144162105_A_G,chr19_15653669_T_C,chr22_18252442_G_A,chr17_56603493_C_T,chr22_28761148_C_T,chr22_39362450_C_T,chr2_24091099_C_T,chr1_92981236_T_G,chr2_31476771_G_C,chr2_46293826_C_T,chr2_74920648_G_A,chr3_128317978_C_T,chr3_128322617_G_A,chr3_17098399_A_G,chr3_184091102_T_G,chr3_46354444_C_T,chr3_71355240_G_C,chr6_41924998_C_T,chr5_1041433_C_T,chr5_1093511_G_A,chr6_41952511_T_G,chr5_75563535_G_A,chr5_35476470_G_T,chr22_50949811_T_C,chr7_100309180_A_G,chr1_92925654_G_C,chr6_82476412_C_T,chr6_82501719_C_T,chr7_100083971_G_A,chr1_29217311_G_A,chr7_100314474_C_T,chr7_101499930_G_A,chr7_139875230_C_T,chr2_219020958_C_T,chr16_155132_C_T,chr8_130429059_G_A,chr7_99760955_G_C,chr1_198680015_G_A,chr8_90995426_C_T,chr12_111844956_C_T,chr15_64349614_G_A,chr9_135874752_G_A,chr9_135920196_C_T,chr9_20576825_C_A,chr9_4533390_C_G,chr9_5079248_T_C')




MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4

 
#SELECTED_VARS=$(echo 'chr7_101499930_G_A')

a=($(echo "$SELECTED_VARS" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        VAR_sel=${i}
       echo "$VAR_sel"
	VAR_ROUTE=$(echo "$MASTER_ROUTE""/""$VAR_sel""/")


	Rscript_FC_calculation_for_GSEA=/nfs/users/nfs_m/mt19/Scripts/R/314_OVERHAUL_DE_FC_for_GSEA_Main_VARS.R

	type=$(echo "FC_calculation_for_GSEA_""$VAR_sel")
	outfile_FC_calculation_for_GSEA=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
	touch $outfile_FC_calculation_for_GSEA
	echo -n "" > $outfile_FC_calculation_for_GSEA
	name_FC_calculation_for_GSEA=$(echo "$type""_job")

	DE_mem=$(expr $mem \* 1)
	DE_pc=$(expr $pc \* 1)


	Results_INTERVAL_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/ALL_BY_ALL_DE_LM_FPKM_results_Main_VARS.tsv")

	echo "$mem""->""$DE_mem"
	echo "$pc""->""$DE_pc"



	#echo "bsub -G team151 -o $outfile_FC_calculation_for_GSEA -M $DE_mem -w\"done($name_PUT_TOGETHER_RESULTS_Main_VARS_ALL_BY_ALL)\" -J $name_FC_calculation_for_GSEA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "bsub -G team151 -o $outfile_FC_calculation_for_GSEA -M $DE_mem -J $name_FC_calculation_for_GSEA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_FC_calculation_for_GSEA \\" >> $output
	echo "--SELECTED_VARS $VAR_sel \\" >> $output
	echo "--Results_INTERVAL_LM $Results_INTERVAL_LM \\" >> $output
	echo "--type $type --out $MASTER_ROUTE\"" >> $output

	######################################### PREPARER_GSEA_AND_ORA
	
	Rscript_PREPARER_GSEA_AND_ORA=/nfs/users/nfs_m/mt19/Scripts/R/315_GSEA_AND_ORA_GW_PREPARER.R

	type=$(echo "PREPARER_GSEA_AND_ORA_""$VAR_sel")
	outfile_PREPARER_GSEA_AND_ORA=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
	touch $outfile_PREPARER_GSEA_AND_ORA
	echo -n "" > $outfile_PREPARER_GSEA_AND_ORA
	name_PREPARER_GSEA_AND_ORA=$(echo "$type""_job")

	DE_mem=$(expr $mem \* 1)
	DE_pc=$(expr $pc \* 1)


	Results_INTERVAL_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/ALL_BY_ALL_DE_LM_FPKM_results_Main_VARS.tsv")

	echo "$mem""->""$DE_mem"
	echo "$pc""->""$DE_pc"



	#echo "bsub -G team151 -o $outfile_PREPARER_GSEA_AND_ORA -M $DE_mem -w\"done($name_PUT_TOGETHER_RESULTS_Main_VARS_ALL_BY_ALL)\" -J $name_PREPARER_GSEA_AND_ORA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "bsub -G team151 -o $outfile_PREPARER_GSEA_AND_ORA -M $DE_mem -J $name_PREPARER_GSEA_AND_ORA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_PREPARER_GSEA_AND_ORA \\" >> $output
	echo "--SELECTED_VARS $VAR_sel \\" >> $output
	echo "--type $type --out $MASTER_ROUTE\"" >> $output

	######################################### GO_GSEA
	
	Rscript_GO_GSEA=/nfs/users/nfs_m/mt19/Scripts/R/316_GSEA.R

	type=$(echo "GO_GSEA_""$VAR_sel")
	outfile_GO_GSEA=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
	touch $outfile_GO_GSEA
	echo -n "" > $outfile_GO_GSEA
	name_GO_GSEA=$(echo "$type""_job")

	DE_mem=$(expr $mem \* 2)
	DE_pc=$(expr $pc \* 2)


	echo "$mem""->""$DE_mem"
	echo "$pc""->""$DE_pc"



	echo "bsub -G team151 -o $outfile_GO_GSEA -M $DE_mem -w\"done($name_PREPARER_GSEA_AND_ORA)\" -J $name_GO_GSEA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_GO_GSEA -M $DE_mem -J $name_GO_GSEA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_GO_GSEA \\" >> $output
	echo "--SELECTED_VARS $VAR_sel \\" >> $output
	echo "--type $type --out $MASTER_ROUTE\"" >> $output



      Threshold_logpval_string=$(echo '1.3_''5_''16_')

      echo $Threshold_logpval_string

      a=($(echo "$Threshold_logpval_string" | tr "_" '\n'))



      for j  in "${a[@]}"
      do


          Threshold_logpval=${j}
          echo "$Threshold_logpval"

	  	######################################### ORA
	
	Rscript_ORA=/nfs/users/nfs_m/mt19/Scripts/R/317_ORA.R

	type=$(echo "ORA_""$Threshold_logpval""_""$VAR_sel")
	outfile_ORA=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
	touch $outfile_ORA
	echo -n "" > $outfile_ORA
	name_ORA=$(echo "$type""_job")

	DE_mem=$(expr $mem \* 2)
	DE_pc=$(expr $pc \* 2)

	

	echo "$mem""->""$DE_mem"
	echo "$pc""->""$DE_pc"



	echo "bsub -G team151 -o $outfile_ORA -M $DE_mem -w\"done($name_PREPARER_GSEA_AND_ORA)\" -J $name_ORA -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_ORA \\" >> $output
	echo "--SELECTED_VARS $VAR_sel \\" >> $output
	echo "--Threshold_logpval $Threshold_logpval \\" >> $output
	echo "--type $type --out $MASTER_ROUTE\"" >> $output

	######################################### ORA_GO
	
	Rscript_ORA_GO=/nfs/users/nfs_m/mt19/Scripts/R/318_ORA_GO.R

	type=$(echo "ORA_GO_""$Threshold_logpval""_""$VAR_sel")
	outfile_ORA_GO=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
	touch $outfile_ORA_GO
	echo -n "" > $outfile_ORA_GO
	name_ORA_GO=$(echo "$type""_job")

	DE_mem=$(expr $mem \* 2)
	DE_pc=$(expr $pc \* 2)



	echo "$mem""->""$DE_mem"
	echo "$pc""->""$DE_pc"



	echo "bsub -G team151 -o $outfile_ORA_GO -M $DE_mem -w\"done($name_PREPARER_GSEA_AND_ORA)\" -J $name_ORA_GO -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_ORA_GO \\" >> $output
        echo "--SELECTED_VARS $VAR_sel \\" >> $output
        echo "--Threshold_logpval $Threshold_logpval \\" >> $output
	echo "--type $type --out $MASTER_ROUTE\"" >> $output

      done
	  
	

	#     # "chr1_101744646_C_G"
#     # "chr8_41589736_T_G"
#     # "chr3_46354444_C_T"
#       # "chr9_135920196_C_T"
	#       # "chr3_71355240_G_C"
	#"chr14_51132622_T_C"
     # if [ $VAR_sel == "chr1_101744646_C_G" ]; then

     #       bash $output
     #        exit
     #  else

     #        echo "Hello_world"" ""$VAR_sel"
     #        fi


done

# ############################### GSEA PUT TOGETHER ####################################################

Rscript_GSEA_PUT_TOGETHER_RESULTS=/nfs/users/nfs_m/mt19/Scripts/R/319_GSEA_GW_PUT_TOGETHER.R

type=$(echo "GSEA_PUT_TOGETHER_RESULTS_GSEA_GW")
outfile_GSEA_PUT_TOGETHER_RESULTS=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_GSEA_PUT_TOGETHER_RESULTS
echo -n "" > $outfile_GSEA_PUT_TOGETHER_RESULTS
name_GSEA_PUT_TOGETHER_RESULTS=$(echo "$type""_job")

DE_mem=$(expr $mem \* 4)
DE_pc=$(expr $pc \* 4)


echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"




echo "bsub -G team151 -o $outfile_GSEA_PUT_TOGETHER_RESULTS -M $DE_mem -w$done_string -J $name_GSEA_PUT_TOGETHER_RESULTS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_GSEA_PUT_TOGETHER_RESULTS -M $DE_mem -J $name_GSEA_PUT_TOGETHER_RESULTS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_GSEA_PUT_TOGETHER_RESULTS \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

# ############################### GSEA graphical_output ####################################################

Rscript_GSEA_GRAPHS=/nfs/users/nfs_m/mt19/Scripts/R/320_GSEA_graphical_output.R

type=$(echo "GSEA_GRAPHS_GSEA_GW")
outfile_GSEA_GRAPHS=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_GSEA_GRAPHS
echo -n "" > $outfile_GSEA_GRAPHS
name_GSEA_GRAPHS=$(echo "$type""_job")

DE_mem=$(expr $mem \* 2)
DE_pc=$(expr $pc \* 2)

GSEA_GLOBAL_RESULTS=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/GSEA_GW_RESULTS.tsv")
Results_INTERVAL_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/ALL_BY_ALL_DE_LM_FPKM_results_Main_VARS.tsv")
echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"




#echo "bsub -G team151 -o $outfile_GSEA_GRAPHS -M $DE_mem -w\"done($name_GSEA_PUT_TOGETHER_RESULTS)\" -J $name_GSEA_GRAPHS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "bsub -G team151 -o $outfile_GSEA_GRAPHS -M $DE_mem -J $name_GSEA_GRAPHS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_GSEA_GRAPHS \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--GSEA_GLOBAL_RESULTS $GSEA_GLOBAL_RESULTS \\" >> $output
echo "--Results_INTERVAL_LM $Results_INTERVAL_LM \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output
 


######################DTU graphical summary ----

### This one comes with function to plot residuals

Rscript_graphical_SUMMARY_DTU_Main_VARS=/nfs/users/nfs_m/mt19/Scripts/R/310_OVERHAUL_DTU_Graphical_summary_Main_VARS.R

type=$(echo "graphical_SUMMARY_DTU_Main_VARS")
outfile_graphical_SUMMARY_DTU_Main_VARS=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_graphical_SUMMARY_DTU_Main_VARS
echo -n "" > $outfile_graphical_SUMMARY_DTU_Main_VARS
name_graphical_SUMMARY_DTU_Main_VARS=$(echo "$type""_job")

DTU_mem=$(expr $mem \* 4)
DTU_pc=$(expr $pc \* 4)


Results_LogLM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/DTU_LM_logRatio_results_Main_VARS.tsv")
Transposed_Isoform_Expression=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/Transposed_Isoform_Expression_df.rds")
Tappas_gff=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/tappAS/Homo_sapiens_GRCh38_Ensembl_86.gff3")


echo "$mem""->""$DTU_mem"
echo "$pc""->""$DTU_pc"



echo "bsub -G team151 -o $outfile_graphical_SUMMARY_DTU_Main_VARS -M $DTU_mem -w\"done($name_PUT_TOGETHER_RESULTS_Main_VARS_DTU)\" -J $name_graphical_SUMMARY_DTU_Main_VARS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
# echo "bsub -G team151 -o $outfile_graphical_SUMMARY_DTU_Main_VARS -M $DTU_mem -J $name_graphical_SUMMARY_DTU_Main_VARS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_graphical_SUMMARY_DTU_Main_VARS \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--Results_LogLM $Results_LogLM \\" >> $output
echo "--Tappas_gff $Tappas_gff \\" >> $output
echo "--Transposed_Isoform_Expression $Transposed_Isoform_Expression \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

bash $output
