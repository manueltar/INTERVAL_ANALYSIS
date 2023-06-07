#!/bin/bash>
 
  
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/309_1.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output

SELECTED_VARS=$(echo 'chr10_9032555_C_G,chr1_101744646_C_G,chr1_158613314_G_A,chr1_158638260_G_A,chr4_1008212_C_T,chr14_74218569_C_T,chr12_110331967_G_A,chr2_169313518_C_T,chr12_112092462_G_A,chr12_22610292_G_C,chr12_53309769_A_G,chr9_135857646_C_T,chr7_50444152_G_T,chr14_23587046_A_C,chr14_51132622_T_C,chr6_7143859_C_A,chr14_74630170_A_G,chr1_202129205_G_A,chr8_130641322_C_T,chr15_65174494_A_G,chr15_65174438_C_A,chr16_24761046_G_A,chr16_333719_G_C,chr16_67250992_C_T,chr16_67690688_C_T,chr16_85595360_G_C,chr16_86016328_C_T,chr16_89094897_C_T,chr17_16949211_C_A,chr17_27197056_G_T,chr17_27778073_C_T,chr17_38764524_T_A,chr17_47856909_A_G,chr17_56339594_A_C,chr13_28604007_T_C,chr17_58602131_G_A,chr17_7106378_G_A,chr18_42041131_T_G,chr18_60880701_T_C,chr18_60920854_C_T,chr18_67856078_G_A,chr19_10676941_G_A,chr8_41589736_T_G,chr6_34947254_A_G,chr1_91606142_A_G,chr9_114663385_T_C,chr19_11210157_C_T,chr19_3203962_C_T,chr1_93482787_G_A,chr19_35776481_C_T,chr20_25409287_A_G,chr20_37544151_C_T,chr20_55990370_A_T,chr2_144084356_C_T,chr2_144162105_A_G,chr19_15653669_T_C,chr22_18252442_G_A,chr17_56603493_C_T,chr22_28761148_C_T,chr22_39362450_C_T,chr2_24091099_C_T,chr1_92981236_T_G,chr2_31476771_G_C,chr2_46293826_C_T,chr2_74920648_G_A,chr3_128317978_C_T,chr3_128322617_G_A,chr3_17098399_A_G,chr3_184091102_T_G,chr3_46354444_C_T,chr3_71355240_G_C,chr6_41924998_C_T,chr5_1041433_C_T,chr5_1093511_G_A,chr6_41952511_T_G,chr5_75563535_G_A,chr5_35476470_G_T,chr22_50949811_T_C,chr7_100309180_A_G,chr1_92925654_G_C,chr6_82476412_C_T,chr6_82501719_C_T,chr7_100083971_G_A,chr1_29217311_G_A,chr7_100314474_C_T,chr7_101499930_G_A,chr7_139875230_C_T,chr2_219020958_C_T,chr16_155132_C_T,chr8_130429059_G_A,chr7_99760955_G_C,chr1_198680015_G_A,chr8_90995426_C_T,chr12_111844956_C_T,chr15_64349614_G_A,chr9_135874752_G_A,chr9_135920196_C_T,chr9_20576825_C_A,chr9_4533390_C_G,chr9_5079248_T_C')






MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4



#rm -rf $MASTER_ROUTE


#### This part of the script finds the SNPS in the WGS vcf files



Rscript_Mapping_VAR_PVAR_to_VCF=/nfs/users/nfs_m/mt19/Scripts/R/293_NEW_FOLLOW_UP_PROXYS_vALL.R

Proxys_VEP_CSQ=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/INTERVAL/Proxys_RRVV/Proxys/VEP_parsed_Proxys.txt')
Proxys_Rosetta=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/INTERVAL/Proxys_RRVV/Preparation_LD_Proxy_VAR.tsv')


type=$(echo "Mapping_VAR_PVAR_to_VCF")

outfile_Mapping_VAR_PVAR_to_VCF=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_Mapping_VAR_PVAR_to_VCF
echo -n "" > $outfile_Mapping_VAR_PVAR_to_VCF
name_Mapping_VAR_PVAR_to_VCF=$(echo "$type""_job")



echo "bsub -G team151 -o $outfile_Mapping_VAR_PVAR_to_VCF -M $mem  -J $name_Mapping_VAR_PVAR_to_VCF  -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_Mapping_VAR_PVAR_to_VCF \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--Proxys_VEP_CSQ $Proxys_VEP_CSQ \\" >> $output
echo "--Proxys_Rosetta $Proxys_Rosetta \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



VCF_files=$(echo 'vcf_file_chr10.vcf,vcf_file_chr13.vcf,vcf_file_chr15.vcf,vcf_file_chr17.vcf,vcf_file_chr19.vcf,vcf_file_chr20.vcf,vcf_file_chr2.vcf,vcf_file_chr4.vcf,vcf_file_chr6.vcf,vcf_file_chr8.vcf,vcf_file_chr12.vcf,vcf_file_chr14.vcf,vcf_file_chr16.vcf,vcf_file_chr18.vcf,vcf_file_chr1.vcf,vcf_file_chr22.vcf,vcf_file_chr3.vcf,vcf_file_chr5.vcf,vcf_file_chr7.vcf,vcf_file_chr9.vcf')

VCF_files_array=($(echo "$VCF_files" | tr "," '\n'))

EGAN_files_path=$(echo "$MASTER_ROUTE""EGAN_files""/")

# rm -rf $EGAN_files_path
# mkdir -p $EGAN_files_path

declare -a arr

for i  in "${VCF_files_array[@]}"
    do
        VCF_files_sel=${i}
        echo "$VCF_files_sel"

	CHROM=$(echo $VCF_files_sel|sed -e 's/^vcf_file_//g')
	CHROM=$(echo $CHROM|sed -e 's/\.vcf$//g')

	echo "$CHROM"

#	path_file=$(echo "/lustre/scratch126/humgen/projects/interval_wgs/release/""$CHROM""/""$CHROM"".intervalwgs_v2_GT_only.eagle_phased.vcf.gz")
	path_file=$(echo "/lustre/scratch126/humgen/projects/interval_wgs/final_release_freeze/gt_phased""/""interval_wgs.""$CHROM"".gt_phased.vcf.gz")


	
	echo "$path_file"
	
	end_of_file=0
	while [[ $end_of_file == 0 ]]
	do
	    read -r line
	    end_of_file=$?
	    FRIEND=$line
#	    echo "LINE: $FRIEND"

	    a=($(echo "$FRIEND" | tr "\t" '\n'))

	    if [[ "$line" =~ [^[:space:]] ]]; then
		# echo "chr: ${a[0]}"
		# echo "pos: ${a[1]}"


		chr=$(echo "chr""${a[0]}")
		pos=${a[1]}
		pos2=${a[1]}

		VAR=${a[2]}

		output_file=$(echo "$EGAN_files_path""$VAR""_EGAN_file"".tsv")

		touch $output_file
		echo -n "" > $output_file

		interval=$(echo "$chr"":""$pos""-""$pos2")

#		echo "$interval"

		tabix -f -h $path_file $interval  >> $output_file

#		exit

	    fi

	done < "$VCF_files_sel"
done

bash $output
