#!/bin/bash

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################



## extract specific samples from VCF file

# VCFExtractSamplesScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/DataFormatting/"
# toolUsed="extarctListSamplesFromVCF.sh"
# Compiler="bash"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData.vcf"
# keepList="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/sampleList.txt"
# outvcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_extracted.vcf"
# # activate the environment of vcftools_env from conda

# # conda activate vcftools_env
# $Compiler $VCFExtractSamplesScriptPath$toolUsed $vcffile $keepList $outvcffile



# #>>>>> filter VCF file
# VCF2GWASScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/SNP_QC/"
# toolUsed="vcf_filter.sh"
# Compiler="bash"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData.vcf"
# outvcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_filtered.vcf"
# # activate the environment of vcftools_env
# source activate vcftools_env
# $Compiler $VCF2GWASScriptPath$toolUsed $vcffile $outvcffile


# # >>>>>> Calculate PCA 
# GEMMAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PCA/"
# toolUsed="PCAbySNPrelate.r"
# Compiler="Rscript"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS/PCA"
# MetData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/MetaData.tsv"
# mkdir -p $outFolder

# $Compiler $GEMMAScriptPath$toolUsed $vcffile $outFolder $outFolder $MetData


# #>>>>>> Calculate Kinship for Field
# GEMMAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="KinshipBySNPrelate.r"
# Compiler="Rscript"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS/Kinship"
# mkdir -p $outFolder

# $Compiler $GEMMAScriptPath$toolUsed $vcffile $outFolder






# #>>>>>> Calculate Kinship for Field
# GEMMAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="KinshipBySNPrelate.r"
# Compiler="Rscript"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/Data/GenotypeData_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/RESULTS/Kinship"
# mkdir -p $outFolder

# $Compiler $GEMMAScriptPath$toolUsed $vcffile $outFolder




# ## Population structure using LEA
# LEAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="LEApopStructure.r"

# Compiler="Rscript"

# # for the field
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/Data/GenotypeData_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/RESULTS/LEA/Run"
# mkdir -p $outFolder

# $Compiler $LEAScriptPath$toolUsed $vcffile $outFolder 12 10000

# # plot the LEA results
# LEAplotScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="plotLEA.r"
# Compiler="Rscript"

# LEAoutFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/RESULTS/LEA/Run/bestKs"
# resultFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/RESULTS/LEA/plots"
# mkdir -p $resultFolder

# passportFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/Data/MetaData_new.tsv"
# orderFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_4_28_11_2023/RESULTS/LEA/Run/TaxaNames.txt"
# bestKnownK=6
# maxK=12

# $Compiler $LEAplotScriptPath$toolUsed $LEAoutFolder $resultFolder $passportFile $orderFile $bestKnownK $maxK

# # Check phenotypic data 
# # convert negative values to 0
# PhenotypicScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PhenotypicData/"
# toolUsed="convertNegativeValuesToZero.py"
# Compiler="python3"

# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/PhenoTypes/Field.tsv"
# $Compiler $PhenotypicScriptPath$toolUsed $phenData

# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/PhenoTypes/Seedling.tsv"
# $Compiler $PhenotypicScriptPath$toolUsed $phenData


# ## Match the vcf file with the phenotypic data
# VCFMatchPhenotypicScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/GWAS/"
# toolUsed="MatchVCFWithPhenotypicData.sh"
# Compiler="bash"

# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/PhenoTypes/Field.tsv"
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/BreadGenotype_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Field"
# mkdir -p $outFolder
# $Compiler $VCFMatchPhenotypicScriptPath$toolUsed $phenData $vcffile $outFolder

# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/PhenoTypes/Seedling.tsv"
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/BreadGenotype_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling"
# mkdir -p $outFolder
# $Compiler $VCFMatchPhenotypicScriptPath$toolUsed $phenData $vcffile $outFolder



# #>>>>>> using vcf2gwas for Adult

# VCF2GWASScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/VCF2GWAS/"
# toolUsed="runVCF2GWAS.sh"
# Compiler="bash"

# # for the field
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Field/Field_filtered.vcf"
# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Field/Field.converted.csv"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Field"
# mkdir -p $outFolder

# $Compiler $VCF2GWASScriptPath$toolUsed $vcffile $phenData $outFolder "usePCA"

# # for the seedling
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling/Seedling_filtered.vcf"
# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling/Seedling.converted.csv"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling"
# mkdir -p $outFolder

# $Compiler $VCF2GWASScriptPath$toolUsed $vcffile $phenData $outFolder "usePCA"



#>>>>>> plot GWAS results
# PLOTGEMMAGWASScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/GWAS/plotGWASGEMMA/"
# toolUsed="plot_GWAS_GEMMA.r"
# Compiler="Rscript"

# # for the field
# gwasoutFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Field/gwas_out"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Field/plots"
# chrTable="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/CharTable.txt"
# mkdir -p $outFolder

# # create the output folder
# mkdir -p $outFolder
# for f in $gwasoutFolder/*.assoc.txt
# do
#     echo "Analysing $f"
#     $Compiler $PLOTGEMMAGWASScriptPath$toolUsed  $f $(basename $f .assoc.txt)  $outFolder 0.1 "fixed" 2.5 $chrTable
# done

# # combine the log files
# cat $outFolder/*.log > $outFolder/combined.log

# # # for the seedling
# gwasoutFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/gwas_out"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/plots"
# chrTable="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/CharTable.txt"
# mkdir -p $outFolder

# # create the output folder
# mkdir -p $outFolder
# for f in $gwasoutFolder/*.assoc.txt
# do
#     echo "Analysing $f"
#     $Compiler $PLOTGEMMAGWASScriptPath$toolUsed  $f $(basename $f .assoc.txt)  $outFolder 0.1 "fixed" 2.5 $chrTable
# done

# # combine the log files
# cat $outFolder/*.log > $outFolder/combined.log







# # for the seedling
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling/Seedling_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/PCA"
# MetData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling/Seedling_Meta.tsv"
# mkdir -p $outFolder

# $Compiler $GEMMAScriptPath$toolUsed $vcffile $outFolder $outFolder $MetData
















# ## Population structure using LEA
# LEAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="LEApopStructure.r"

# Compiler="Rscript"

# # for the seedling
# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling/Seedling_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/LEA/Run"
# mkdir -p $outFolder

# $Compiler $LEAScriptPath$toolUsed $vcffile $outFolder 12 10000

# # plot the LEA results
# LEAplotScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="plotLEA.r"
# Compiler="Rscript"
# mkdir -p $outFolder/plots


# LEAoutFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/LEA/Run/bestKs"
# resultFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/LEA/plots"
# passportFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Seedling/Seedling_Meta.tsv"
# orderFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Seedling/LEA/Run/TaxaNames.txt"
# bestKnownK=6
# maxK=12
# mkdir -p $resultFolder
# $Compiler $LEAplotScriptPath$toolUsed $LEAoutFolder $resultFolder $passportFile $orderFile $bestKnownK $maxK














#>>>>>> using vcf2gwas for Seedling

# VCF2GWASScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/VCF2GWAS/"
# toolUsed="runVCF2GWAS.sh"
# Compiler="bash"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/Seedling/Seedling_filtered.vcf"
# phenData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/Seedling/SeedlingTraits_Ordered2.csv"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Seedling"
# mkdir -p $outFolder

# $Compiler $VCF2GWASScriptPath$toolUsed $vcffile $phenData $outFolder "usePCA"


# #>>>>>> plot GWAS results

# PLOTGEMMAGWASScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/GWAS/plotGWASGEMMA/"
# toolUsed="plot_GWAS_GEMMA.r"
# Compiler="Rscript"

# gwasoutFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Seedling/gwas_out"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Seedling/plots"

# chrTable="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/ChrTable.txt"
# # create the output folder
# mkdir -p $outFolder

# for f in $gwasoutFolder/*.assoc.txt
# do
#     echo "Analysing $f"
#     #   print(paste("Usage: Rscript plot_GWAS_GEMMA.r <gwasfile> <traitName> <ResultFolder> <fdrcutoff> <thresholdType> <qthreshold> <chrTable>",sep=""))
#     $Compiler $PLOTGEMMAGWASScriptPath$toolUsed  $f $(basename $f .assoc.txt)  $outFolder 0.1 "fixed" 2.5 $chrTable
# done

# # combine the log files
# cat $outFolder/*.log > $outFolder/combined.log



#>>>>>> convert tsv to xlsx
# FormattingBasicFilesTypesPath=$AlsammanFormattingBasicFilesTypesPath"FormattingBasicFilesTypes/"
# toolUsed="tsv_csv2xls_xlsx.sh"
# Compiler="bash"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Seedling/plots/"

# for f in $outFolder/*.tsv
# do
#     echo "Analysing $f"
#     $Compiler $FormattingBasicFilesTypesPath$toolUsed $f $f.xlsx
# done

#>>>>>> convert tsv to xlsx
# FormattingBasicFilesTypesPath=$AlsammanFormattingBasicFilesTypesPath"FormattingBasicFilesTypes/"
# toolUsed="tsv_csv2xls_xlsx.sh"
# Compiler="bash"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Adults/plots/"

# for f in $outFolder/*.tsv
# do
#     echo "Analysing $f"
#     $Compiler $FormattingBasicFilesTypesPath$toolUsed $f $f.xlsx
# done



# #>>>>>> Calculate PCA for Adult
# GEMMAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PCA/"
# toolUsed="PCAbySNPrelate.r"
# Compiler="Rscript"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/Adult/Adult_filtered.vcf"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Adults/PCA"
# MetData="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Adults/PopK_Meta.txt"
# mkdir -p $outFolder

# $Compiler $GEMMAScriptPath$toolUsed $vcffile $outFolder $outFolder $MetData



# #>>>>>> Calculate Kinship for Adult
# GEMMAScriptPath="/home/samman/Documents/MyGitHub/PopulationGenetics/PopulationStructure/"
# toolUsed="KinshipBySNPrelate.r"
# Compiler="Rscript"

# vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/Adult/Adult_filtered.vcf"

# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Adult/Kinship"

# mkdir -p $outFolder

# $Compiler $GEMMAScriptPath$toolUsed $vcffile $outFolder


#>>>>>> convert tsv to xlsx
# FormattingBasicFilesTypesPath=$AlsammanFormattingBasicFilesTypesPath"FormattingBasicFilesTypes/"
# toolUsed="tsv_csv2xls_xlsx.sh"
# Compiler="bash"
# outFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Seedling/gwas_out"

# for f in $outFolder/*.txt
# do
#     echo "Analysing $f"
#     $Compiler $FormattingBasicFilesTypesPath$toolUsed $f $f.xlsx
# done



# #################### Report ####################

### Data Filtering
# Data were filtered using vcftools with the following parameters:
# MAF 0.05 Max missing 50%

### GWAS
# The filtered data were used to run GWAS using VCF2GWAS with the following parameters:
# GEMMA
# Association Tests with Univariate Linear Mixed Models. (LMM)
# PCA was calculated to correct for GWAS as a covariate
# GWAS was filtered by FDR < 0.1 and if it was failed with 0.001 of the quantile of the p-value distribution


### Population Clustering
# PCA was calculated using SNPrelate
# Kinship was calculated using SNPrelate
# Population structure was calculated using LEA for K1 to K12