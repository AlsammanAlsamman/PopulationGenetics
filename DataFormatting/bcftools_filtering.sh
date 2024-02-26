# using bcftools to filter the vcf file

# filter MAF
bcftools filter -i 'MAF[0] > 0.05' GOATs.recode.vcf -o GOATs.recode.MAF.vcf