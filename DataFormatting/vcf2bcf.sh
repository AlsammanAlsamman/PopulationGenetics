inputfile=$1
outputfile=$2

# convert vcf to gz using bgzip and create index using tabix
bgzip -c $inputfile > $outputfile.gz
tabix -p vcf $outputfile.gz

# convert vcf to bcf using bcftools
bcftools view -O b $outputfile.gz > $outputfile.bcf

# compress bcf using bcftools
gzip $outputfile.bcf

# calculate  the md5sum of the compressed bcf file
md5sum $outputfile.bcf.gz > $outputfile.bcf.gz.md5

# print the size of the compressed bcf file by gigabytes
ourtfileSzie=$(du -h $outputfile.bcf.gz | awk '{print $1}')
echo "The size of the compressed bcf file is $ourtfileSzie"