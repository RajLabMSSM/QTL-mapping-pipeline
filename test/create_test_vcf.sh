ml bcftools
ml tabix
# creates a small semi-QC'd VCF file for the 827 test genes on chr1
# removes indels and multiallelic sites
#VCF=/sc/orga/projects/als-omics/raw_data/Project_CGND_Freeze01_311samples_GRM_WGS/jgwd/joint_vcf/CGND_311JG_GRM_WGS_2019-06-19_chr1.recalibrated_variants.vcf.gz
#REGION="chr1:108181178-112364535"

# use regions.bed as regions file

test_vcf1=/sc/arion/projects/H_ad-omics/data/references/QTL-mapping-pipeline-references/CGND_311JG_GRM_WGS_2019-06-19_chr1.recalibrated_variants.vcf.gz
test_vcf2=/sc/arion/projects/H_ad-omics/data/references/QTL-mapping-pipeline-references/CGND_311JG_GRM_WGS_2019-06-19_chr2.recalibrated_variants.vcf.gz

test_vcf1=/sc/arion/projects/als-omics/WGS_QC/NYGC_Freeze02_European_Feb2020/WGS_QC_Pipeline/NYGC_Freeze02_European_Feb2020/output/chrAll_QCFinished_MAF0.01.anno.vcf.gz
test_vcf2=$test_vcf1

bcftools view -R regions.bed --min-af 0.05 --exclude-types indels --exclude-types indels $test_vcf1 | bcftools sort -Oz -o test_chr1.vcf.gz

tabix test_chr1.vcf.gz

bcftools view -R regions.bed --min-af 0.05 --exclude-types indels --exclude-types indels $test_vcf2 | bcftools sort -Oz -o test_chr2.vcf.gz

tabix test_chr2.vcf.gz

bcftools concat -n test_chr1.vcf.gz test_chr2.vcf.gz | bcftools sort -Oz -o test_all_chr.vcf.gz

tabix test_all_chr.vcf.gz
