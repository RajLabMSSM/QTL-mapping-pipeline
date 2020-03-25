ml bcftools
ml tabix
# creates a small semi-QC'd VCF file for the 827 test genes on chr1
# removes indels and multiallelic sites
#VCF=/sc/orga/projects/als-omics/raw_data/Project_CGND_Freeze01_311samples_GRM_WGS/jgwd/joint_vcf/CGND_311JG_GRM_WGS_2019-06-19_chr1.recalibrated_variants.vcf.gz
#REGION="chr1:108181178-112364535"

# use regions.bed as regions file

test_vcf1=/sc/hydra/projects/ad-omics/data/references/QTL-mapping-pipeline-references/CGND_311JG_GRM_WGS_2019-06-19_chr1.recalibrated_variants.vcf.gz
test_vcf2=/sc/hydra/projects/ad-omics/data/references/QTL-mapping-pipeline-references/CGND_311JG_GRM_WGS_2019-06-19_chr2.recalibrated_variants.vcf.gz


bcftools view -R regions.bed -f PASS --min-af 0.05 --exclude-types indels --exclude-types indels $test_vcf1 | bgzip -c > test_chr1.vcf.gz

tabix test_chr1.vcf.gz

bcftools view -R regions.bed -f PASS --min-af 0.05 --exclude-types indels --exclude-types indels $test_vcf2 | bgzip -c > test_chr2.vcf.gz

tabix test_chr2.vcf.gz

bcftools concat -n -o test_all_chr.vcf.gz test_chr1.vcf.gz test_chr2.vcf.gz

tabix test_all_chr.vcf.gz 
