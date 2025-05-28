#GWAS

input:

#Example dataset (Hepatitis C infection)
gwas.bed = binary file with genotypes
gwas.fam = information about individuals
gwas.bim = information about SNP markers

QC (SAMPLES)

1. Check genotyping rate per individual
In general, 1 - 5% of missing data per individual and SNP are allowed

#!/bin/bash
ml PLINK/1.9b_6.21-x86_64 

# START THE APPLICATION
plink --bfile HA --missing --out sampleQC

-----------------------------------------------------------------------------------------------
2. Check individual heterozygosity rate
Less heterozygosity than expected : possible inbreeding
Higher than expected : DNA contamination

# START THE APPLICATION
plink --bfile HA --het --out sampleQC

----------------------------------------------------------------------------------------------

3. Discordant sex 

# START THE APPLICATION
plink --bfile HA --check-sex --out sampleQC

--In our samples: ** --requires at least one polymorphic x chromosome locus 
-----------------------------------------------------------------------------------------------

4. Sample Relatedness/Duaplicated samples

#Create a set of ‘independent’ SNPs (i.e. not in linkage disequilibrium)

# START THE APPLICATION
plink --bfile HA  --indep-pairwise 50 5 0.2 --out LD_HA
plink --bfile HA --extract LD_HA.prune.in --make-bed --out HApruned

#Estimate pairwise genetic identity coefficient (Generate pairwise IBS estimates)

# START THE APPLICATION
plink --bfile HApruned --genome --min 0.125 --out Related (Generates pairwise IBS estimates/ --genome= to exclude 2nd degree relatedness)

#Then remove samples from each pairs of related individuals 

# START THE APPLICATION
plink --bfile HAinfection --remove dup.txt --make-bed --out HAfullQC (create the .txt file with the samples -FID and IID)

#Remove bad samples
plink --bfile HAinfection --remove bad_samples.txt --make-bed --out HAsamplesQC (create the .txt file with the samples -FID and IID)

-------------------------------------------------------------------------------------------------------------------------------------------

QC (SNPs)

1. Genotyping quality: call rate per SNP

# START THE APPLICATION
plink --bfile HAsamplesQC --missing -out SNPQC

2. Differential missingness

# START THE APPLICATION
plink --bfile HAsamplesQC --test-missing -out SNPQC

---------AWK one-liner to see if any SNPs reach cutoff p-value------
awk ' {print $5} ' hapmapEA_sampleQC.missing | sort -g | head -n20

# START THE APPLICATION
plink --bfile HAsamplesQC --hardy -out SNPQC
plink --bfile HAsampleQC --freq -out SNPQC (writes a MAF report)
plink --bfile HAsampleQC --maf 0.05 --exclude badSNPs.txt --make-bed --out HAfullQC (filter out SNPs with MAF lower than 0.05)
plink --bfile HAfullQC --remove dup_samples.txt --make-bed --out HAfullQCnodup
plink --bfile HAfullQC --logistic --out HAnocov
plink --bfile HAfullQC --pca header tabs --out PCA
plink --bfile HAfullQC --covar PCA.eigenvec --covar-name PC1, PC2, PC3, PC4, PC5 --logistic --out HAPC
