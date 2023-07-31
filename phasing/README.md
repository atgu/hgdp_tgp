## Overview

This README provides a high-level overview of the steps taken to generate the phased HGDP1KG data.
All the input, notebooks, intermediate, and final files used in the analyses are in the following
bucket: [gs://hgdp-1kg/phasing/]()

## 1. QC and subsetting to one VCF per chromosome

As shown in the [Notebook 1](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb1.ipynb), we applied gnomAD sample, variant, and genotype QC filters. The QC'ed MT
had 159,795,273 variants and 4,120 samples. Then, we
[subsetted](https://github.com/atgu/hgdp_tgp/tree/master/phasing/prepare_data_phasing.py) the MT
by chromosome and wrote out a VCF for each chromosome.


## 2. Relatedness checks, removing outliers, and creating pedigree file
### 2.1 Relatedness checks
To check if there are any parent-child relationships in the data, we used
[PC-Relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) and
[IBD](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.identity_by_descent). We cross-checked
the PC-Relate results with the publicly available
[1KG pedigree file](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt)
and found that all parent-child relationships picked by PC-Relate are reported in the 1KG fam pedigree file.
We then looked at possible duplicates/monozygotic twins. The following 5 samples were flagged as duplicates
and one sample from each pair was removed ['HGDP00758', 'HGDP00472', 'HGDP00452', 'LP6005441-DNA_G02', 'HGDP01279']:

| i        | j | kin    |
|----------|----------|----------|
| HGDP00720  | HGDP00758 | 0.488895  |
| HGDP00981	 | HGDP00472 | 0.451554  |
| HGDP01087    | HGDP00452 | 0.438848  |
| HGDP01092 | LP6005441-DNA_G02 | 0.454102  |
| HGDP01273    | HGDP01279 | 0.441555  |

Three out of the 5 pairs have been flagged as possible duplicate before by [Joanna L Mountain & Ramakrishnan, 2005](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3525116/)

### 2.2 Sample ID mismatches
Relatedness checks were followed by checking if there were no sampleID mismatches between the
gnomAD metadata and 1KG pedigree file for the 1KG samples. The following samples IDs were found to be mismatched:

| 1KG ID     | gnomAD meta ID | final ID |
|------------|----------------|----------|
| NA12546  | NA12546B      | NA12546B |
| NA12830	 | NA12830A      | NA12830A |
| NA18874  | NA18874A      | NA18874A |

### 2.3 Filtering and final pedigree file
After addressing the issues above, the 5 duplicates and 24 PCA outliers from
[Notebook 2](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb2.ipynb) were
removed. The final pedigree file has 4091 samples in total.

## 3. Removing duplicates and PCA outliers from VCF and converting to BCF
Because the VCFs were exported before any filtering, we went back and
[filtered out](https://github.com/atgu/hgdp_tgp/tree/master/phasing/filter_and_convert_to_bcf.py)
the duplicates and PCA outliers

## 4. Applying QC filters as recommended by SHAPEIT5 developers
We then applied extra [QC](https://github.com/atgu/hgdp_tgp/tree/master/phasing/qc.py) filters as [recommended by the SHAPEIT5](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/) developers. Only variants with: (1) HWE>=1e-30, (2) F_MISSING<=0.1, and (3) ExcHet>=0.5 && ExcHet<=1.5 were included.

## 5. Phasing
Phasing was run using SHAPEIT5, following the [recommendations from SHAPEIT5](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#phasing) developers.
Large chunks of [20cM](https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38/20cM) and small chunks of [4cM](https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38/4cM)
were used to phase common and rare variants respectively.

## 6. Imputation performance comparison with 1kGP
We compared  imputation accuracy between using HGDP+1kGP vs. 1kGP NYGC in array (GSA, H3Africa, MEGA) and
low coverage (0.5X, 1.0X, 2.0X, 4.0X) datasets. The 1kGP reference panel used in the comparisons is from [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV) and
was [downloaded using this script](https://github.com/atgu/hgdp_tgp/tree/master/phasing/comparisons/donwload_1kgp.sh). Below are the steps
we followed:

### 6.1. Array dataset
We looked at 3 different arrays: GSA, H3Africa, MEGA, and 3 reference panels: HGDP+1kGP, NYGC 1kGP, and TOPMed. For TOPMed,
we already had imputed data for the 3 arrays, and the only thing we did was
[split the VCF files by chromosome and convert to BCF](https://github.com/atgu/hgdp_tgp/tree/master/phasing/comparisons/split_topmed_by_chrom.py)
to avoid loading the large imputed VCF file for every chromosome. Next, for each array data, we ran [phasing
using SHAPEIT5](https://github.com/atgu/hgdp_tgp/tree/master/phasing/comparisons/impute_array.py) following the [recommendations from here](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_snp_array/#phasing)
and imputation using IMPUTE5 twice, using HGDP+1kGP and NYGC 1kGP as reference panels. 

### 6.2 Low-coverage dataset
For low-coverage data we only looked at HGDP+1kGP vs. NYGC 1kGP and looked at the following coverages: 0.5X, 1.0X, 2.0X, 4.0X.
We used GLIMPSE1 for [phasing+imputation](https://github.com/atgu/hgdp_tgp/tree/master/phasing/comparisons/impute_lowcov.py).

### 6.3 Concordance
To look at how well each reference panel is doing, we looked at the concordance between the imputed sets (estimation) using the references
mentioned above and the original 30X NeuroGAP set (truth) using the [GLIMPSE2_concordance](https://odelaneau.github.io/GLIMPSE/docs/documentation/concordance/) module.
Lastly, for each MAF bin, we computed the average r<sup>2</sup> (across chr1-22).
