
### SMuRF vignette
by [Huang Weitai](https://www.researchgate.net/profile/Weitai_Huang) 

13th Nov 2018

#### <br/>Introduction

SMuRF is an R package that contains functions for the prediction of a consensus set of somatic mutation calls based on a Random Forest machine learning approach. SMuRF generates a set of point mutations and insertions/deletions (indels) trained based on the latest community-curated tumor whole genome sequencing data. Our method is fast and accurate that could be applied to data from different cancer types as well as whole genome or exome data. 

For more information see our BioRxiv preprint doi: https://doi.org/10.1101/270413   

#### <br/>Table of contents

[1. Input from bcbio-nextgen pipeline](#input-bcbio)
</br>[1a. Input directly from VCF Callers (optional)](#input-alt)
</br>[2. Test Dataset](#test)
</br>[3. Requirements for package](#requirements)
</br>[4. Installation instructions](#installation)
</br>[5. Output file description/legend](#output)
</br>[6. Extracting Gene Annotations for somatic mutations in the coding transcripts](#annotation)
</br>[7. Running on multiple samples](#multiple-samples)


<a name="input-bcbio"></a>

#### <br/>1. Input from bcbio-nextgen pipeline

Before you run SMuRF, you will need output data from the [bcbio-nextgen pipeline](http://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#cancer-variant-calling) containing the VCF output for algorithms MuTect2, FreeBayes, VarDict and VarScan. Note that the tabix (.tbi) files will be required for retrieving gene annotations in SMuRF. We would recommend the bcbio-nextgen pipeline for a better user experience.  

<a name="input-alt"></a>

#### <br/>1a. Input directly from VCF Callers (optional)

**For Users not running bcbio-nextgen pipeline:**
Alternatively, you may run these callers individually. The VCF outputs from each caller (.vcf.gz) is required for SMuRF to run. Refer to the installation and instructions for each caller:
<br/>- [VarDict](https://github.com/AstraZeneca-NGS/VarDict)
<br/>- [VarScan](https://github.com/dkoboldt/varscan)
<br/>- [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php)
<br/>- [FreeBayes](https://github.com/ekg/freebayes)

<a name="test"></a>

#### <br/>2. Test Dataset

In this vignette, we will be using a [partial output dataset](https://github.com/skandlab/SMuRF/tree/master/test) derived from the chronic lymphocytic leukemia (CLL) data downloaded from the European Genome-phenome Archive (EGA) under the accession number EGAS00001001539. You may download the test set for testing SMuRF's functions.

<a name="requirements"></a>

#### <br/>3. Requirements for package

Dependencies: R >=3.3.1
Packages: data.table
          VariantAnnotation
          h2o 3.10.3.3 (must be this version)
_These packages will be installed the first time you run SMuRF._          

<a name="installation"></a>

#### <br/>4. Installation instructions

1. The latest version of the package is updated on Github https://github.com/skandlab/SMuRF
<br/>2. You can install the current SMuRF directly from Github via the following R command: 
```r
#devtools is required
install.packages("devtools")
library(devtools)
install_github("skandlab/SMuRF", subdir="smurf")
```
(*Alternative option*) SMuRF installation via downloading of the package from Github: 
```r
#Clone or download package from Github https://github.com/skandlab/SMuRF/tree/master/smurf and save to your local directory
install.packages("my/current/directory/smurf", repos = NULL, type = "source")
```

Download the test files into your designated directory: https://github.com/skandlab/SMuRF/tree/master/test 
```r
download.file('https://github.com/skandlab/SMuRF/raw/master/test/varscan.vcf.gz','varscan.vcf.gz')
download.file('https://github.com/skandlab/SMuRF/raw/master/test/vardict.vcf.gz','vardict.vcf.gz')
download.file('https://github.com/skandlab/SMuRF/raw/master/test/mutect2.vcf.gz','mutect2.vcf.gz')
download.file('https://github.com/skandlab/SMuRF/raw/master/test/freebayes.vcf.gz','freebayes.vcf.gz')
```

Before we start using the package's functions, set your designated file directory containing your sample/test files
```r
mydir <- getwd() #get current directory
# alternative option
# mydir <- setwd("my/local/directory/for/test/files")

```
_SMuRF_ will predict both single somatic nucleotide variants (SNV) as well as small insertions and deletions (indels). In this example we will be using the "combined" import functionality.
```r
library("smurf") #load SMuRF package
myresults <- smurf(mydir, "combined") #save output into 'myresults' variable

#this will run SMuRF and generate predictions based on input files in 'mydir'

```
_The first time you run SMuRF, required packages that are missing will be installed._

Output files saved includes:
1. Variant statistics (_stats_) 
2. Predicted reads (_predicted_)
3. Parsed-raw file (_parse_)
4. Predicted reads with annotations (_annotated_)* #for smurf's "annotation" function only
5. Time taken (_time_)

```r
# myresults <- smurf(mydir, "combined") 

myresults$smurf_indel$stats_indel

<!--             Passed_Calls -->
<!-- Mutect2             1546 -->
<!-- FreeBayes            339 -->
<!-- VarDict              515 -->
<!-- VarScan             2228 -->
<!-- Atleast1            4343 -->
<!-- Atleast2             244 -->
<!-- Atleast3              37 -->
<!-- All4                   4 -->
<!-- SMuRF_INDEL            4 -->

myresults$smurf_indel$predicted_indel

<!-- Chr START_POS_REF END_POS_REF REF ALT  REF_MFVdVs ALT_MFVdVs FILTER_Mutect2 FILTER_Freebayes FILTER_Vardict FILTER_Varscan -->
<!--   1      17820432    17820433  AT   A AT/AT/AT/AT    A/A/A/A           TRUE            FALSE           TRUE           TRUE -->
<!--   1      81654021    81654022  CA   C CA/CA/CA/CA    C/C/C/C           TRUE             TRUE           TRUE           TRUE -->
<!--   1      91134042    91134043  CT   C CT/CT/CT/CT    C/C/C/C           TRUE             TRUE           TRUE           TRUE -->
<!--   1      32639063    32639064  TA   T TA/TA/NA/TA   T/T/NA/T          FALSE             TRUE          FALSE           TRUE -->
   
<!-- Sample_Name Alt_Allele_Freq N_refDepth N_altDepth T_refDepth T_altDepth SMuRF_score -->
<!--    icgc_cll           0.565         24          1         12         13   0.8727273 -->
<!--    icgc_cll           0.464         26          0         15         13   0.5454545 -->
<!--    icgc_cll           0.485         32          0         17         16   0.7090909 -->
<!--    icgc_cll           0.381         21          0         13          8   0.6701292 -->
     
myresults$smurf_snv$stats_snv

<!--           Passed_Calls -->
<!-- Mutect2           4906 -->
<!-- FreeBayes          247 -->
<!-- VarDict            315 -->
<!-- VarScan           5101 -->
<!-- Atleast1         10302 -->
<!-- Atleast2           170 -->
<!-- Atleast3            58 -->
<!-- All4                39 -->
<!-- SMuRF_SNV          417 -->

head(myresults$smurf_snv$predicted_snv)

<!-- Chr START_POS_REF END_POS_REF REF ALT REF_MFVdVs ALT_MFVdVs FILTER_Mutect2 FILTER_Freebayes FILTER_Vardict FILTER_Varscan -->
<!--   1       5035185     5035185   C   T    C/C/C/C    T/T/T/T           TRUE             TRUE           TRUE           TRUE -->
<!--   1       8929624     8929624   A   G    A/A/A/A    G/G/G/G           TRUE             TRUE           TRUE           TRUE -->
<!--   1      11398873    11398873   T   C  T/NA/NA/T  C/NA/NA/C           TRUE            FALSE          FALSE           TRUE -->
<!--   1      12207135    12207135   G   A    G/G/G/G    A/A/A/A           TRUE             TRUE           TRUE           TRUE -->
<!--   1      14955425    14955425   C   A    C/C/C/C    A/A/A/A           TRUE             TRUE           TRUE           TRUE -->
<!--   1      22385823    22385823   A   G    A/A/A/A    G/G/G/G           TRUE             TRUE           TRUE           TRUE -->
    
<!--  Sample_Name Alt_Allele_Freq N_refDepth N_altDepth T_refDepth T_altDepth SMuRF_score -->
<!-- icgc_cll           0.309         86          2         42         18   0.9861111 -->
<!-- icgc_cll           0.375         60          1         13          8   0.9722222 -->
<!-- icgc_cll           0.147         76          1         59         10   0.9136776 -->
<!-- icgc_cll           0.250         74          1         37         14   0.9861111 -->
<!-- icgc_cll           0.385         33          1         12          7   0.9305556 -->
<!-- icgc_cll           0.559         40          2         16         20   1.0000000 -->
       
```

<a name="output"></a>

#### </br>5. Output file description/legend

Column Name | Description
----------- | -------------------------------------------------------------------------------------------------
Chr         | Chromosome number
START_POS_REF/END_POS_REF         | Start and End nucleotide position of the somatic mutation
REF/ALT     | Consensus Ref and Alt nucleotide changes of the highest likelihood
REF_MFVdVs/ALT_MFVdVs        | Reference and Alternative nucleotide changes from each caller; Mutect2 (M), Freebayes (F), Vardict (Vd), Varscan (Vs)
FILTER | Passed (TRUE) or Reject (FALSE) [boolean] mutation calls from the individual callers
Sample_Name | Sample name is extracted based on your labeled samples in the vcf files
Alt_Allele_Freq | Mean Variant allele frequency calculated from the tumor reads of the callers
Depth ref/alt N/T | Mean read depth from the N/T sample for ref/alt alleles
SMuRF_score      | SMuRF confidence score of the predicted mutation


</br>You may also retrieve the time taken for your run.
```r
myresults$time.taken

<!-- Time difference of 20.52405 secs -->
```

You can check the parsed output used for the prediction:
```r
myresults$smurf_indel$parse_indel

myresults$smurf_snv$parse_snv

```
Proceed to save the output in your desired formats
Example:
```r
a<- myresults$smurf_indel$stats_indel
write.table(a , file = "indel-stats.txt", sep = "\t", quote = FALSE, row.names = TRUE, na = ".")

a<- myresults$smurf_indel$predicted_indel
write.table(a , file = "indel-predicted.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

a<- myresults$smurf_indel$parse_indel
write.table(a , file = "indel-parse.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

a<- myresults$smurf_snv$stats_snv
write.table(a , file = "snv-stats.txt", sep = "\t", quote = FALSE, row.names = TRUE, na = ".")

a<- myresults$smurf_snv$predicted_snv
write.table(a , file = "snv-predicted.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

a<- myresults$smurf_snv$parse_snv
write.table(a , file = "snv-parse.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

a<- myresults$smurf_snv_annotation$annotated
if(!is.null(a)){
  write.table(a , file = "snv-annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
}

a<- myresults$smurf_indel_annotation$annotated
if(!is.null(a)){
  write.table(a , file = "indel-annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
}

a<- myresults$time.taken
write(a, file = "time.txt")
```
<a name="annotation"></a>

#### <br/>6. Extracting Gene Annotations for somatic mutations in the coding transcripts

As an additional function, we have added gene annotations using SnpEff (from bcbio) and SMuRF extracts the coding annotations from the canonical transcripts with the highest impact for your convenience. Note that your vcf.gz files should ready been tab-indexed (.tbi files required).
```r
myresults <- smurf(mydir, "cdsannotation") #runs SMuRF for SNV and indels + generate annotations

```

You may check the output files generated by the test samples in this section to the expected results we provided located in the _results_ folder https://github.com/skandlab/SMuRF/tree/master/test/results.


<a name="multiple-samples"></a>

#### <br/>7. Running on multiple samples

Use our R package to efficiently do somatic mutation predictions on multiple matched tumor-normal samples by providing the list of directories of where your sample files are located. 
```r
#Example

sample_directories <- list("my/dir/sample_A", "my/dir/sample_B", "my/dir/sample_C")

myresults <- list()

for(i in 1:length(sample_directories))
 {
 myresults[[i]] <- smurf(sample_directories[i], "combined")
 } 
 
#myresults[[1]]$time.taken
#Time difference of 9.973997 secs

#myresults[[2]]$time.taken
#Time difference of 11.1712 secs

#myresults[[3]]$time.taken
#Time difference of 15.18325 secs
```

<br/>
For errors and bugs, please report on our Github page.

  
  
  
  
  
  
  

