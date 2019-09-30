### SMuRF vignette
by [Huang Weitai](https://www.researchgate.net/profile/Weitai_Huang) 

24th Sept 2019

#### <br/>Introduction

_SMuRF_ R package predicts a consensus set of somatic mutation calls using RandomForest machine learning. SMuRF generates a set of point mutations and insertions/deletions (indels) trained on the latest community-curated tumor whole genome sequencing data (Alioto _et. al._, 2015, Nat. Comms.). Our method is fast and accurate and analyses both whole-genome and whole-exome sequencing data from different cancer types. 

For more information see our Bioinformatics paper: https://doi.org/10.1093/bioinformatics/btz018   

**Citation** 
</br>Huang, W., et al., SMuRF: Portable and accurate ensemble prediction of somatic mutations. Bioinformatics, 2019: p. 3157-3159

#### <br/>Table of contents
[Input from bcbio-nextgen pipeline](#input-bcbio)
</br>[Input directly from VCF Callers (optional)](#input-alt)
</br>[Test Dataset](#test)
</br>[Requirements](#requirements)
</br>[Installation](#installation)
</br>[Functions](#functions)
</br>-Get SNV and indel predictions
</br>-Annotate genes
</br>-Subset region-of-interest
</br>-Change prediction cutoffs
</br>[Output format](#output)
</br>[Running on multiple samples](#multiple-samples)

____________________________________________________

<a name="input-bcbio"></a>

#### <br/>Input from bcbio-nextgen pipeline

Before running _SMuRF_, you require output data from the [bcbio-nextgen pipeline](http://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#cancer-variant-calling) that generates the VCF output for the variant callers: MuTect2, FreeBayes, VarDict and VarScan. Note that your vcf.gz files need to be tab-indexed (.tbi files required) for retrieving gene annotations in SMuRF. We would recommend the bcbio-nextgen pipeline for a better user experience.  

_SMuRF_ requires the VCF output from each caller (.vcf.gz) to be placed in the same directory and files tagged with the caller (eg. sample1-mutect.vcf.gz, sample1-freebayes.vcf.gz, sample1-vardict.vcf.gz, sample1-varscan.vcf.gz)


<a name="input-alt"></a>

#### <br/>Input directly from VCF Callers (optional)

**For Users not running bcbio-nextgen pipeline:**
Alternatively, install and execute the individual callers. 

Refer to the installation and instructions for each caller:
<br/>- [VarDict](https://github.com/AstraZeneca-NGS/VarDict)
<br/>- [VarScan](https://github.com/dkoboldt/varscan)
<br/>- [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php)
<br/>- [FreeBayes](https://github.com/ekg/freebayes)


<a name="test"></a>

#### <br/>Test Dataset

In this vignette, we utilise a [partial output dataset](https://github.com/skandlab/SMuRF/tree/master/test) derived from the chronic lymphocytic leukemia (CLL) data downloaded from the European Genome-phenome Archive (EGA) under the accession number EGAS00001001539. This test dataset is provided in the SMuRF package.


<a name="requirements"></a>

#### <br/>Requirements

| R 3.3 & 3.4 | R >=3.5 |
|:-:|:-:|
|Java 7 or 8|Java 7 or higher (tested on Java 11)|
|h2o-3.10.3.3 (included in package)|h2o-3.26.0.2 (requires download*)|
|bioconductor::VariantAnnotation|BiocManager::VariantAnnotation|

_If h2o-3.26.0.2 takes a long time to install, try manually installing from their [AWS page.](https://h2o-release.s3.amazonaws.com/h2o/rel-yau/2/index.html)_

<a name="installation"></a>

#### <br/>Installation

<br/>1. The latest version of the package is updated on Github https://github.com/skandlab/SMuRF

<br/>2. You can install the current SMuRF directly from Github via the following R commands: 
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

_SMuRF_ concurrently predicts single somatic nucleotide variants (SNV) as well as small insertions and deletions (indels) and saves time by parsing the VCF files once. 
```r
library("smurf") #load SMuRF package

smurf() #check version and parameters

# "SMuRFv1.6 (9th Sept 2019)"
# smurf(directory=NULL, model=NULL, nthreads = -1,
                 save.files=F,  output.dir=NULL, parse.dir=NULL, roi.dir=NULL,
                 snv.cutoff = 'default', indel.cutoff = 'default',
                 t.label='-T', build='hg19', change.build=F,
                 check.packages=T)
myresults = smurf(mydir, "combined") #save output into 'myresults' variable

```

<a name="functions"></a>

#### <br/>Functions

|model|Description|
|:-----|:--------------------------------------------------------|
|combined|SNVs + indels|
|cdsannotation|SNVs + indels + Gene annotations|
|newcutoff.cdsannotation|Re-adjusts cutoffs for (SNVs + indels) + Gene annotations|
|ROI|Analyse a region based on a set of genomic coordinates (.BED or .txt)|
|featureselectionall| _Debug mode_: retrieve parse file without _SMuRF_ predictions|

_Missing packages will be installed the first time you run _SMuRF_._

</br>**_Debug options_**

options|Description
----|--------------------------------
build='hg19'|Specifies genomic build='hg38' for correct gene annotations
change.build=F|Choose the conversion from 'hg19.to.hg38' or 'hg38.to.hg19' for your snv- and indel- annotated output
check.packages=T|Developer mode

For more information on the parameters see R documentation:
```r
help(smurf)
```
</br>
Examples:
```r
library("smurf") #load SMuRF package

myresults = smurf(directory = paste0(find.package("smurf"), "/data"), #test dataset
                  model="combined") #model
                  
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  model="cdsannotation",
                  nthreads=1, #no. of cores
                  save.files = T, #save output as .txt files
                  output.dir = 'C:/Users/admin/myresults') #output dir
```

</br>_SMuRF_ is fine-tuned to achieve higher sensitivity. Re-adjust the stringency of the prediction with a specific cut-off value. 
</br>Use parameters 'snv.cutoff' or 'indel.cutoff' to adjust the thresholds.
</br>To re-adjust the cut-off value of an **existing** SMuRF output, simply provide the path to the snv-parse and indel-parse files for immediate re-processing. 

```r
#start with a specific cutoff
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  model="combined", #change cutoffs
                  snv.cutoff=0.2, indel.cutoff=0.1, #specify new cutoffs
                  save.files = T, 
                  output.dir = 'C:/Users/admin/myresults') 

#modify cutoff from existing SMuRF parse files
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  model="newcutoff.cdsannotation", #change cutoffs only
                  snv.cutoff=0.2, indel.cutoff=0.1,
                  parse.dir = 'C:/Users/admin/myresults', #SMuRF output path to existing parse.txt files
                  save.files = T, 
                  output.dir = 'C:/Users/admin/myresults2') 

#Plot histogram
hist(as.numeric(myresults$smurf_indel$predicted_indel[,'SMuRF_score']), main = 'Re-adjusted predicted indels', xlab = 'SMuRF_score', col = 'grey50')
```

</br>SMuRF can extract variant calls from a BED or (.txt) file containing the genomic coordinates of your region of interest.
Use genomic coordinates referencing hg19/GRCh37. Specify the BED file path using 'roi.dir'

```r
# roi.df = read.delim(roi.dir, header = F)
# roi.df
#   chrom    start      end
# 1     1        1    25000
# 2     1     1000   100000
# 3     1 77000000 78000000

myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  model="ROI", #retrieve SNVs + indels from regions-of-interest
                  roi.dir = paste0(find.package("smurf"), "/data/roi.bed"), #BED file containing ROIs
                  save.files = T, 
                  output.dir = 'C:/Users/admin/myresults3')
```


<a name="output"></a>

#### </br>Output format

Output files saved include:

1. Variant statistics (_stats_) 

2. Predicted reads (_predicted_)

3. Parsed-raw file (_parse_)

4. Predicted reads with annotations (_annotated_)* #for smurf's "cdsannotation" function only

5. Time taken (_time_)

</br>
```r
#Main predicted file (SNV & indel)

myresults$smurf_snv$predicted_snv

myresults$smurf_indel$predicted_indel

```

Column | Description
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


</br>
```r
myresults$smurf_indel$stats_indel

#             Passed_Calls
# Mutect2             1546
# FreeBayes            339
# VarDict              515
# VarScan             2228
# Atleast1            4343
# Atleast2             244
# Atleast3              37
# All4                   4
# SMuRF_INDEL            6


myresults$smurf_snv$stats_snv

#           Passed_Calls
# Mutect2           4906
# FreeBayes          247
# VarDict            315
# VarScan           5101
# Atleast1         10302
# Atleast2           170
# Atleast3            58
# All4                39
# SMuRF_SNV         1345

```

</br>We added gene annotations using SnpEff (from bcbio) and _SMuRF_ extracts the coding annotations from the canonical transcripts with the highest impact. Take note that your vcf.gz files should be tab-indexed (.tbi files required).
```r
myresults = smurf(mydir, "cdsannotation") #runs SMuRF for SNV and indels + generate annotations

myresults$smurf_snv_annotation$annotated
#   Chr START_POS_REF END_POS_REF REF ALT REF_MFVdVs ALT_MFVdVs FILTER_Mutect2 FILTER_Freebayes FILTER_Vardict FILTER_Varscan Sample_Name Alt_Allele_Freq
# 1   1      77712621    77712621   G   A    G/G/G/G    A/A/A/A           TRUE             TRUE           TRUE           TRUE    icgc_cll           0.296
# 2   1      77806132    77806132   G   A    G/G/G/G    A/A/A/A           TRUE             TRUE           TRUE           TRUE    icgc_cll           0.483

#   T_altDepth T_refDepth N_refDepth N_altDepth Allele       Annotation   Impact Gene_name         Gene_ID Feature_Type      Feature_ID Transcript_BioType
# 1          8         19         40          0   <NA>             <NA>     <NA>      <NA>            <NA>         <NA>            <NA>               <NA>
# 2         15         16         22          0      A missense_variant MODERATE       AK5 ENSG00000154027   transcript ENST00000354567     protein_coding

#   Rank   HGVS.c      HGVS.p  cDNA.pos  CDS.pos  AA.pos Distance     REGION SMuRF_score
# 1 <NA>     <NA>        <NA>      <NA>     <NA>    <NA>     <NA> Non-coding   0.8148148
# 2 6/14 c.770G>A p.Arg257His 1033/3248 770/1689 257/562        .        CDS   0.7777778

```

</br>Time taken for your run:
```r
myresults$time.taken

<!-- Time difference of 20.52405 secs -->
```

</br>The raw parsed output:
```r
myresults$smurf_indel$parse_indel

myresults$smurf_snv$parse_snv

```
</br>Proceed to save the output in your desired formats
Use 'save.files=T' to generate .txt outputs in your output directory
Example:
```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                           model ="cdsannotation", nthreads = 1,
                           output.dir = 'C:/Users/admin/myresults', save.files = T)

```

<a name="multiple-samples"></a>

#### <br/>Running on multiple samples

Iterate over multiple samples by providing the list of directories of where your sample files are located. 

```r
#Example

project.dir = 'path/to/my/dir'
samples = c('sample_A', 'sample_B', 'sample_C') #sample dir where vcf files are located

for(i in 1:length(samples)) {
 smurf(directory=paste0(project.dir, '/', samples[i]),
        model="cdsannotation",
        output.dir = paste0('C:/Users/admin/myresults-',sample[i]), save.files = T))
 } 
```

<br/>
For errors and bugs, please report on our Github page.

  
  
  
  
  
  
  

