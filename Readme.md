# SMuRF v3.0
By [Skandlab](https://github.com/skandlab) 

Genome Institute of Singapore, A*STAR

Check out the [latest SMuRF version here](https://github.com/skandlab/SMuRF/releases)

<a name="home"></a>

#### <br/>Introduction

_SMuRF_ R package predicts a consensus set of somatic mutation calls using RandomForest machine learning. _SMuRF_ generates a set of point mutations and insertions/deletions (indels) trained on the latest community-curated tumor whole genome sequencing data (Alioto _et. al._, 2015, Nat. Comms.). Our method is fast and accurate and analyses both whole-genome and whole-exome sequencing data from different cancer types. 

For more information see our Bioinformatics paper: https://doi.org/10.1093/bioinformatics/btz018   

**Citation** 
</br>Huang W, Guo YA, Chang MM and Skanderup AJ. Ensemble-Based Somatic Mutation Calling in Cancer Genomes. In: Boegel S, editor. Bioinformatics for Cancer Immunotherapy: Methods and Protocols. New York, NY: Springer US; 2020. p. 37-46.

Huang W, Guo YA, Muthukumar K, Baruah P, Chang MM and Skanderup AJ. SMuRF: Portable and accurate ensemble prediction of somatic mutations. Bioinformatics (Oxford, England). 2019:btz018-btz. doi:10.1093/bioinformatics/btz018.

#### <br/>Table of contents
[Input from bcbio-nextgen pipeline](#input-bcbio)
</br>[Input directly from VCF Callers (optional)](#input-alt)
</br>[Test Dataset](#test)
</br>[Requirements](#requirements)
</br>[Installation](#installation)
</br>[Parameters](#functions)
</br>[Running SMuRF: Selecting the correct input vcfs](#input)
</br>[Running SMuRF: Detecting and changing genome build](#build)
</br>[Running SMuRF: Tweaking SMuRF score cut-off](#cutoff)
</br>[Output format](#output)
</br>[Running on multiple samples](#multiple-samples)

____________________________________________________

<a name="input-bcbio"></a>

#### <br/>Input from bcbio-nextgen pipeline

Before running _SMuRF_, you require output data from the [bcbio-nextgen pipeline](http://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#cancer-variant-calling) that generates the VCF output for the variant callers: MuTect2, FreeBayes, VarDict, VarScan and the latest Strelka2. An additional caller Strelka2, has been added since SMuRF 2.0  and the information is documented on our [wiki page](https://github.com/skandlab/SMuRF/wiki/SMuRF-3.0). 

SMuRF v1.6.4 is still available here: [SMuRFv1.6.4](https://github.com/skandlab/SMuRF/releases/tag/SMuRFv1.6.4)
</br>SMuRF v1.6.4 wiki page: [readme file](https://github.com/skandlab/SMuRF/wiki/SMuRF-v1.6.4-vignette)

Note that your vcf.gz files need to be tab-indexed (.tbi files required) for retrieving gene annotations in SMuRF. We would recommend the bcbio-nextgen pipeline for a better user experience.  See [Running SMuRF: Selecting the correct input vcfs](#input) for more information.

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
<br/>- [Strelka2](https://github.com/Illumina/strelka)


<a name="test"></a>

#### <br/>Test Dataset

In this vignette, we utilise a [partial output dataset](https://github.com/skandlab/SMuRF/tree/master/test) derived from the chronic lymphocytic leukemia (CLL) data downloaded from the European Genome-phenome Archive (EGA) under the accession number EGAS00001001539. The dataset for testing the package is provided in the SMuRF package.


<a name="requirements"></a>

#### <br/>Requirements

**R 3.3 & 3.4** : bioconductor::VariantAnnotation

**R >=3.5** : BiocManager::VariantAnnotation

**h2o package** : 
_If h2o package takes some time to download/install (~350MB), try manually installing from their [AWS page.](https://h2o-release.s3.amazonaws.com/h2o/rel-yau/2/index.html)_

<a name="installation"></a>

#### <br/>Installation

<br/>1. The latest version of the package is updated on Github https://github.com/skandlab/SMuRF

2. You can install the current SMuRF directly from Github via the following R commands: 
```r
#devtools is required
install.packages("devtools")
library(devtools)
install_github("skandlab/SMuRF", subdir="smurf")
```

</br>(*Alternative option*) SMuRF installation via downloading of the package from Github: 
```r
#Clone or download package from Github https://github.com/skandlab/SMuRF/tree/master/smurf and save to your local directory
install.packages("my/current/directory/smurf", repos = NULL, type = "source")
```

</br> _SMuRF_ concurrently predicts single somatic nucleotide variants (SNV) as well as small insertions and deletions (indels) and saves time by parsing the VCF files once. 

_Missing packages will be installed the first time you run _SMuRF_._

```r
library("smurf") #load SMuRF package

smurf() #check version and parameters

# "SMuRFv3.0 (16th Jan 2024)"
smurf(directory=NULL, mode=NULL, nthreads = -1,
                 annotation=F, output.dir=NULL,  parse.dir=NULL,
                 snv.cutoff = 'default', indel.cutoff = 'default',
                 build=NULL, change.build=F, find.build=F,
                 t.label=NULL, re.tabIndex=F,
                 check.packages=T, file.exclude=NULL)

myresults = smurf(mydir, 'combined', build='hg19') #save output into 'myresults' variable

```
[back to top](#home)

<a name="functions"></a>

#### <br/>Parameters

Arguments|Description
:-:|:--------------------------------------------------------
directory|Choose directory where the Variant Caller Format(VCF) files are located
output.dir|Path to output directory (if saving files as .txt)
parse.dir|Specify if changing SMuRF default cutoffs. Path to the location of existing snv-parse.txt and indel-parse.txt files generated by SMuRF
mode|Choose "snv", "indel" or "combined" (snv+indel). "combined" provides a separate list of SNVs and indels.
annotation|TRUE or FALSE (default). Provide gene annotations for each variant call.
nthreads|Number of cores used for RandomForest prediction. Default (-1) for maximum number of cores. _For 32-bit Windows, only 1 core is allowed (nthreads=1)._
t.label|(Optional) Provide the sample name for your tumour sample to ease the identification of the normal and tumour sample names in your vcf
file.exclude|(Optional) Additional keywords in file directory names to be filtered.
build|Specify your human genome build: build="hg19" or build="hg38"
change.build|TRUE or FALSE (default). For conversion of your genomic coordinates
find.build|TRUE or FALSE (default). Additional genome build check for the annotation step.
snv.cutoff|Default SMuRF_score cutoff for the SNV model unless a number between 0 to 1 is stated
indel.cutoff|Default SMuRF_score cutoff for the INDEL model unless a number between 0 to 1 is stated
re.tabIndex|TRUE or FALSE (default). Set to TRUE to create tab-indexed (.tbi) files for each vcf
check.packages=T|Developer mode


For more information on the parameters see R documentation:
```r
help(smurf)
```
</br>
Examples:

```r
library("smurf") #load SMuRF package

 myresults = smurf(directory="/path/to/directory..",
                   mode="snv", #snv only
                   output.dir="/path/to/output", #saving your output
                   build='hg19')
 
 #Include gene annotations for coding regions in output
 myresults = smurf(directory="/path/to/directory..",
                   mode="combined", #snv and indel predictions
                   annotation=T, #generate gene annotations
                   build='hg19')

```

[back to top](#home)

<a name="input"></a>

#### </br>Running SMuRF: Selecting the correct input vcfs

</br>_SMuRF_ requires 5 caller VCF (vcf.gz) files as input stated under the "directory" parameter. Provide a path to a directory containing all 5 caller VCF files. **caller.vcf.gz** (compressed) and **caller.vcf** are accepted formats. 

The tab-indexed (.tbi) files for each caller are required for the parsing step. If the **.tbi** files are missing,  specify using _re.tabIndex=T_ on SMuRF to generate these files.

```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode ="snv", nthreads = 1, annotation = T, build = 'hg19',
                  re.tabIndex = T)  #generate .tbi files in directory
#"Generating .tbi files in directory..."
# Connection successful!

#If the vcf files are in different directories:

#Specify directories manually
 dir.list = list(mutect='/path1/to/mutect.vcf.gz',
                 freebayes='/path2/to/freebayes.vcf.gz',
                 vardict='/path3/to/vardict.vcf.gz',
                 varscan='/path4/to/varscan.vcf.gz',
                 strelka='/path5/to/strelka.vcf.gz')
 myresults = smurf(directory=dir.list, 
                   mode="combined", build='hg19')

```
</br>In some cases, your input directory may contain other VCF files generated by bcbio. For example, germline VCF files, copy-number related files, older version VCFs. An exclusion _file.exclude_ can be added to make sure that SMuRF selects the correct VCF files.

```r
list.files(directory)
# sample1.mutect.vcf.gz
# sample1.mutect-germline.vcf.gz #to be excluded
# sample1.freebayes.vcf.gz
# sample1.vardict.vcf.gz
# sample1.varscan.vcf.gz
# sample1.varscan-version1.vcf.gz #to be excluded
# sample1.strelka.vcf.gz
# sample1.strelka-archive.vcf.gz #to be excluded

myresults = smurf(directory="/path/to/directory..", 
                  file.exclude = c("germline","version1","archive") #keywords in file name to be excluded
                  mode="snv",
                  output.dir="/path/to/output", build='hg19')
```

</br>It is optional to indicate your normal and tumour sample labels. _SMuRF_ detects your normal and tumour sample names in order to generate variant allele frequency information. If this information is missing in your VCF headers, _SMuRF_ will terminate with an error. State your unique tumour sample label using _t.label_.

Possible normal/tumour sample labels:

sample1-N, sample1-T
</br>sample1_normal, sample1_tumour
</br>sample1.healthy, sample1.cancer

```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode ="combined", nthreads = 1, build = 'hg19',
                  t.label = 'tumour' #optional if labels were detected from vcf headers correctly
                  )
```

[back to top](#home)

<a name="build"></a>

#### </br>Running SMuRF: Detecting and changing genome build

</br> The genome build for your sample must be specified ( _build='hg19'_ or _build='hg38'_ ). 

hg19 also refers to the Genome Reference Consortium Human Build 37 (GRCh37) 
</br>hg38 also refers to the Genome Reference Consortium Human Build 38 (GRCh38)

The genome build stated in _SMuRF_ will be cross-checked with the build used in your VCF files. 

```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode ="combined", nthreads = 1, annotation = T, 
                  build = 'hg38' #wrong build stated
                  )
# "Genome build stated in SMuRF:"
# "hg38"
# "Ref genome used in vcf:"
# "file:///home/projects/13001264/softwares/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"
# "Warning: build provided does not match ref genome used in vcf. SMuRF CDS annotation may not run properly if genome build is incorrect."
# "Final genome build used for analysis: hg38"
# 
# Warning message
```

</br>If you are unsure of the genome build used in your analysis, specify _find.build=T_.

```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode ="combined", nthreads = 1, annotation = T, 
                  build = 'hg38', #wrong build stated
                  find.build = T, #if unsure of genome build
                  )
# "Genome build stated in SMuRF:"
# "hg38"
# "Ref genome used in vcf:"
# "file:///home/projects/13001264/softwares/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"
# "Warning: build provided does not match ref genome used in vcf. SMuRF CDS annotation may not run properly if genome build is incorrect."
# "Changing build variable provided"
# "hg38 -> hg19"
# "Final genome build used for analysis: hg19"

# No errors
```

</br>Samples from different batches may be aligned to a different genome reference build. In order to standardize your gene annotations and output, specify _change.build_ for genome build conversion.

```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode ="combined", nthreads = 1, annotation = T, 
                  build = 'hg19',
                  change.build = T, #genome build conversion
                  )
# "Genome build stated in SMuRF:"
# "hg19"
# "Ref genome used in vcf:"
# "file:///home/projects/13001264/softwares/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"
# "Final genome build used for analysis: hg19"

# "Compiling annotations"
# "Changing annotations from hg19 to hg38"
```


[back to top](#home)

<a name="cutoff"></a>

#### </br>Running SMuRF: Tweaking SMuRF score cut-off
</br>_SMuRF_ v3.0 is fine-tuned to achieve the max f1 score in our test set. 

Re-adjust the stringency of the prediction with a specific cut-off value. 
Use parameters _snv.cutoff_ or _indel.cutoff_ to adjust the thresholds (higher cut-off provide a smaller set of calls with better confidence).

To re-adjust the cut-off value of an **existing** SMuRF output, simply provide the _parse.dir_ to the snv-parse and indel-parse files for re-processing. 

```r
#start with default cutoffs
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode="combined", 
                  snv.cutoff='default', indel.cutoff='default',
                  output.dir = 'C:/Users/admin/myresults') 

#modify cutoff from existing SMuRF parse files
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode="combined",
                  snv.cutoff=0.2, indel.cutoff=0.1, #specify new cutoffs
                  parse.dir = 'C:/Users/admin/myresults', #SMuRF path existing parse.txt files
                  output.dir = 'C:/Users/admin/myresults2' #new output) 

#Plot histogram
hist(as.numeric(myresults$smurf_indel$predicted_indel[,'SMuRF_score']), main = 'Re-adjusted predicted indels', xlab = 'SMuRF_score', col = 'grey50')
```
[back to top](#home)

<a name="output"></a>

#### </br>Output format

Output files available include:

1. Parsed-raw file (_parse_)

2. Predicted positive mutations (_predicted_)

3. Predicted positive mutations with annotations (_annotated_)* #for smurf's "cdsannotation" function only

4. Variant statistics (_stats_) 

5. Time taken (_time_)

```r
#Viewing predicted output in R

myresults$smurf_snv$predicted_snv

myresults$smurf_indel$predicted_indel

#see column description below
```

Column | Description
----------- | -------------------------------------------------------------------------------------------------
Chr         | Chromosome number
START_POS_REF/END_POS_REF         | Start and End nucleotide position of the somatic mutation
REF/ALT     | Consensus Ref and Alt nucleotide changes of the highest likelihood
REF_MFVdVs/ALT_MFVdVs        | Reference and Alternative nucleotide changes from each caller; Mutect2 (M), Freebayes (F), Vardict (Vd), Varscan (Vs) and Strelka2 (not abbreviated to preserve column name)
FILTER | Pass (TRUE) or Reject (FALSE) [boolean] mutation calls from the individual callers
Sample_Name | Sample name is extracted based on your labeled samples in the vcf files
Alt_Allele_Freq | Mean Variant allele frequency calculated from the tumor reads of the callers
Depth ref/alt N/T | Mean read depth from the N/T sample for ref/alt alleles
SMuRF_score      | SMuRF confidence score of the predicted mutation


</br>

```r
myresults$smurf_indel$stats_indel

#             Passed_Calls
# Strelka2             466
# Mutect2              232
# FreeBayes            306
# VarDict              483
# VarScan             1273
# Atleast1            2431
# Atleast2             278
# Atleast3              43
# Atleast4               7
# All5                   1
# SMuRF_INDEL           88

myresults$smurf_snv$stats_snv

#           Passed_Calls
# Strelka2          1362
# Mutect2           1539
# FreeBayes          216
# VarDict            239
# VarScan           1734
# Atleast1          4017
# Atleast2           928
# Atleast3            60
# Atleast4            48
# All5                37
# SMuRF_SNV         1043
```

</br>We added gene annotations using SnpEff (from bcbio) and _SMuRF_ extracts the coding annotations from the canonical transcripts with the highest fucntional impact. Take note that your vcf.gz files should be tab-indexed (see [Running SMuRF: re.tabIndex](#input)).

```r
myresults = smurf(mydir, "cdsannotation") #runs SMuRF for SNV and indels + generate annotations

myresults$smurf_snv_annotation$annotated[order(myresults$smurf_snv_annotation$annotated$REGION)[1:2],]
#    Chr START_POS_REF END_POS_REF REF ALT   REF_MFVdVs   ALT_MFVdVs FILTER_Mutect2 FILTER_Freebayes FILTER_Vardict
# 52   1      77806132    77806132   G   A    G/G/G/G/G    A/A/A/A/A           TRUE             TRUE           TRUE
# 81   1     170961432   170961432   C   T C/NA/NA/NA/C T/NA/NA/NA/T           TRUE            FALSE          FALSE
#    FILTER_Varscan FILTER_Strelka2     Sample_Name Alt_Allele_Freq N_refDepth N_altDepth T_refDepth T_altDepth Allele
# 52           TRUE            TRUE icgc_cll_tumour           0.500         14          0         15         15      A
# 81          FALSE            TRUE icgc_cll_tumour           0.467         33          0         16         14      T
#          Annotation   Impact Gene_name         Gene_ID Feature_Type      Feature_ID Transcript_BioType  Rank    HGVS.c
# 52 missense_variant MODERATE       AK5 ENSG00000154027   transcript ENST00000354567     protein_coding  6/14  c.770G>A
# 81 missense_variant MODERATE     MROH9 ENSG00000117501   transcript ENST00000367759     protein_coding 12/22 c.1156C>T
#         HGVS.p  cDNA.pos   CDS.pos  AA.pos Distance REGION SMuRF_score
# 52 p.Arg257His 1033/3248  770/1689 257/562        .    CDS   0.9083840
# 81 p.Arg386Cys 1310/3165 1156/2586 386/861        .    CDS   0.8107475
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

</br>Indicate the _output.dir_ to save the _SMuRF_ output as tab-delimited .txt files in your targeted directory.

```r
myresults = smurf(directory = paste0(find.package("smurf"), "/data"),
                  mode="combined", 
                  output.dir = 'C:/Users/admin/myresults' #path to output directory
                  ) 
```

[back to top](#home)

<a name="multiple-samples"></a>

#### <br/>Running on multiple samples

Iterate over multiple samples by providing the list of directories of where your sample files are located. 

```r
project.dir = 'path/to/my/dir'
list.files(project.dir)
# sample_A
# sample_B
# sample_C
samples = c('sample_A', 'sample_B', 'sample_C') #sample dir where vcf files are located

for(i in 1:length(samples)) {
 smurf(directory=paste0(project.dir, '/', samples[i]),
        mode="combined", build='hg19', annotation = T,
        output.dir = paste0('C:/Users/admin/myresults/',samples[i]))
 } 
```

Running SMuRF on multiple samples on a cluster (parallel multi-core instance)  

```r
install.packages("foreach")
install.packages("doParallel")
install.packages("doSNOW")

library(foreach)
library(doParallel)
library(doSNOW)
library(smurf)

project.dir = 'path/to/my/dir'
samples = Sys.glob(paste0(project.dir,'/*'))

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

foreach(i=1:length(samples), .packages = 'smurf', .verbose = F) %dopar% {
print(i)
  smurf(directory = paste0(project.dir, '/', samples[i]),
      mode ="combined", nthreads = 1, build = 'hg19',
      output.dir = paste0('C:/Users/admin/myresults/',samples[i]))
)
}
stopCluster(cl)
h2o.shutdown()
```

<br/>
For errors and bugs, please report on our Github page.

[back to top](#home)
  
  
  
  
  
  
  

