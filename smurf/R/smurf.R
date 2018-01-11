#' SMuRF v1.3
#'
#' Somatic mutation consensus calling based on four callers:
#' MuTect2, Freebayes, VarDict, VarScan
#' using a RandomForest model to consolidate a list of high accuracy calls.
#'
#' @note
#' Input files containing variant calls should be ".vcf.gz" format of each caller.
#' Supported for R (>=3.3.1)
#' 
#' @param directory Choose directory where the files Variant Caller Format(VCF) files are located.
#'
#' @param model Choose either "snv" or "indel" or "combined"
#' The appropriate parsing and prediction model will be performed
#' to obtain a list of somatic mutation calls. 
#' Choose "combined" to generate both SNV and Indel outputs
#' @examples
#' smurf("/path/to/directory..","combined")
#' smurf("/path/to/directory..","cdsannotation")
#' 
#' @export
smurf = function(directory, model, nthreads = -1, ncores = 1){
  directory <-paste(directory,"/", sep="")
  if(dir.exists(directory)==TRUE){
    

    #check for existing and required packages
    
    
    if("data.table" %in% rownames(installed.packages()) == FALSE){
      install.packages("data.table")
    }
    if("VariantAnnotation" %in% rownames(installed.packages()) == FALSE){
      source("https://bioconductor.org/biocLite.R")
      biocLite("VariantAnnotation")
    }
    
    #check for h2o version 3.10.3.3
    if("h2o" %in% rownames(installed.packages()) == FALSE){
        print("h2o version 3.10.3.3 not found. Installing required h2o package.")
        install.packages("versions")
        library(versions)
        install.versions("h2o", versions = "3.10.3.3")
        print("Running h2o version 3.10.3.3")
  
    } else {
      if((packageVersion("h2o") == '3.10.3.3') == FALSE){
        print("h2o version 3.10.3.3 not found. Installing required h2o package.")
        install.packages("versions")
        library(versions)
        install.versions("h2o", versions = "3.10.3.3")
        print("h2o version 3.10.3.3")
      }  
      #print("h2o version 3.10.3.3 found")
    }
    
    
  
    #load packages
    suppressWarnings(suppressMessages(library(VariantAnnotation)))
    #suppressMessages(library(GenomicRanges))
    
    suppressWarnings(suppressMessages(library(h2o)))
    
    if (exists("nthreads")==FALSE){
      nthreads=-1
    }

    #suppressWarnings(h2o.init(nthreads = -1))
    suppressWarnings(h2o.init(nthreads = nthreads))
    
    
    #Retrieving files from the directory
    
    mutect2 <- paste(directory, Sys.glob("*mutect*.vcf.gz"), sep = "")
    freebayes <- paste(directory, Sys.glob("*freebayes*.vcf.gz"), sep = "")
    varscan <- paste(directory, Sys.glob("*varscan*.vcf.gz"), sep = "")
    vardict <- paste(directory, Sys.glob("*vardict*.vcf.gz"), sep = "")
    
    x<-list(mutect2,freebayes,varscan,vardict)
    
    write("Accessing files:", stderr())
    write("Accessing files:", stdout())
    write(mutect2, stderr())
    write(mutect2, stdout())
    write(freebayes, stderr())
    write(freebayes, stdout())
    write(varscan, stderr())
    write(varscan, stdout())
    write(vardict, stderr())
    write(vardict, stdout())
    
  
    if(length(grep("mutect", mutect2))==1 & length(grep("freebayes", freebayes))==1 & length(grep("varscan", varscan))==1 & length(grep("vardict", vardict))==1){
      
      #Executing smurf if correct parameters stated
      
      if (model == "combined" || model == "cdsannotation" || model == "featureselection" || model == "parseraw" || model == "annotationonly" || model == "parserawmutect") {
          start.time <- Sys.time()
          
          if (model == "combined") {
            parsevcf<-parsevcf_allfeatures(x)
            snvpredict<-snvRFpredict(parsevcf)
            indelpredict<-indelRFpredict(parsevcf)
                      end.time <- Sys.time()
                      time.taken <- end.time - start.time
                      print(time.taken)
                      
                      return(list("smurf_snv"=snvpredict, "smurf_indel"=indelpredict, "time.taken"=time.taken))
          }

          else if (model == "cdsannotation") {
            
            parsevcf<-parsevcf_allfeatures(x)
            
            snvpredict<-snvRFpredict(parsevcf)
            indelpredict<-indelRFpredict(parsevcf)
            
            if (length(snvpredict)>1) {
              print("SNV annotation")
              snvannotation<-CDSannotation(x,snvpredict)
            } else {
              print("No snv predictions. Skipping snv annotation step.")
              snvannotation=NULL
            }
            
            if (length(indelpredict)>1) {
              print("Indel annotation")
              indelannotation<-CDSannotation(x,indelpredict)
            } else {
              print("No indel predictions. Skipping indel annotation step.")
              indelannotation=NULL
            }
            
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            
            return(list("smurf_snv"=snvpredict, "smurf_indel"=indelpredict, "smurf_snv_annotation"=snvannotation, "smurf_indel_annotation"=indelannotation, "time.taken"=time.taken))
 
          }
          
          else if (model == "featureselection") {
            parsevcf<-parsevcf_allfeatures(x)
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            
            return(list("parsevcf_featureselection"=parsevcf, "time.taken"=time.taken))
            
          }
          
          else if (model == "parseraw") { #for Mutectv1
            parsevcf<-parsevcf_raw(x)
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            
            return(list("parsevcf_raw"=parsevcf, "time.taken"=time.taken))
            
          }
          
          else if (model == "parserawmutect") { #for Mutectv1
            parsevcf<-parsevcf_raw_mutect(x)
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            
            return(list("parsevcf_raw"=parsevcf, "time.taken"=time.taken))
            
          }
          
          else if (model == "annotationonly") { #for Mutectv1
            
            suppressWarnings(suppressMessages(library(data.table)))
            
            print("reading snv-raw.txt")
            snvpredict <- fread(paste(directory, "snv-raw.txt", sep = ""), stringsAsFactors=F, skip = 1, drop = 1) #row.names=T
            colnames(snvpredict) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                                      "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                      "Alt_Allele_Freq",
                                      "N_refDepth","N_altDepth","T_refDepth","T_altDepth")
            
            print("reading indel-raw.txt")
            indelpredict <- fread(paste(directory, "indel-raw.txt", sep = ""), stringsAsFactors=F) #row.names=F
            colnames(indelpredict) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                                      "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                      "Alt_Allele_Freq",
                                      "N_refDepth","N_altDepth","T_refDepth","T_altDepth")
            
            #generating vcf subset from bed-granges
            print("merge snv-raw.txt")
            smurfdir <- find.package("smurf")
            beddir <- paste0(smurfdir, "/data/gastric_RF.rds")
            bedfile <- readRDS(beddir)
            sample.name <- snvpredict$Sample_Name[1]
            samplebed <- subset(bedfile, V6==sample.name)
            samplebed <- as.data.frame(samplebed)
            samplebed <- samplebed[,c("V1","V2")]
            colnames(samplebed) <- c("Chr","START_POS_REF")
            snvpredict=merge(samplebed,snvpredict,by=c("Chr","START_POS_REF"),all.x=TRUE)
            
            print("merge indel-raw.txt")
            beddir <- paste0(smurfdir, "/data/gastric_RF_indel.rds")
            bedfile <- readRDS(beddir)
            sample.name <- indelpredict$Sample_Name[1]
            samplebed <- subset(bedfile, V6==sample.name)
            samplebed <- as.data.frame(samplebed)
            samplebed <- samplebed[,c("V1","V2")]
            colnames(samplebed) <- c("Chr","START_POS_REF")
            indelpredict=merge(samplebed,indelpredict,by=c("Chr","START_POS_REF"),all.x=TRUE)
            
            
            if (nrow(snvpredict)>1) {
              print("SNV annotation")
              snvannotation<-CDSannotation_raw(x,snvpredict)
            } else {
              print("No snv predictions. Skipping snv annotation step.")
              snvannotation=NULL
            }
            
            if (nrow(indelpredict)>1) {
              print("Indel annotation")
              indelannotation<-CDSannotation_raw(x,indelpredict)
            } else {
              print("No indel predictions. Skipping indel annotation step.")
              indelannotation=NULL
            }
            
            return(list("smurf_snv_annotation"=snvannotation, "smurf_indel_annotation"=indelannotation))
            
          }
          
          
          
          
  
      }
      else{
        print("Error: Model unrecognized.")
        write("Error: Model unrecognized.", stderr())
        
      }
    

      
    }
    else {
    print("Error: Input file check failed. One or more files may be missing or duplicated.")
    write("Error: Input file check failed. One or more files may be missing or duplicated.", stderr())
      
    }
  
  }
  else{
    print("Error: Entered directory doesn't exists (or) misspelled (or) Directory section skipped")
    write("Error: Entered directory doesn't exists (or) misspelled (or) Directory section skipped", stderr())
    
  }
  
  #rm(first_time, envir = globalenv())
  
}  
