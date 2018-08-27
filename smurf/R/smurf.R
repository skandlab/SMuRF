#' SMuRF v1.4
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
#' @param model Choose either "combined" or "cdsannotation"
#' The appropriate parsing and prediction model will be performed
#' to obtain a list of somatic mutation calls. 
#' Choose "combined" to generate both SNV and Indel outputs
#' Choose "cdsannotation" to run "combined" + add annotations to coding regions (from the respective coding transcripts)
#' 
#' @param nthreads Default as "-1", where all available cores will be used for RandomForest prediction. 
#' Specify any integer from 1 to x, depending on your resources available.
#' @examples
#' smurf("/path/to/directory..","combined")
#' smurf("/path/to/directory..","cdsannotation")
#' 
#' @export
smurf = function(directory, model, nthreads = -1){
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
      
      if (model == "combined" || model == "cdsannotation" || model == "featureselectionall") {
          start.time <- Sys.time()
          
          if (model == "combined") {
            parsevcf<-parsevcf_allfeaturesall(x)
            snvpredict<-snvRFpredict(parsevcf)
            indelpredict<-indelRFpredict(parsevcf)
                      end.time <- Sys.time()
                      time.taken <- end.time - start.time
                      print(time.taken)
                      
                      return(list("smurf_snv"=snvpredict, "smurf_indel"=indelpredict, "time.taken"=time.taken))
          }

          else if (model == "cdsannotation") {
            
            parsevcf<-parsevcf_allfeaturesall(x)
            
            snvpredict<-snvRFpredict(parsevcf)
            indelpredict<-indelRFpredict(parsevcf)
            
            if (length(snvpredict)>1) {
              print("SNV annotation")
              snvannotation<-CDSannotation_snv(x,snvpredict)
            } else {
              print("No snv predictions. Skipping snv annotation step.")
              snvannotation=NULL
            }
            
            if (length(indelpredict)>1) {
              print("Indel annotation")
              indelannotation<-CDSannotation_indel(x,indelpredict)
            } else {
              print("No indel predictions. Skipping indel annotation step.")
              indelannotation=NULL
            }
            
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            
            return(list("smurf_snv"=snvpredict, "smurf_indel"=indelpredict, "smurf_snv_annotation"=snvannotation, "smurf_indel_annotation"=indelannotation, "time.taken"=time.taken))
 
          }
          

          else if (model == "featureselectionall") {
            parsevcf<-parsevcf_allfeaturesall(x)
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            
            return(list("parsevcf_featureselection"=parsevcf, "time.taken"=time.taken))
            
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
