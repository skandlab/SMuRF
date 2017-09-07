#' SMuRF v1.1
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
#' smurf("/path/to/directory..","snv")
#' smurf("/path/to/directory..","indel")
#' smurf("/path/to/directory..","combined")
#' 
#' @export
smurf = function(directory, model, nthreads = -1){
  directory <-paste(directory,"/", sep="")
  if(dir.exists(directory)==TRUE){
    
    #setwd(directory)
    
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
      #if((packageVersion("h2o") == '3.10.3.3') == FALSE){
        print("h2o version 3.10.3.3 not found. Installing required h2o package.")
        # h2odir <- find.package("smurf1.0")
        # h2olibdir <- paste0(h2odir, "/data/h2o_3.10.3.3.tar.gz")
        # install.packages(h2olibdir)
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
    #suppressMessages(library(randomForest))
    suppressMessages(library(VariantAnnotation))
    
    suppressMessages(library(h2o))
    #suppressMessages(library(h2o))
    #suppressMessages(library(h2o, lib.loc="C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.0/data"))
    
    
    #suppressWarnings(h2o.init(nthreads = -1))
    suppressWarnings(h2o.init(nthreads = nthreads))
    

	#Retrieving files from the directory
    
    #m2 <- paste(directory,list.files(directory,pattern="CDS_mutect2.vcf", full.names=F), sep = "")
    #fb <- paste(directory,list.files(directory,pattern="CDS_freebayes.vcf", full.names=F), sep = "")
    #vs <- paste(directory,list.files(directory,pattern="CDS_varscan.vcf", full.names=F), sep = "")
    #vd <- paste(directory,list.files(directory,pattern="CDS_vardict.vcf", full.names=F), sep = "")
    #coding<-list(m2,fb,vs,vd)
    
    mutect2 <- paste(directory, Sys.glob("*mutect2*.vcf.gz"), sep = "")
    freebayes <- paste(directory, Sys.glob("*freebayes*.vcf.gz"), sep = "")
    varscan <- paste(directory, Sys.glob("*varscan*.vcf.gz"), sep = "")
    vardict <- paste(directory, Sys.glob("*vardict*.vcf.gz"), sep = "")
    
    # mutect2 <- paste(directory,list.files(directory,pattern="*mutect2*.vcf.gz$", full.names=F), sep = "")
    # freebayes <- paste(directory,list.files(directory,pattern="*freebayes*.vcf.gz$", full.names=F), sep = "")
    # varscan <- paste(directory,list.files(directory,pattern="*varscan*.vcf.gz$", full.names=F), sep = "")
    # vardict <- paste(directory,list.files(directory,pattern="*vardict*.vcf.gz$", full.names=F), sep = "")
    #x<-list(mutect2,freebayes,varscan,vardict)
    
    #x<-list(mutect2,freebayes,varscan,vardict,m2,fb,vs,vd)
    x<-list(mutect2,freebayes,varscan,vardict)
    
  
    #if(length(mutect2)==1 & length(freebayes)==1 & length(varscan)==1 & length(vardict)==1){
    if(length(grep("mutect2", mutect2))==1 & length(grep("freebayes", freebayes))==1 & length(grep("varscan", varscan))==1 & length(grep("vardict", vardict))==1){
      
      #Executing smbio if correct parameters stated
      
      if (model == "snv" || model == "indel" || model == "combined") {
          start.time <- Sys.time()
          a<-parsevcf(x)
          if (model == "snv") {
                      y<-snvRFparse(a)
                      end.time <- Sys.time()
                      time.taken <- end.time - start.time
                      #smurf_sm<<-list("smurf_snv"=y)
                      return(list("smurf_snv"=y, "time.taken"=time.taken))
          }

          else if (model == "indel") { 
                      z<-indelRFparse(a)
                      end.time <- Sys.time()
                      time.taken <- end.time - start.time
                      #smurf_sm<<-list("smurf_indel"=z)
                      return(list("smurf_indel"=z, "time.taken"=time.taken))
          }
    
          else if (model == "combined") {
                      y<-snvRFparse(a)
                      z<-indelRFparse(a)
                      end.time <- Sys.time()
                      time.taken <- end.time - start.time
					  #smurf_sm<<-list("smurf_snv"=y,"smurf_indel"=z)
                      return(list("smurf_snv"=y,"smurf_indel"=z,"time.taken"=time.taken))
 
          }
          
  
        # end.time <- Sys.time()
        # time.taken <<- end.time - start.time
      }
      else{
        print("Error: Model unrecognized.")
      }
    
      if (model == "totalfeatures") {
        print("Initializing total feature extraction.")
        start.time <- Sys.time()
        
        #to include CDS calls
        m2 <- paste(directory,list.files(directory,pattern="CDS_mutect2.vcf", full.names=F), sep = "")
        fb <- paste(directory,list.files(directory,pattern="CDS_freebayes.vcf", full.names=F), sep = "")
        vs <- paste(directory,list.files(directory,pattern="CDS_varscan.vcf", full.names=F), sep = "")
        vd <- paste(directory,list.files(directory,pattern="CDS_vardict.vcf", full.names=F), sep = "")
        
        x<-list(mutect2,freebayes,varscan,vardict,m2,fb,vs,vd)
        
        a<-parsevcfall(x)
        y<-snvRFparseall(a)
        z<-indelRFparseall(a)
        end.time <- Sys.time()
        time.taken <<- end.time - start.time
        #smbio_sm<<-list("smbio_snv"=y,"smbio_indel"=z)
        return(list("smurf_snv"=y,"smurf_indel"=z,"time.taken"=time.taken))
        
      }
      
      if (model == "totalfeaturesrejects") {
        print("Initializing total feature extraction.")
        start.time <- Sys.time()
        suppressMessages(library(data.table))
        
        #to include CDS calls
        m2 <- paste(directory,list.files(directory,pattern="CDS_mutect2.vcf", full.names=F), sep = "")
        fb <- paste(directory,list.files(directory,pattern="CDS_freebayes.vcf", full.names=F), sep = "")
        vs <- paste(directory,list.files(directory,pattern="CDS_varscan.vcf", full.names=F), sep = "")
        vd <- paste(directory,list.files(directory,pattern="CDS_vardict.vcf", full.names=F), sep = "")
        
        x<-list(mutect2,freebayes,varscan,vardict,m2,fb,vs,vd)
        
        a<-parsevcfall(x)
        y<-snvRFparseallrejects(a)
        z<-indelRFparseallrejects(a)
        end.time <- Sys.time()
        time.taken <<- end.time - start.time
        #smbio_sm<<-list("smbio_snv"=y,"smbio_indel"=z)
        return(list("smurf_snv"=y,"smurf_indel"=z,"time.taken"=time.taken))
        
      }
      

      if (model == "features") {  #new model with REGION column
        print("Initializing CDS and feature extraction.")
        start.time <- Sys.time()
        
        #to include CDS calls
        m2 <- paste(directory,list.files(directory,pattern="CDS_mutect2.vcf", full.names=F), sep = "")
        fb <- paste(directory,list.files(directory,pattern="CDS_freebayes.vcf", full.names=F), sep = "")
        vs <- paste(directory,list.files(directory,pattern="CDS_varscan.vcf", full.names=F), sep = "")
        vd <- paste(directory,list.files(directory,pattern="CDS_vardict.vcf", full.names=F), sep = "")
        
        x<-list(mutect2,freebayes,varscan,vardict,m2,fb,vs,vd)
        
        a<-parsevcffeatures(x)
        y<-snvRFparsefeatures(a)
        z<-indelRFparsefeatures(a)
        end.time <- Sys.time()
        time.taken <<- end.time - start.time
        #smbio_sm<<-list("smbio_snv"=y,"smbio_indel"=z)
        return(list("smurf_snv"=y,"smurf_indel"=z,"time.taken"=time.taken))
        
      }
      
      if (model == "annotation") {  #re-annotated ANN field, output includes VAF and callers PASSED 
        print("Initializing VAF and annotation extraction.")
        start.time <- Sys.time()
        
        a<-parsevcfANN(x)
        y<-snvRFparseANN(a)
        z<-indelRFparseANN(a)
        end.time <- Sys.time()
        time.taken <<- end.time - start.time
        #smbio_sm<<-list("smbio_snv"=y,"smbio_indel"=z)
        return(list("smurf_snv"=y,"smurf_indel"=z,"time.taken"=time.taken))
        
      }
      
      if (model == "cdsannotation") {  #re-annotated ANN field, output includes VAF and callers PASSED 
        print("Initializing CDs and ANN extraction.")
        start.time <- Sys.time()
        
        #to include CDS calls
        m2 <- paste(directory,list.files(directory,pattern="CDS_mutect2.vcf", full.names=F), sep = "")
        fb <- paste(directory,list.files(directory,pattern="CDS_freebayes.vcf", full.names=F), sep = "")
        vs <- paste(directory,list.files(directory,pattern="CDS_varscan.vcf", full.names=F), sep = "")
        vd <- paste(directory,list.files(directory,pattern="CDS_vardict.vcf", full.names=F), sep = "")
        
        x<-list(mutect2,freebayes,varscan,vardict,m2,fb,vs,vd)
        
        a<-parsevcfANNcds(x)
        y<-snvRFparseANNcds(a)
        z<-indelRFparseANNcds(a)
        end.time <- Sys.time()
        time.taken <<- end.time - start.time
        #smbio_sm<<-list("smbio_snv"=y,"smbio_indel"=z)
        return(list("smurf_snv"=y,"smurf_indel"=z,"time.taken"=time.taken))
        
      }
      
      
    }
    else {
    print("Error: Input file check failed. One or more files may be missing or duplicated.")
    }
  
  }
  else{
    print("Error: Entered directory doesn't exists (or) misspelled (or) Directory section skipped")
  }
  
  #rm(first_time, envir = globalenv())
  
}  
