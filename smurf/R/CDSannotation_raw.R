#' snvRF-Annotations
#'
#' 
#' Step 3 Annotations
#'
#'  
#' @examples
#' 
#' 
#' @export
CDSannotation_raw = function(x,predicted){
  
  print("Adding CDS and annotations")
  
  bed.to.granges <- function(fname) {
    roit <- read.table(fname,header=F)[,1:3]
    colnames(roit) <- c('chr','start','end')
    # create GRanges object
    roi.gr <- with(roit, GRanges(chr, IRanges(start, end)))
  }

  annotate.vcf=function(filename,mut.subset){

    # read in mutect2/freebayes/varscan/vardict
    vcf=readVcf(filename,"hg19",param=mut.subset)
    gr=rowRanges(vcf)

    if(length(gr)==length(mut.subset)){
      gr$ANN=info(vcf)$EFF #EFF
    } else if (length(unlist(gr$ALT))==length(gr)) {
      z=which(gr$FILTER=="PASS" & nchar(as.character(gr$REF))==1 & nchar(as.character(unlist(gr$ALT)))==1)
      gr=gr[z]
      gr$ANN=info(vcf)$EFF[z]
    } else {
      print("Error part 1") # need to fix this if there's error
      ann=matrix(,nrow = nrow(mut.subset), ncol = 12)
      ann=as.data.frame(ann)
      # colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
      #                    "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
      colnames(ann) <- c("Annotation","Impact", "Class", "Codon.Change", "AA.Change", "AA.length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                         "Feature_ID", "Exon.rank", "Allele")
      
      mut.subset <- cbind(mut.subset,ann)
      return(mut.subset)
    }

    ovl=findOverlaps(mut.subset,gr)
    mut.subset=mut.subset[queryHits(ovl)]
    mut.subset$ANN=gr[subjectHits(ovl)]$ANN
    ANN=mut.subset$ANN
    
    ##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 
    #'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | 
    #'Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )' ">
    if (length(ANN)!=0) {
      
    ANN=lapply(ANN,FUN=function(x){
      #ann=strsplit(x,split="[|]")
      ann=gsub("(.*)\\)","\\1",x) # remove last bracket
      ann=strsplit(ann,split="[(]") # extract first column
      ann <- lapply(ann, FUN=function(k) {
        c(k[1],strsplit(k[2],split="[|]")[[1]])
      }) # extract other columns
      #ann=strsplit(x,split="(")
      #ann=strsplit(ann,split=")")
      ann <- lapply(ann, FUN=function(k) {
        k=sub("^$", ".", k) # replace empty entries with [.]
        k[1:12]
      }
      )
      ann <- do.call(rbind,ann)
      ann=as.data.frame(ann)
      # colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
      #                    "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
      colnames(ann) <- c("Annotation","Impact", "Class", "Codon.Change", "AA.Change", "AA.length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                         "Feature_ID", "Exon.rank", "Allele")
      
      #colnames(ann)=colnames(mutations.orig)[16:30]
      # if (sum(ann$Feature_Type=="transcript")!=0){
        #ann=ann[which(ann$Feature_Type=="transcript"),]
        in.uniprot=which(ann$Feature_ID %in% ensembl2uni$Transcript.stable.ID)
        selected=NULL
        if (length(in.uniprot)==0){
          selected=1
        } else {
          selected=in.uniprot[1]
        }
        ann[selected,]
      # } else {
      #   print("No transcript")
      #   ann[1,]
      # }
    }
    )

    mut.subset=as.data.frame(mut.subset)
    mut.subset=mut.subset[,-which(colnames(mut.subset)=="ANN")]
    mut.subset=cbind(mut.subset,do.call(rbind,ANN))
    } else {
      ann=matrix(,nrow = length(mut.subset), ncol = 12)
      ann=as.data.frame(ann)
      # colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
      #                    "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
      colnames(ann) <- c("Annotation","Impact", "Class", "Codon.Change", "AA.Change", "AA.length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                         "Feature_ID", "Exon.rank", "Allele")
      
      mut.subset <- cbind(mut.subset,ann)
      
    }

    return(mut.subset)
  }


  # annotate ensembl ids that are present in uniprot
  smurfdir <- find.package("smurf")
  uniprotdir <- paste0(smurfdir, "/data/ensembl2uni.rds")
  #uniprotdir <- "C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/ensembl2uniprot_0917.txt"
  #ensembl2uni=read.delim(uniprotdir,header=TRUE,stringsAsFactors = FALSE)
  #ensembl2uni <- readRDS("C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/ensembl2uni.rds")
  ensembl2uni <- readRDS(uniprotdir)

  #cdsdir <- paste0(smurfdir, "/data/Ensembl75.CDS.bed")
  #cdsdir <- "C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/Ensembl75.CDS.bed"
  #cds=bed.to.granges(cdsdir)
  #roit <- read.table(cdsdir,header=F)[,1:3]
  #colnames(roit) <- c('chr','start','end')
  # create GRanges object 
  #cds <- with(roit, GRanges(chr, IRanges(start, end)))
  #saveRDS(cds,"C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/cds.rds")
  cdsdir <- paste0(smurfdir, "/data/cds.rds")
  cds <- readRDS(cdsdir)
  #cds <- readRDS("C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/cds.rds")
  
  seqlevelsStyle(cds)="NCBI"
  
  
    #tryCatch({
      

  mutations.orig=predicted
  
  #mutations.orig=predicted[[2]]
  
      #mutations.orig=snvpredict[[2]]
      #mutations.orig=indelpredict[[2]]
      
      # all predicted files should have the following fields
      mutations=with(mutations.orig,GRanges(Chr,IRanges(START_POS_REF,END_POS_REF),REF=REF,ALT=ALT,REF_MFVdVs=REF_MFVdVs,ALT_MFVdVs=ALT_MFVdVs,
                                            FILTER_Mutect2=FILTER_Mutect2,FILTER_Freebayes=FILTER_Freebayes,FILTER_Vardict=FILTER_Vardict,FILTER_Varscan=FILTER_Varscan,
                                            Sample_Name=Sample_Name,Alt_Allele_Freq=Alt_Allele_Freq,T_altDepth=T_altDepth,
                                            T_refDepth=T_refDepth,N_refDepth=N_refDepth,N_altDepth=N_altDepth))
      

      # annotate CDS region by overlapping with ensembl cds file
      ovl=findOverlaps(mutations,cds)
      mutations=mutations[unique(queryHits(ovl))]
      mutations.orig$REGION<-"Non-coding"
      mutations.orig[unique(queryHits(ovl)),"REGION"]<-"CDS"
      
      if (length(mutations)!=0){
        # mutations=mutations
        # mut.mutect2=mutations[mutations$FILTER_Mutect2==TRUE]
        # mut.vardict=mutations[mutations$FILTER_Vardict==TRUE & mutations$FILTER_Mutect2==FALSE]
        # mut.varscan=mutations[mutations$FILTER_Varscan==TRUE & mutations$FILTER_Vardict==FALSE & mutations$FILTER_Mutect2==FALSE]
        # mut.freebayes=mutations[mutations$FILTER_Freebayes==TRUE & mutations$FILTER_Varscan==FALSE & mutations$FILTER_Vardict==FALSE & mutations$FILTER_Mutect2==FALSE]
        
        mutations2=mutations
        
        mut.mutect2=mutations[mutations$FILTER_Mutect2==TRUE & nchar(gsub("/.*.","",mutations$ALT_MFVdVs))==1]
        if (length(mut.mutect2)!=0){
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.mutect2)))]
        }
        
        mut.vardict=mutations[mutations$FILTER_Vardict==TRUE & nchar(sub(".*/(.*)/.*","\\1",mutations$ALT_MFVdVs))==1]
        if (length(mut.vardict)!=0){
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.vardict)))]
        }
        
        mut.varscan=mutations[mutations$FILTER_Varscan==TRUE & nchar(gsub(".*./","",mutations$ALT_MFVdVs))==1]
        if (length(mut.varscan)!=0){
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.varscan)))]
        }
        
        mut.freebayes=mutations[mutations$FILTER_Freebayes==TRUE & nchar(sub(".*/(.*)/(.*)/.*","\\1",mutations$ALT_MFVdVs))==1]
        if (sum(length(mut.mutect2),length(mut.vardict),length(mut.varscan),length(mut.freebayes))!=length(mutations2)){
          print("Warning: missing annotations")
        }
        
        if (sum(length(mut.mutect2),length(mut.vardict),length(mut.varscan),length(mut.freebayes))==0){
          print("Warning: missing annotations, no cds transcripts found for annotations")
          
          ann=matrix(,nrow = nrow(mutations.orig), ncol = 12)
          ann=as.data.frame(ann)
          # colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
          #                    "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
          colnames(ann) <- c("Annotation","Impact", "Class", "Codon.Change", "AA.Change", "AA.length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                             "Feature_ID", "Exon.rank", "Allele")
          
          mutuations.orig <- cbind(mutations.orig[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                                     "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                                     "Sample_Name","Alt_Allele_Freq",
                                                     "T_altDepth","T_refDepth","N_refDepth","N_altDepth")],
                                   ann,
                                   mutations.orig[,c("REGION")])
          
          return(list(annotated=mutations.orig))
          
        }
          
        if (length(mut.mutect2)!=0){
          
          print("reading mutect2")
          mut.mutect2=annotate.vcf(x[[1]],mut.mutect2)

        }
        
        if (length(mut.varscan)!=0){
          
          print("reading varscan")
          mut.varscan=annotate.vcf(x[[3]],mut.varscan)

        }
        
        if (length(mut.vardict)!=0){
          
          print("reading vardict")
          mut.vardict=annotate.vcf(x[[4]],mut.vardict)

        }
        
        if (length(mut.freebayes)!=0) {
          
          print("reading freebayes")
          mut.freebayes=annotate.vcf(x[[2]],mut.freebayes)

        }
        
        print("Compiling annotations")
        start.time=Sys.time()
        
        mut=rbind(mut.mutect2,mut.vardict,mut.varscan,mut.freebayes)
        mut=mut[,-which(colnames(mut) %in% c("width","strand"))]
        mutations.orig=mutations.orig[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                         "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                         "Sample_Name","Alt_Allele_Freq",
                                         "T_altDepth","T_refDepth","N_refDepth","N_altDepth",
                                         "REGION")]
        colnames(mut)[1:17]=colnames(mutations.orig)[1:17]
        mut$Chr=as.character(mut$Chr)
        mutations.orig$Chr=as.character(mutations.orig$Chr)
        mutations.orig=merge(mut,mutations.orig,by=c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                                     "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                                     "Sample_Name","Alt_Allele_Freq",
                                                     "T_altDepth","T_refDepth","N_refDepth","N_altDepth"),all.y=TRUE)
        end.time=Sys.time() 
        
        print(end.time-start.time)
        
      } else {

        print("Warning: no cds transcripts found for annotations")
        ann=matrix(,nrow = nrow(mutations.orig), ncol = 12)
        ann=as.data.frame(ann)
        # colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
        #                    "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
        colnames(ann) <- c("Annotation","Impact", "Class", "Codon.Change", "AA.Change", "AA.length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                           "Feature_ID", "Exon.rank", "Allele")
        
        mutuations.orig <- cbind(mutations.orig[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                                   "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                                   "Sample_Name","Alt_Allele_Freq",
                                                   "T_altDepth","T_refDepth","N_refDepth","N_altDepth")],
                                 ann,
                                 mutations.orig[,c("REGION")])
        
      }

  return(list(annotated=mutations.orig))
  
  
  
}
