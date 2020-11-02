#' snvRF-Annotations
#' Step 3 Annotations
#'
#'@param x List object containing the four vcf.gz files from 
#' callers Strelka2, MuTect2, Freebayes, VarDict and VarScan.
#' @param tbi TabixFile intance of VCF file
#' @param predicted List object containing 1.stats, 2.predicted-snv and 3.parse-snv matrices
#' @param build Your current genomic build
#' @param change.build Choose TRUE to convert from 'hg19.to.hg38' or 'hg38.to.hg19' annotation output
#'  
#'@examples
#' 
#' 
#'@export
CDSannotation_snv = function(x, tbi, predicted, build, change.build){
  
  print("Adding CDS and annotations")
  
  bed.to.granges <- function(fname) {
    roit <- read.table(fname,header=F)[,1:3]
    colnames(roit) <- c('chr','start','end')
    # create GRanges object
    roi.gr <- with(roit, GRanges(chr, IRanges(start, end)))
  }

  annotate.vcf=function(filename,build,mut.subset){

    # read in mutect2/freebayes/varscan/vardict/strelka2
    # vcf=readVcf(filename,"hg19",param=mut.subset)
    vcf=readVcf(filename,genome = build,param=mut.subset)
    gr=rowRanges(vcf)

    if(length(gr)==length(mut.subset)){
      gr$ANN=info(vcf)$ANN
    } else if (length(unlist(gr$ALT))==length(gr)) {
      z=which(gr$FILTER=="PASS" & nchar(as.character(gr$REF))==1 & nchar(as.character(unlist(gr$ALT)))==1)
      gr=gr[z]
      gr$ANN=info(vcf)$ANN[z]
    } else {
      print("Error part 1") # need to fix this if there's error
      ann=matrix(,nrow = length(mut.subset), ncol = 15)
      ann=as.data.frame(ann)
      colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
                         "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
      mut.subset <- cbind(mut.subset,ann)
      return(mut.subset)
    }

    ovl=findOverlaps(mut.subset,gr)
    mut.subset=mut.subset[queryHits(ovl)]
    mut.subset$ANN=gr[subjectHits(ovl)]$ANN
    ANN=mut.subset$ANN

    ANN=lapply(ANN,FUN=function(x){
      ann=strsplit(x,split="[|]")
      ann <- lapply(ann, FUN=function(k) {
        k=sub("^$", ".", k)
        k[1:15]
      }
      )
      ann <- do.call(rbind,ann)
      ann=as.data.frame(ann)
      colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
                         "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
      #colnames(ann)=colnames(mutations.orig)[16:30]
      if (sum(ann$Feature_Type=="transcript")!=0){
        ann=ann[which(ann$Feature_Type=="transcript"),]
        in.uniprot=which(ann$Feature_ID %in% ensembl2uni$Transcript.stable.ID)
        selected=NULL
        if (length(in.uniprot)==0){
          selected=1
        } else {
          selected=in.uniprot[1]
        }
        ann[selected,]
      } else {
        print("No transcript")
        ann[1,]
      }
    }
    )

    mut.subset=as.data.frame(mut.subset)
    mut.subset=mut.subset[,-which(colnames(mut.subset)=="ANN")]
    mut.subset=cbind(mut.subset,do.call(rbind,ANN))

    return(mut.subset)
  }


    #tryCatch({
      
      mutations.orig=predicted[[2]]
      #mutations.orig=snvpredict[[2]]
      #mutations.orig=indelpredict[[2]]
      #mutations.orig=snv.roi[[2]]
      #mutations.orig=test.smurf$smurf_snv$predicted_snv
      
      # annotate ensembl ids that are present in uniprot
      smurfdir <- find.package("smurf")
      uniprotdir <- paste0(smurfdir, "/data/ensembl2uni.rds")
      ensembl2uni <- readRDS(uniprotdir)
      
      #hg19 or hg38 for cds
      if(build=='hg19') {
        
        cdsdir <- paste0(smurfdir, "/data/cds.rds") #HG19
        cds <- readRDS(cdsdir)
        seqlevelsStyle(cds)="NCBI"
        
      } else if(build=='hg38') {
        
        cdsdir <- paste0(smurfdir, "/data/cds-hg38.rds")
        cds <- readRDS(cdsdir)
        seqlevelsStyle(cds)="UCSC"
        mutations.orig$Chr = paste0('chr',mutations.orig$Chr)
        
      }
      
      # all predicted files should have the following fields
      mutations=with(mutations.orig,GRanges(Chr,IRanges(START_POS_REF,END_POS_REF),REF=REF,ALT=ALT,REF_MFVdVs=REF_MFVdVs,ALT_MFVdVs=ALT_MFVdVs,
                                            FILTER_Mutect2=FILTER_Mutect2,FILTER_Freebayes=FILTER_Freebayes,FILTER_Vardict=FILTER_Vardict,FILTER_Varscan=FILTER_Varscan,FILTER_Strelka2=FILTER_Strelka2,
                                            Sample_Name=Sample_Name,Alt_Allele_Freq=Alt_Allele_Freq,
                                            N_refDepth=N_refDepth,N_altDepth=N_altDepth,
                                            T_refDepth=T_refDepth,T_altDepth=T_altDepth))
      
      # annotate CDS region by overlapping with ensembl cds file
      ovl=suppressWarnings(findOverlaps(mutations,cds)) #Warning message: In .Seqinfo.mergexy(x, y) 
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
        
        suppressWarnings(suppressMessages(library(dplyr)))
        
        REF = strsplit(mutations$REF_MFVdVs, split='/')
        names(REF) = c(1:length(REF))
        REF = t(bind_rows(REF))
        REF[REF=='NA'] = ''
        REF=as.data.frame(REF)
        
        ALT = strsplit(mutations$ALT_MFVdVs, split='/')
        names(ALT) = c(1:length(ALT))
        ALT = t(bind_rows(ALT))
        ALT[ALT=='NA'] = ''
        ALT=as.data.frame(ALT)
        
        mut.mutect2=mutations[mutations$FILTER_Mutect2==TRUE & nchar(as.vector(as.matrix(REF[,1])))==1 & nchar(as.vector(as.matrix(ALT[,1])))==1]
        # mut.mutect2=mutations[mutations$FILTER_Mutect2==TRUE & nchar(gsub("/.*.","",mutations$ALT_MFVdVs))==1 & nchar(gsub("/.*.","",mutations$REF_MFVdVs))==1]
        # mut.mutect2=mutations[nchar(gsub("/.*.","",mutations$ALT_MFVdVs))==1 & nchar(gsub("/.*.","",mutations$REF_MFVdVs))==1]
        if (length(mut.mutect2)!=0){
          #remove using index from REF and ALT matrix 
          # REF[unique(queryHits(findOverlaps(mutations,mut.mutect2))),] = ''
          # ALT[unique(queryHits(findOverlaps(mutations,mut.mutect2))),] = ''
          REF=REF[-unique(queryHits(findOverlaps(mutations,mut.mutect2))),]
          ALT=ALT[-unique(queryHits(findOverlaps(mutations,mut.mutect2))),]
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.mutect2)))]
        }
        
        mut.strelka2=mutations[mutations$FILTER_Strelka2==TRUE & nchar(as.vector(as.matrix(REF[,5])))==1 & nchar(as.vector(as.matrix(ALT[,5])))==1]
        if (length(mut.strelka2)!=0){
          #remove using index from REF and ALT matrix 
          REF=REF[-unique(queryHits(findOverlaps(mutations,mut.strelka2))),]
          ALT=ALT[-unique(queryHits(findOverlaps(mutations,mut.strelka2))),]
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.strelka2)))]
        }
        
        mut.vardict=mutations[mutations$FILTER_Vardict==TRUE & nchar(as.vector(as.matrix(REF[,3])))==1 & nchar(as.vector(as.matrix(ALT[,3])))==1]
        # mut.vardict=mutations[mutations$FILTER_Vardict==TRUE & nchar(sub(".*/(.*)/.*","\\1",mutations$ALT_MFVdVs))==1 & nchar(sub(".*/(.*)/.*","\\1",mutations$REF_MFVdVs))==1]
        # mut.vardict=mutations[nchar(sub(".*/(.*)/.*","\\1",mutations$ALT_MFVdVs))==1 & nchar(sub(".*/(.*)/.*","\\1",mutations$REF_MFVdVs))==1]
        if (length(mut.vardict)!=0){
          #remove using index from REF and ALT matrix 
          REF=REF[-unique(queryHits(findOverlaps(mutations,mut.vardict))),]
          ALT=ALT[-unique(queryHits(findOverlaps(mutations,mut.vardict))),]
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.vardict)))]
        }
        
        mut.varscan=mutations[mutations$FILTER_Varscan==TRUE & nchar(as.vector(as.matrix(REF[,4])))==1 & nchar(as.vector(as.matrix(ALT[,4])))==1]
        # mut.varscan=mutations[mutations$FILTER_Varscan==TRUE & nchar(gsub(".*./","",mutations$ALT_MFVdVs))==1 & nchar(gsub(".*./","",mutations$REF_MFVdVs))==1]
        # mut.varscan=mutations[nchar(gsub(".*./","",mutations$ALT_MFVdVs))==1 & nchar(gsub(".*./","",mutations$REF_MFVdVs))==1]
        if (length(mut.varscan)!=0){
          #remove using index from REF and ALT matrix 
          REF=REF[-unique(queryHits(findOverlaps(mutations,mut.varscan))),]
          ALT=ALT[-unique(queryHits(findOverlaps(mutations,mut.varscan))),]
          mutations=mutations[-unique(queryHits(findOverlaps(mutations,mut.varscan)))]
        }
        
        mut.freebayes=mutations[mutations$FILTER_Freebayes==TRUE & nchar(as.vector(as.matrix(REF[,2])))==1 & nchar(as.vector(as.matrix(ALT[,2])))==1]
        # mut.freebayes=mutations[mutations$FILTER_Freebayes==TRUE & nchar(sub(".*/(.*)/(.*)/.*","\\1",mutations$ALT_MFVdVs))==1 & nchar(sub(".*/(.*)/(.*)/.*","\\1",mutations$REF_MFVdVs))==1]
        # mut.freebayes=mutations[nchar(sub(".*/(.*)/(.*)/.*","\\1",mutations$ALT_MFVdVs))==1 & nchar(sub(".*/(.*)/(.*)/.*","\\1",mutations$REF_MFVdVs))==1]
        
        if (sum(length(mut.mutect2),length(mut.strelka2),length(mut.vardict),length(mut.varscan),length(mut.freebayes))!=length(mutations2)){
          print("Warning: missing annotations")
        }
          
        if (length(mut.mutect2)!=0){
          print("reading mutect2")
          # mut.mutect2=suppressWarnings(annotate.vcf(x[[1]],build,mut.mutect2))
          mut.mutect2=suppressWarnings(annotate.vcf(tbi[[1]],build,mut.mutect2))
          # mut.mutect2=suppressWarnings(annotate.vcf(tbi[[1]],build=seqinfo(scanVcfHeader(x[[1]])),mut.mutect2))
          
        }
        
        if (length(mut.strelka2)!=0){
          print("reading strelka2")
          mut.strelka2=suppressWarnings(annotate.vcf(tbi[[5]],build,mut.strelka2))
        }
        
        if (length(mut.varscan)!=0){
          print("reading varscan")
          mut.varscan=suppressWarnings(annotate.vcf(tbi[[3]],build,mut.varscan))
        }
        
        if (length(mut.vardict)!=0){
          print("reading vardict")
          mut.vardict=suppressWarnings(annotate.vcf(tbi[[4]],build,mut.vardict))
        }
        
        if (length(mut.freebayes)!=0) {
          print("reading freebayes")
          mut.freebayes=suppressWarnings(annotate.vcf(tbi[[2]],build,mut.freebayes))
        }
        
        print("Compiling annotations")
        start.time=Sys.time()
        
        mut=rbind(mut.mutect2,mut.strelka2,mut.vardict,mut.varscan,mut.freebayes)
        mut=mut[,-which(colnames(mut) %in% c("width","strand"))]
        mutations.orig=mutations.orig[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                         "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                                         "Sample_Name","Alt_Allele_Freq",
                                         "N_refDepth","N_altDepth","T_refDepth","T_altDepth",
                                         "REGION","SMuRF_score")]
        colnames(mut)[1:18]=colnames(mutations.orig)[1:18]
        mut$Chr=as.character(mut$Chr)
        mutations.orig$Chr=as.character(mutations.orig$Chr)
        mutations.orig=merge(mut,mutations.orig,by=c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                                     "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                                                     "Sample_Name","Alt_Allele_Freq",
                                                     "N_refDepth","N_altDepth","T_refDepth","T_altDepth"),all.y=TRUE)
        
        mutations.orig$Chr = gsub('chr','',mutations.orig$Chr)
        
        if (change.build==T & build=='hg19') {
          print('Changing annotations from hg19 to hg38')
          hg.gr = with(mutations.orig,GRanges(Chr,IRanges(START_POS_REF,END_POS_REF)))
          seqlevelsStyle(hg.gr)="UCSC"
          suppressWarnings(suppressMessages(library(rtracklayer)))
          chainObject <- import.chain(paste0(smurfdir,"/data/hg19ToHg38.over.chain"))
          hg38 <- as.data.frame(liftOver(hg.gr, chainObject))
          mutations.orig[,c('Chr','START_POS_REF','END_POS_REF')] = hg38[,c('seqnames','start','end')]
          mutations.orig$Chr = gsub('chr','',mutations.orig$Chr)
        } else if (change.build==T & build=='hg38') {
          print('Changing annotations from hg38 to hg19')
          hg.gr = with(mutations.orig,GRanges(Chr,IRanges(START_POS_REF,END_POS_REF)))
          seqlevelsStyle(hg.gr)="UCSC"
          suppressWarnings(suppressMessages(library(rtracklayer)))
          chainObject <- import.chain(paste0(smurfdir,"/data/hg38ToHg19.over.chain"))
          hg19 <- as.data.frame(liftOver(hg.gr, chainObject))
          mutations.orig[,c('Chr','START_POS_REF','END_POS_REF')] = hg19[,c('seqnames','start','end')]
          mutations.orig$Chr = gsub('chr','',mutations.orig$Chr)
        }
        
        
        end.time=Sys.time() 
        
        print(end.time-start.time)
        
      } else {

        print("Warning: no cds transcripts found for annotations")
        ann=matrix(,nrow = nrow(mutations.orig), ncol = 15)
        ann=as.data.frame(ann)
        colnames(ann) <- c("Allele","Annotation", "Impact", "Gene_name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
                           "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance")
        mutations.orig <- cbind(mutations.orig[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                                  "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                                                   "Sample_Name","Alt_Allele_Freq",
                                                   "N_refDepth","N_altDepth","T_refDepth","T_altDepth")],
                                 ann,
                                 mutations.orig[,c("REGION", "SMuRF_score")])
        
        mutations.orig$Chr = gsub('chr','',mutations.orig$Chr)
        
        if (change.build==T & build=='hg19') {
          print('Changing annotations from hg19 to hg38')
          hg.gr = with(mutations.orig,GRanges(Chr,IRanges(START_POS_REF,END_POS_REF)))
          seqlevelsStyle(hg.gr)="UCSC"
          suppressWarnings(suppressMessages(library(rtracklayer)))
          chainObject <- import.chain(paste0(smurfdir,"/data/hg19ToHg38.over.chain"))
          hg38 <- as.data.frame(liftOver(hg.gr, chainObject))
          mutations.orig[,c('Chr','START_POS_REF','END_POS_REF')] = hg38[,c('seqnames','start','end')]
          mutations.orig$Chr = gsub('chr','',mutations.orig$Chr)
        } else if (change.build==T & build=='hg38') {
          print('Changing annotations from hg38 to hg19')
          hg.gr = with(mutations.orig,GRanges(Chr,IRanges(START_POS_REF,END_POS_REF)))
          seqlevelsStyle(hg.gr)="UCSC"
          suppressWarnings(suppressMessages(library(rtracklayer)))
          chainObject <- import.chain(paste0(smurfdir,"/data/hg38ToHg19.over.chain"))
          hg19 <- as.data.frame(liftOver(hg.gr, chainObject))
          mutations.orig[,c('Chr','START_POS_REF','END_POS_REF')] = hg19[,c('seqnames','start','end')]
          mutations.orig$Chr = gsub('chr','',mutations.orig$Chr)
        }
        
        
      }

  return(list(annotated=mutations.orig))
  
  
  
}
