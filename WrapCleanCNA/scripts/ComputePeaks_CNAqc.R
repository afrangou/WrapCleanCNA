# compute peaks in copy-number-state-specific VAF distributions
#
# usage:
# Rscript computePeaksVAF.R [FILENANE_SMALL_VARIANT_VCF] [FILENAME_BATTENBERG_SEGS] [FILENAME_BATTENBERG_PURITY_PLOIDY] [FILENAME_CONFIG] [FILENAME_OUTPUT]
#
#### testing ####
#filename.ssm <- "/re_gecip/cancer_renal_cell/analysisResults/1.variantCalls/A.Strelka/finalCalls/[redacted].tumo[redacted]_norm[redacted].ssm.vcf.gz"
#filename.segs <- "/re_gecip/cancer_renal_cell/analysisResults/3.chromosomeAberrations/G.CakeTin/intermediateFiles/[redacted]/tumo[redacted]_norm[redacted]/Postprocessing1/[redacted]_subclones_hg38.txt"
#filename.purity.ploidy <- "/re_gecip/cancer_renal_cell/analysisResults/3.chromosomeAberrations/G.CakeTin/intermediateFiles/[redacted]/tumo[redacted]_norm[redacted]/CallSubclones1/[redacted]_cellularity_ploidy.txt"
#filename.config <- "/home/acornish/public/code/CakeTin.nf/config/peaksVAF.config"
#filename.output <- "VAFPeaks1/[redacted].tumo[redacted]_norm[redacted].peak_data.RData"

#filename.ssm="/mnt/lustre/users/bschuster/OAC_Trial_WGS_Tissue/Results/AnnotatedVCFs/strelka2.069-001:ScrBsl:duodenum-ScrBsl.annotated.vcf.gz"
#filename.snvalleles="/mnt/lustre/users/afrangou/public/results/OACimmuno_orig/intermediateFiles/069-001/tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum/Run1/E-DPClust1/069-001:ScrBsl:tumour.output.txt"
#filename.segs="/mnt/lustre/users/afrangou/public/results/OACimmuno_orig/intermediateFiles/069-001/tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum/Run1/D2-CallSubclones1/069-001:ScrBsl:tumour_subclones.txt"
#filename.purity.ploidy="/mnt/lustre/users/afrangou/public/results/OACimmuno_orig/intermediateFiles/069-001/tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum/Run1/D2-CallSubclones1/069-001:ScrBsl:tumour_cellularity_ploidy.txt"
#filename.output="/mnt/lustre/users/afrangou/public/results/OACimmuno_orig/intermediateFiles/069-001/tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum/Run1/F-VAFPeaks1/069-001.tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum.peak_data.RData"
#output_dir="/mnt/lustre/users/afrangou/public/results/OACimmuno_orig/intermediateFiles/069-001/tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum/Run1/F-VAFPeaks1/"
#filename.output.vafpeaks="/mnt/lustre/users/afrangou/public/results/OACimmuno_orig/intermediateFiles/069-001/tumo069-001:ScrBsl:tumour_norm069-001:ScrBsl:duodenum/Run1/F-VAFPeaks1/069-001:ScrBsl:tumour_vafpeaks.pdf"
#
#069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum
#
#filename.ssm="/mnt/lustre/users/bschuster/OAC_Trial_WGS_Tissue/Results/AnnotatedVCFs/strelka2.069-013:ScrBsl:duodenum-ScrBsl.annotated.vcf.gz"
#filename.snvalleles="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run1/E-DPClust1/069-013:ScrBsl:tumour.output.txt"
#filename.segs="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run1/D2-CallSubclones1/069-013:ScrBsl:tumour_subclones.txt"
#filename.purity.ploidy="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run1/D2-CallSubclones1/069-013:ScrBsl:tumour_cellularity_ploidy.txt"
#filename.output="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run1/F-VAFPeaks1/069-001.tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum.peak_data.RData"
#output_dir="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run1/F-VAFPeaks1/"
#filename.output.vafpeaks="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run1/F-VAFPeaks1/069-013:ScrBsl:tumour_vafpeaks.pdf"
#
#filename.ssm="/mnt/lustre/users/bschuster/OAC_Trial_WGS_Tissue/Results/AnnotatedVCFs/strelka2.069-013:ScrBsl:duodenum-ScrBsl.annotated.vcf.gz"
#filename.snvalleles="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/H-DPClust2/069-013:ScrBsl:tumour.output.txt"
#filename.segs="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/G2-CallSubclones2/069-013:ScrBsl:tumour_subclones.txt"
#filename.purity.ploidy="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/G2-CallSubclones2/069-013:ScrBsl:tumour_cellularity_ploidy.txt"
#filename.output="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/I-VAFPeaks2/069-001.tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum.peak_data.RData"
#output_dir="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/I-VAFPeaks2/"
#filename.output.vafpeaks="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/I-VAFPeaks2/069-013:ScrBsl:tumour_vafpeaks.pdf"


#filename.snvalleles="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/H-DPClust2/069-013:ScrBsl:tumour.output.txt"
#filename.segs="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/G2-CallSubclones2/069-013:ScrBsl:tumour_subclones.txt"
#filename.purity.ploidy="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/G2-CallSubclones2/069-013:ScrBsl:tumour_cellularity_ploidy.txt"
#filename.output="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/I-VAFPeaks2/069-001.tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum.peak_data.RData"
#output_dir="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/I-VAFPeaks2/"
#filename.output.vafpeaks="/mnt/lustre/users/afrangou/public/results/OACimmuno_extras3/intermediateFiles/069-013/tumo069-013:ScrBsl:tumour_norm069-013:ScrBsl:duodenum/Run2/I-VAFPeaks2/069-013:ScrBsl:tumour_vafpeaks.pdf"


#### subroutines ####
log.message <- function(x, level=0) message(paste0("[", Sys.time(), "] ", paste(rep(" ", times=level * 4), collapse=""), "- ", x))

if (sys.nframe() == 0L) {

    args <- commandArgs(T)
    print(args)

    library(CleanCNA)
    library(peakPick)
    library(ggrepel)
    library(pio)
    library(tidyverse)
    library(CNAqc)
    library(crayon)

    filename.ssm=args[1]
    filename.snvalleles=args[2]
    filename.segs=args[3]
    filename.purity.ploidy=args[4]
    #filename.config=args[5]
    filename.output=args[6]
    output_dir=args[7]
    filename.output.vafpeaks=args[8]

    # load snv counts
    snvcounts=read.table(filename.ssm)

    # load allele file used for dpclust
    snvalleles=read.table(filename.snvalleles,hea=T,stringsAsFactors=F)

    # match positions to those used in dpclust
    snvcounts_unique=paste(snvcounts[,1],snvcounts[,2],sep="_")
    snvalleles_unique=paste(snvalleles[,1],snvalleles[,2],sep="_")

    touse=match(snvalleles_unique,snvcounts_unique)
    snvalleles=snvalleles[which(!is.na(touse)),]
    if (length(which(!is.na(touse)))>0) {touse=touse[which(!is.na(touse))]}
    snvcounts=snvcounts[touse,]
    snvcounts[,1]=paste0('chr',snvcounts[,1])

    # reshape data for input
    snvs=cbind(snvcounts[,c(1,2,2,4,5)],snvalleles[,3:7])
    colnames(snvs)[1:5]=c("chr","from","to","ref","alt")
    
    # recode alleles
    snvs$refnum[snvs$ref=="A"]=1
    snvs$refnum[snvs$ref=="C"]=2
    snvs$refnum[snvs$ref=="G"]=3
    snvs$refnum[snvs$ref=="T"]=4

    snvs$altnum[snvs$alt=="A"]=1
    snvs$altnum[snvs$alt=="C"]=2
    snvs$altnum[snvs$alt=="G"]=3
    snvs$altnum[snvs$alt=="T"]=4

    # col for ref and alt
    snvs$refnum=as.integer(snvs$refnum)+5
    snvs$altnum=as.integer(snvs$altnum)+5

    # counts of ref and alt
    snvs$refcounts=snvs[cbind(seq_along(snvs$refnum),snvs$refnum)]
    snvs$altcounts=snvs[cbind(seq_along(snvs$altnum),snvs$altnum)]
 
    # DP=total depth, NV=count of alt
    snvs$DP=as.numeric(snvs[,10])
    snvs$NV=as.numeric(snvs$altcounts)
    snvs$VAF=snvs$NV/snvs$DP

    # required format
    snvs=snvs[,c(1:5,15:17)]
    for (i in c(1,4,5)) {snvs[,i]=as.character(snvs[,i])}
    for (i in c(2,3,6,7)) {snvs[,i]=as.numeric(snvs[,i])}

    # remove indels
    if (length(which(nchar(snvs$ref)>1))>0) {snvs=snvs[-which(nchar(snvs$ref)>1),]}
    if (length(which(nchar(snvs$alt)>1))>0) {snvs=snvs[-which(nchar(snvs$alt)>1),]}

    # get segs info
    cna=read.table(filename.segs,hea=T,stringsAsFactors=F)[,1:13]
    cna=as.data.frame(cbind(paste0("chr",cna$chr),cna[,c(2,3,8,9,10)]))
    colnames(cna)=c("chr","from","to","Major","minor","CCF")
    for (i in c(2:5)) {cna[,i]=as.integer(as.character(cna[,i]))}
    cna[,6]=as.numeric(as.character(cna[,6]))
    cna[,1]=as.character(cna[,1])

    # get purity
    purity=read.table(filename.purity.ploidy,sep="\t",hea=T,stringsAsFactors=F)[1,1]
    # define ref
    ref="GRCh37"

    # create CNAqc data frame
    data=init(snvs,cna,purity,ref)

    peaks=analyze_peaks(data,
			 karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'),
                         min_karyotype_size = 0.05,
                         min_absolute_karyotype_mutations = 50,
                         p_binsize_peaks = 0.005,
                         purity_error = 0.025,
			 n_bootstrap = 1, 
    			 kernel_adjust = 1, 
			 matching_strategy = "closest")

    # Fix tlesluyes
    if (is.null(peaks$peaks_analysis)) peaks$peaks_analysis$QC='FLAG'
    
    save(peaks,file=filename.output)
    if (peaks$peaks_analysis$QC=="FLAG") {write.csv("FLAG",paste0(output_dir,"/FLAG"),quote=F,row.names=F)}
    if (peaks$peaks_analysis$QC=="PASS") {write.csv("PASSPEAKS",paste0(output_dir,"/PASSPEAKS"),quote=F,row.names=F)}
    if (peaks$peaks_analysis$QC=="FAIL") {write.csv("FAILPEAKS",paste0(output_dir,"/FAILPEAKS"),quote=F,row.names=F)}

    png(filename.output.vafpeaks,width=25,height=15,units='cm',res=300,pointsize=6)
    plot(plot_peaks_analysis(peaks))
    dev.off()
}
