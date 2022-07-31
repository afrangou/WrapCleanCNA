# cluster using DPClust
#
# usage:
# Rscript DPClustRun.R [PLATEID_TUMOUR] [FILENAME_CELLULARITY] [FILENAME_DPCLUST_INPUT] [FILENAME_CONFIG] [DIR_OUTPUT]
#
#### code notes ####
# script required R/3.4


#### subroutines ####
log.message <- function(x, level=0) message(paste0("[", Sys.time(), "] ", paste(rep(" ", times=level * 4), collapse=""), "- ", x)) 


if (sys.nframe() == 0L) {
    args <- commandArgs(T)
    print(args)
    library(CleanCNA)
    library(lattice)#,lib.loc="/re_gecip/shared_allGeCIPs/Rpackages/3.3.0")
    library(ks)#,lib.loc="/re_gecip/shared_allGeCIPs/Rpackages/3.4.0")
    library(ggplot2)#,"/re_gecip/shared_allGeCIPs/Rpackages/3.3.0")
    library(gridExtra)#,"/re_gecip/shared_allGeCIPs/Rpackages/3.3.0")
    library(VariantAnnotation)#,"/home/afrangou/R/x86_64-pc-linux-gnu-library/3.4")
    library(DPClust)#,"/home/afrangou/R/x86_64-pc-linux-gnu-library/3.4")
    CleanCNA:::RunDPClust_Beagle_nf(
        tumourplatekey=args[1],
        normalplatekey=args[2],
        rhoandpsifilepath=args[3],
	dpinfofilepath=args[4],
        output_folder=args[5])
}


