# complete QC of Battenberg output
#
# usage:
# Rscript QC.R [RUN] [DIR_BATTENBERG] [SUBDIR_CALLSUBCLONES] [SUBDIR_POSTPROCESSING] [SUBDIR_DPCLUST] [SUBDIR_VAFPEAKS] [FILENAME_SAMPLELIST] [FILENAME_CONFIG] [FILENAME_CHRSIZES] [FILENAME_QC]

#### main script ####
if (sys.nframe() ==0L) {

        library(CleanCNA)
        args <- commandArgs(T)

        #CleanCNA:::qc(
        CleanCNA:::qc_CNAqc(
                run=as.numeric(args[1]),
                dir.battenberg=args[2],
                subdir.callsubclones=args[3],
                subdir.postprocessing=args[4],
                subdir.dpclust=args[5],
                subdir.vafpeaks=args[6],
                filename.samplelist=args[7],
                filename.config=args[8],
                filename.chr.sizes=args[9],
                filename.qc=args[10]
        )

}