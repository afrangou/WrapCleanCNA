# fit copy number, call subclones and switch from hg37 to hg38
#
# usage:
# Rscript BattenbergFitCopyNumber.R [PLATEKEY_TUMOUR] [SEX] [PREFIX_OUTPUT_BAF] [PREFIX_OUTPUT_GC] [PREFIX_OUTPUT_SEGMENT] [PURITY_PRESET] [PLOIDY_PRESET] [PLOIDY_MIN] [PLOIDY_MAX] [FILENAME_BATTENBERG_CONFIG] [FILENAME_IMPUTE_INFO] [PREFIX_OUTPUT_FIT] [PREFIX_OUTPUT_SUBCLONES]
#
#### functions ####
log.message <- function(x, level=0) message(paste0("[", Sys.time(), "] ", paste(rep(" ", times=level * 4), collapse=""), "- ", x))


#### main script ####
if (sys.nframe() == 0L) {
    # handle parameters
    args <- commandArgs(T)
    plateid.tumo <- args[1]
    sex <- args[2]
    prefix.output.baf <- args[3]
    prefix.output.gc <- args[4]
    prefix.output.segment <- args[5]
    preset.purity <- as.numeric(args[6]) #ifelse(args[6] == "NA", NA, as.numeric(args[6]))
    preset.ploidy <- as.numeric(args[7]) #ifelse(args[7] == "NA", NA, as.numeric(args[7]))
    min.ploidy <- as.numeric(args[8])
    max.ploidy <- as.numeric(args[9])
    filename.config <- args[10]
    filename.impute.info <- args[11]
    prefix.output.fit <- args[12]
    prefix.output.subclones <- args[13]
    prior_breakpoints_file = args[14]

    # required libraries and parameters
    library(Battenberg)
    source(filename.config)

    print(preset.purity)
    print(preset.ploidy)

    if (!file.exists(prior_breakpoints_file)) {
        print(paste0('prior_breakpoints_file (',prior_breakpoints_file,') does not exist'))
        prior_breakpoints_file=NULL
    } else {
        print(paste0('prior_breakpoints_file (',prior_breakpoints_file,') exists'))
    }

    # fit the copy number
    log.message("fitting copy number")
    file.copy(paste0(prefix.output.segment, ".BAFsegmented.txt"),paste0(prefix.output.subclones, ".BAFsegmented.txt"),overwrite=T)
    
    fit.copy.number(
	samplename=plateid.tumo,
	outputfile.prefix=paste0(prefix.output.fit, "_"),
	inputfile.baf.segmented=paste0(prefix.output.subclones, ".BAFsegmented.txt"),
	inputfile.baf=paste0(prefix.output.baf, "_mutantBAF.tab"),
	inputfile.logr=paste0(prefix.output.gc, "_mutantLogR_gcCorrected.tab"),
	use_preset_rho_psi=all(!is.na(c(preset.purity, preset.ploidy))),
	preset_rho=preset.purity,
	preset_psi=preset.ploidy,
	min.ploidy=min.ploidy,
	max.ploidy=max.ploidy,
	min.rho=min.purity,
	max.rho=max.purity,
	dist_choice=clonality.dist.metric,
	ascat_dist_choice=ascat.dist.metric,
	min.goodness=min.goodness.of.fit,
	uninformative_BAF_threshold=balanced.threshold,
	gamma_param=platform.gamma,
	read_depth=read.depth
    )

    # call subclones
    log.message("calling subclones")
    callSubclones(
	sample.name=plateid.tumo,
	baf.segmented.file=paste0(prefix.output.subclones, ".BAFsegmented.txt"),
	logr.file=paste0(prefix.output.gc, "_mutantLogR_gcCorrected.tab"),
	rho.psi.file=paste0(prefix.output.fit, "_rho_and_psi.txt"),
	output.file=paste0(prefix.output.subclones, "_subclones.txt"),
	output.figures.prefix=paste0(prefix.output.subclones, "_subclones_chr"),
	output.gw.figures.prefix=paste0(prefix.output.subclones, "_BattenbergProfile"),
	chr_names=get.chrom.names(filename.impute.info, sex == "MALE"),
	masking_output_file=paste0(prefix.output.subclones, "_segment_masking_details.txt"),
	prior_breakpoints_file=prior_breakpoints_file,
	max_allowed_state=max.cn.state,
	gamma=platform.gamma,
	segmentation.gamma=seg.gamma,
	siglevel=siglevel.subclones,
	maxdist=maxdist.subclones,
	noperms=noperms.subclones,
	seed=seed)
    system(paste0('ln -s ',prefix.output.subclones,'_purity_ploidy.txt ',prefix.output.subclones,'_cellularity_ploidy.txt'))
}
