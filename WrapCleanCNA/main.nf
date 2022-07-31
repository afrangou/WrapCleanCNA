#!/usr/bin/env nextflow
//**********************************************************
// Nextflow pipeline for running Battenberg and DPClust and conducting automated QC
//**********************************************************

//**********************************************************
// Help message
//**********************************************************

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run ./CNVCalling.nf --samplelist [samplelist.tsv] --dir_output [dir_output]

    Mandatory arguments:
      --samplelist      Tab-separated file containing information on each tumour-normal pair. Must containing the following columns: participant_id, sex [MALE/FEMALE], tumour_sample_platekey, germline_sample_platekey, somatic_small_variants_qced_vcf, tumour_bam, germline_bam. Platekeys are used to select columns from VCFs.
      --dir_output  	Directory to output results.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit Process
}



//**********************************************************
// Parameters
//**********************************************************

file_samplelist = file(params.samplelist)
dir_output_base = params.dir_output



//**********************************************************
// Get sample info from tab separated file 
//**********************************************************

// read in samplelist
Channel
    .from(file_samplelist)
    .splitCsv(header: true, sep: '\t')
    .map { row -> [
        dir_output_base + "/intermediateFiles/" + row.participant_id + "/tumo" + row.tumour_sample_platekey + "_norm" + row.germline_sample_platekey,
	row.participant_id,
	row.sex,
	row.tumour_sample_platekey,
	row.germline_sample_platekey,
	row.somatic_small_variants_qced_vcf,
	row.tumour_bam,
	row.germline_bam
    ] }
    .set { samplelist1_ch }


///////////////////////////////////////////////////////////////////////
// first DPClust run
///////////////////////////////////////////////////////////////////////
process initialrun_stage2_dpclust {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    //beforeScript { task.attempt<=1 ? 'sleep 0' : 'sleep 600' }
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv, tumour_bam, germline_bam from samplelist1_ch

    output:
    val dir_output into dpclust1_ch
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv into fit1_ch2
    file "*" optional true into dpclust1_files_ch

    shell:
    '''
    module load R/3.6.2-foss-2019b
    if [ ! -s "!{dir_output}/Run1/E-DPClust1/!{plateid_tumo}_optimaInfo.txt" ]; then
         mkdir -p !{dir_output}/Run1/E-DPClust1
         Rscript !{baseDir}/scripts/DPClustPrep.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{sex}\
            !{filename_snv}\
            !{dir_output}/Run1/D1-FitCopyNumber1/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run1/D2-CallSubclones1/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run1/E-DPClust1/!{plateid_tumo}.loci.txt\
            !{dir_output}/Run1/E-DPClust1/!{plateid_tumo}.output.txt\
            !{dir_output}/Run1/E-DPClust1/!{plateid_tumo}.DPinput.txt

        Rscript !{baseDir}/scripts/DPClustRun.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{dir_output}/Run1/D1-FitCopyNumber1/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run1/E-DPClust1/!{plateid_tumo}.DPinput.txt\
            !{dir_output}/Run1/E-DPClust1

    fi
    '''

}



///////////////////////////////////////////////////////////////////////
// assess peaks in copy-number-state-specific VAF distributions
///////////////////////////////////////////////////////////////////////
process initialrun_stage3_vaf_peaks {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv from fit1_ch2 

    output:
    val dir_output into peaksvaf1_ch
    file "*" optional true into peaksvaf1_files_ch

    shell:
    '''
    module load R/3.6.2-foss-2019b
    basename_output=!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData
    if [ ! -s "!{dir_output}/Run1/F-VAFPeaks1/${basename_output}" ]; then
        mkdir -p !{dir_output}/Run1/F-VAFPeaks1
            Rscript !{baseDir}/scripts/ComputePeaks_CNAqc.R\
                !{filename_snv}\
                !{dir_output}/Run1/E-DPClust1/!{plateid_tumo}.output.txt\
                !{dir_output}/Run1/D2-CallSubclones1/!{plateid_tumo}_subclones.txt\
                !{dir_output}/Run1/D2-CallSubclones1/!{plateid_tumo}_cellularity_ploidy.txt\
                !{baseDir}/config/peaksVAF.config\
                !{dir_output}/Run1/F-VAFPeaks1/!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData\
                !{dir_output}/Run1/F-VAFPeaks1\
                !{dir_output}/Run1/F-VAFPeaks1/!{plateid_tumo}_vafpeaks.png
    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// first QC run
///////////////////////////////////////////////////////////////////////
process initialrun_stage4_qc {
    errorStrategy 'terminate'
    publishDir "${dir_output_base}", mode: "copy", overwrite: true

    input:
    val x from dpclust1_ch.collect()
    val x from peaksvaf1_ch.collect()

    output:
    file "QualityControl_noBB.tsv" into qc1_table_file_ch

    shell:
    '''
    module load R/3.6.2-foss-2019b

    Rscript !{baseDir}/scripts/QualityControl.R\
        1\
        !{dir_output_base}\
        /Run1/D2-CallSubclones1\
        /Run1/D2-CallSubclones1\
        /Run1/E-DPClust1\
        /Run1/F-VAFPeaks1\
        !{file_samplelist}\
        !{baseDir}/config/qc.config\
        !{baseDir}/data/hg19.chrom.sizes\
        QualityControl_noBB.tsv
    '''
}



///////////////////////////////////////////////////////////////////////
// FIRST RERUN
///////////////////////////////////////////////////////////////////////
// refit copy number of samples failing first round QC
qc1_table_file_ch
    .splitCsv(header: true, sep: '\t')
    .map { row -> [
        dir_output_base + "/intermediateFiles/" + row.participant_id + "/tumo" + row.tumour_sample_platekey + "_norm" + row.germline_sample_platekey,
        row.participant_id,
        row.sex,
        row.tumour_sample_platekey,
        row.germline_sample_platekey,
        row.run1_reestimated_purity,
        row.run1_reestimated_ploidy,
        row.somatic_small_variants_qced_vcf,
        row.run1_filter_overallfilter
    ] }
    .filter { it[8] == "FAIL" || it[8] == "FLAG" }
    .set { samplelist2_ch }

process rerun1_stage1_fit_copy_number {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, purity, ploidy, filename_snv, filter_overall from samplelist2_ch

    output:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv into fit2_ch1
    file "*" optional true into fit2_files_ch

    shell:
    '''
    module load R/3.6.2-foss-2019b

    if [ ! -s "!{dir_output}/Run2/G2-CallSubclones2/!{plateid_tumo}_cellularity_ploidy.txt" ]; then
        mkdir -p !{dir_output}/Run2 !{dir_output}/Run2/G1-FitCopyNumber2 !{dir_output}/Run2/G2-CallSubclones2
        Rscript !{baseDir}/scripts/BattenbergFitCopyNumber.R\
             !{plateid_tumo}\
             !{sex}\
             !{dir_output}/Run1/B1-RunBAFLogR/!{plateid_tumo}\
             !{dir_output}/Run1/B1-RunGCcorrect/!{plateid_tumo}\
             !{dir_output}/Run1/C-SegmentBAF/!{plateid_tumo}\
             !{purity}\
             !{ploidy}\
             !{params.min_ploidy_full}\
             !{params.max_ploidy_full}\
             !{baseDir}/config/battenberg.config\
             !{baseDir}/data/imputation/impute_info.txt\
             !{dir_output}/Run2/G1-FitCopyNumber2/!{plateid_tumo}\
             !{dir_output}/Run2/G2-CallSubclones2/!{plateid_tumo}\
             !{dir_output}/!{plateid_tumo}_prior_breakpoints.tsv
    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// second DPClust run 
///////////////////////////////////////////////////////////////////////
process rerun1_stage2_dpclust {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv from fit2_ch1

    output:
    val dir_output into dpclust2_ch
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv into fit2_ch2
    file "*" optional true into dpclust2_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    if [ ! -s "!{dir_output}/Run2/H-DPClust2/!{plateid_tumo}_optimaInfo.txt" ]; then
        mkdir -p !{dir_output}/Run2/H-DPClust2

        Rscript !{baseDir}/scripts/DPClustPrep.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{sex}\
            !{filename_snv}\
            !{dir_output}/Run2/G1-FitCopyNumber2/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run2/G2-CallSubclones2/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run2/H-DPClust2/!{plateid_tumo}.loci.txt\
            !{dir_output}/Run2/H-DPClust2/!{plateid_tumo}.output.txt\
            !{dir_output}/Run2/H-DPClust2/!{plateid_tumo}.DPinput.txt

        Rscript !{baseDir}/scripts/DPClustRun.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{dir_output}/Run2/G1-FitCopyNumber2/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run2/H-DPClust2/!{plateid_tumo}.DPinput.txt\
            !{dir_output}/Run2/H-DPClust2

    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// assess peaks in copy-number-state-specific VAF distributions
///////////////////////////////////////////////////////////////////////
process rerun1_stage3_vaf_peaks {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv from fit2_ch2

    output:
    val dir_output into peaksvaf2_ch
    file "*" optional true into peaksvaf2_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    basename_output=!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData
    if [ ! -s "!{dir_output}/Run2/I-VAFPeaks2/${basename_output}" ]; then
        mkdir -p !{dir_output}/Run2/I-VAFPeaks2
        Rscript !{baseDir}/scripts/ComputePeaks_CNAqc.R\
            !{filename_snv}\
            !{dir_output}/Run2/H-DPClust2/!{plateid_tumo}.output.txt\
            !{dir_output}/Run2/G2-CallSubclones2/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run2/G2-CallSubclones2/!{plateid_tumo}_cellularity_ploidy.txt\
            !{baseDir}/config/peaksVAF.config\
            !{dir_output}/Run2/I-VAFPeaks2/!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData\
            !{dir_output}/Run2/I-VAFPeaks2\
            !{dir_output}/Run2/I-VAFPeaks2/!{plateid_tumo}_vafpeaks.png
    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// second QC run
///////////////////////////////////////////////////////////////////////
process rerun1_stage4_qc {
    errorStrategy 'terminate'
    publishDir "${dir_output_base}", mode: "copy", overwrite: true

    input:
    val x from dpclust2_ch.collect()
    val x from peaksvaf2_ch.collect()

    output:
    file "QualityControl_noBB.tsv" into qc2_table_file_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    cp !{dir_output_base}/QualityControl_noBB.tsv .
    Rscript !{baseDir}/scripts/QualityControl.R\
        2\
        !{dir_output_base}\
        /Run2/G2-CallSubclones2\
        /Run2/G2-CallSubclones2\
        /Run2/H-DPClust2\
        /Run2/I-VAFPeaks2\
        !{file_samplelist}\
        !{baseDir}/config/qc.config\
        !{baseDir}/data/hg19.chrom.sizes\
        QualityControl_noBB.tsv
    '''
}


///////////////////////////////////////////////////////////////////////
// SECOND RERUN
///////////////////////////////////////////////////////////////////////
// refit copy number in samples failing second round QC
qc2_table_file_ch
    .splitCsv(header: true, sep: '\t')
    .map { row -> [
        dir_output_base + "/intermediateFiles/" + row.participant_id + "/tumo" + row.tumour_sample_platekey + "_norm" + row.germline_sample_platekey,
        row.participant_id,
        row.sex,
        row.tumour_sample_platekey,
        row.germline_sample_platekey,
        row.run2_reestimated_purity,
        row.run2_reestimated_ploidy,
        row.somatic_small_variants_qced_vcf,
        row.run2_filter_overallfilter
    ] }
    .filter { it[8] == "FAIL" || it[8] == "FLAG" }
    .set { samplelist3_ch }

process rerun2_stage1_fit_copy_number {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, purity, ploidy, filename_snv, filter_overall from samplelist3_ch

    output:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv into fit3_ch1
    //set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv into fit3_ch2
    file "*" optional true into fit3_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b
    if [ ! -s "!{dir_output}/Run3/J2-CallSubclones3/!{plateid_tumo}_cellularity_ploidy.txt" ]; then
      mkdir -p !{dir_output}/Run3 !{dir_output}/Run3/J1-FitCopyNumber3 !{dir_output}/Run3/J2-CallSubclones3
      Rscript !{baseDir}/scripts/BattenbergFitCopyNumber.R\
             !{plateid_tumo}\
             !{sex}\
             !{dir_output}/Run1/B1-RunBAFLogR/!{plateid_tumo}\
             !{dir_output}/Run1/B1-RunGCcorrect/!{plateid_tumo}\
             !{dir_output}/Run1/C-SegmentBAF/!{plateid_tumo}\
             !{purity}\
             !{ploidy}\
             !{params.min_ploidy_full}\
             !{params.max_ploidy_full}\
             !{baseDir}/config/battenberg.config\
             !{baseDir}/data/imputation/impute_info.txt\
             !{dir_output}/Run3/J1-FitCopyNumber3/!{plateid_tumo}\
             !{dir_output}/Run3/J2-CallSubclones3/!{plateid_tumo}\
             !{dir_output}/!{plateid_tumo}_prior_breakpoints.tsv
    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// third DPClust run
///////////////////////////////////////////////////////////////////////
process rerun2_stage2_dpclust {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    //beforeScript { task.attempt<=1 ? 'sleep 0' : 'sleep 600' }
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv from fit3_ch1

    output:
    val dir_output into dpclust3_ch
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv into fit3_ch2
    file "*" optional true into dpclust3_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    if [ ! -s "!{dir_output}/Run3/K-DPClust3/!{plateid_tumo}_optimaInfo.txt" ]; then
        mkdir -p !{dir_output}/Run3/K-DPClust3

        Rscript !{baseDir}/scripts/DPClustPrep.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{sex}\
            !{filename_snv}\
            !{dir_output}/Run3/J1-FitCopyNumber3/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run3/J2-CallSubclones3/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run3/K-DPClust3/!{plateid_tumo}.loci.txt\
            !{dir_output}/Run3/K-DPClust3/!{plateid_tumo}.output.txt\
            !{dir_output}/Run3/K-DPClust3/!{plateid_tumo}.DPinput.txt


        Rscript !{baseDir}/scripts/DPClustRun.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{dir_output}/Run3/J1-FitCopyNumber3/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run3/K-DPClust3/!{plateid_tumo}.DPinput.txt\
            !{dir_output}/Run3/K-DPClust3

    fi
    '''

}


///////////////////////////////////////////////////////////////////////
// assess peaks in copy-number-state-specific VAF distributions
///////////////////////////////////////////////////////////////////////
process rerun2_stage3_vaf_peaks {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv from fit3_ch2

    output:
    val dir_output into peaksvaf3_ch
    file "*" optional true into peaksvaf3_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    basename_output=!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData
    if [ ! -s "!{dir_output}/Run3/L-VAFPeaks3/${basename_output}" ]; then
        mkdir -p !{dir_output}/Run3/L-VAFPeaks3
        Rscript !{baseDir}/scripts/ComputePeaks_CNAqc.R\
            !{filename_snv}\
            !{dir_output}/Run3/K-DPClust3/!{plateid_tumo}.output.txt\
            !{dir_output}/Run3/J2-CallSubclones3/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run3/J2-CallSubclones3/!{plateid_tumo}_cellularity_ploidy.txt\
            !{baseDir}/config/peaksVAF.config\
            !{dir_output}/Run3/L-VAFPeaks3/!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData\
            !{dir_output}/Run3/L-VAFPeaks3\
            !{dir_output}/Run3/L-VAFPeaks3/!{plateid_tumo}_vafpeaks.png
    fi
    '''
}

///////////////////////////////////////////////////////////////////////
// third QC run
///////////////////////////////////////////////////////////////////////
process rerun2_stage4_qc {
    errorStrategy 'terminate'
    publishDir "${dir_output_base}", mode: "copy", overwrite: true

    input:
    val x from dpclust3_ch.collect()
    val x from peaksvaf3_ch.collect()

    output:
    file "QualityControl_noBB.tsv" into qc3_table_file_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    cp !{dir_output_base}/QualityControl_noBB.tsv .
    Rscript !{baseDir}/scripts/QualityControl.R\
        3\
        !{dir_output_base}\
        /Run3/J2-CallSubclones3\
        /Run3/J2-CallSubclones3\
        /Run3/K-DPClust3\
        /Run3/L-VAFPeaks3\
        !{file_samplelist}\
        !{baseDir}/config/qc.config\
        !{baseDir}/data/hg19.chrom.sizes\
        QualityControl_noBB.tsv
    '''
}


///////////////////////////////////////////////////////////////////////
// THIRD RERUN
///////////////////////////////////////////////////////////////////////
// refit copy number of samples failing third round QC
qc3_table_file_ch
    .splitCsv(header: true, sep: '\t')
    .map { row -> [
        dir_output_base + "/intermediateFiles/" + row.participant_id + "/tumo" + row.tumour_sample_platekey + "_norm" + row.germline_sample_platekey,
        row.participant_id,
        row.sex,
        row.tumour_sample_platekey,
        row.germline_sample_platekey,
        row.run3_reestimated_purity,
        row.run3_reestimated_ploidy,
        row.somatic_small_variants_qced_vcf,
        row.run3_filter_overallfilter
    ] }
    .filter { it[8] == "FAIL" || it[8] == "FLAG" }
    .set { samplelist4_ch }

process rerun3_stage1_fit_copy_number {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, purity, ploidy, filename_snv, filter_overall from samplelist4_ch

    output:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv into fit4_ch1
    //set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv into fit4_ch2
    file "*" optional true into fit4_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    if [ ! -s "!{dir_output}/Run4/M2-CallSubclones4/!{plateid_tumo}_cellularity_ploidy.txt" ]; then
        mkdir -p !{dir_output}/Run4 !{dir_output}/Run4/M1-FitCopyNumber4 !{dir_output}/Run4/M2-CallSubclones4
        Rscript !{baseDir}/scripts/BattenbergFitCopyNumber.R\
             !{plateid_tumo}\
             !{sex}\
             !{dir_output}/Run1/B1-RunBAFLogR/!{plateid_tumo}\
             !{dir_output}/Run1/B1-RunGCcorrect/!{plateid_tumo}\
             !{dir_output}/Run1/C-SegmentBAF/!{plateid_tumo}\
             !{purity}\
             !{ploidy}\
             !{params.min_ploidy_full}\
             !{params.max_ploidy_full}\
             !{baseDir}/config/battenberg.config\
             !{baseDir}/data/imputation/impute_info.txt\
             !{dir_output}/Run4/M1-FitCopyNumber4/!{plateid_tumo}\
             !{dir_output}/Run4/M2-CallSubclones4/!{plateid_tumo}\
             !{dir_output}/!{plateid_tumo}_prior_breakpoints.tsv
    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// fourth DPClust run
///////////////////////////////////////////////////////////////////////
process rerun3_stage2_dpclust {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, sex, plateid_tumo, plateid_norm, filename_snv from fit4_ch1

    output:
    val dir_output into dpclust4_ch
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv into fit4_ch2
    file "*" optional true into dpclust4_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    if [ ! -s "!{dir_output}/Run4/N-DPClust4/!{plateid_tumo}_optimaInfo.txt" ]; then
        mkdir -p !{dir_output}/Run4/N-DPClust4

        Rscript !{baseDir}/scripts/DPClustPrep.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{sex}\
            !{filename_snv}\
            !{dir_output}/Run4/M1-FitCopyNumber4/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run4/M2-CallSubclones4/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run4/N-DPClust4/!{plateid_tumo}.loci.txt\
            !{dir_output}/Run4/N-DPClust4/!{plateid_tumo}.output.txt\
            !{dir_output}/Run4/N-DPClust4/!{plateid_tumo}.DPinput.txt

        Rscript !{baseDir}/scripts/DPClustRun.R\
            !{plateid_tumo}\
            !{plateid_norm}\
            !{dir_output}/Run4/M1-FitCopyNumber4/!{plateid_tumo}_rho_and_psi.txt\
            !{dir_output}/Run4/N-DPClust4/!{plateid_tumo}.DPinput.txt\
            !{dir_output}/Run4/N-DPClust4

    fi
    '''
}


///////////////////////////////////////////////////////////////////////
// assess peaks in copy-number-state-specific VAF distributions
///////////////////////////////////////////////////////////////////////
process rerun3_stage3_vaf_peaks {
    errorStrategy 'finish'
    tag "tumo:${plateid_tumo};norm:${plateid_norm}"
    publishDir "${dir_output}", mode: "copy", overwrite: true

    input:
    set dir_output, part_id, plateid_tumo, plateid_norm, filename_snv from fit4_ch2

    output:
    val dir_output into peaksvaf4_ch
    file "*" optional true into peaksvaf4_files_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    basename_output=!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData
    if [ ! -s "!{dir_output}/Run4/O-VAFPeaks4/${basename_output}" ]; then
        mkdir -p !{dir_output}/Run4/O-VAFPeaks4
        Rscript !{baseDir}/scripts/ComputePeaks_CNAqc.R\
            !{filename_snv}\
            !{dir_output}/Run4/N-DPClust4/!{plateid_tumo}.output.txt\
            !{dir_output}/Run4/M2-CallSubclones4/!{plateid_tumo}_subclones.txt\
            !{dir_output}/Run4/M2-CallSubclones4/!{plateid_tumo}_cellularity_ploidy.txt\
            !{baseDir}/config/peaksVAF.config\
            !{dir_output}/Run4/O-VAFPeaks4/!{part_id}.tumo!{plateid_tumo}_norm!{plateid_norm}.peak_data.RData\
            !{dir_output}/Run4/O-VAFPeaks4\
            !{dir_output}/Run4/O-VAFPeaks4/!{plateid_tumo}_vafpeaks.png
    fi
    '''

}


///////////////////////////////////////////////////////////////////////
// fourth QC run
///////////////////////////////////////////////////////////////////////
process rerun3_stage4_qc {
    errorStrategy 'terminate'
    publishDir "${dir_output_base}", mode: "copy", overwrite: true

    input:
    val x from dpclust4_ch.collect()
    val x from peaksvaf4_ch.collect()

    output:
    file "QualityControl_noBB.tsv" into qc4_table_file_ch

    shell:
    '''

    module load R/3.6.2-foss-2019b

    cp !{dir_output_base}/QualityControl_noBB.tsv .
    Rscript !{baseDir}/scripts/QualityControl.R\
        4\
        !{dir_output_base}\
        /Run4/M2-CallSubclones4\
        /Run4/M2-CallSubclones4\
        /Run4/N-DPClust4\
        /Run4/O-VAFPeaks4\
        !{file_samplelist}\
        !{baseDir}/config/qc.config\
        !{baseDir}/data/hg19.chrom.sizes\
        QualityControl_noBB.tsv
    '''

}
