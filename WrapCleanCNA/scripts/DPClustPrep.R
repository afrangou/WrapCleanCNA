if (sys.nframe() == 0L) {
  args <- commandArgs(T)
  print(args)
  
  tumourplatekey=args[1]
  normalplatekey=args[2]
  gender=args[3]
  vcffilepath=args[4]
  rhoandpsifilepath=args[5]
  subcloneshg38filepath=args[6]
  output_loci_file=args[7]
  output_file=args[8]
  output_DPinput_file=args[9]
  
  suppressPackageStartupMessages(library(VariantAnnotation))
  suppressPackageStartupMessages(library(dpclust3p))
  print(tumourplatekey)
  print(normalplatekey)
  print(gender)
  print(vcffilepath)
  print(rhoandpsifilepath)
  print(subcloneshg38filepath)
  print(output_loci_file)
  print(output_file)
  print(output_DPinput_file)
  nucleotides = c("A", "C", "G", "T")
  samplename = paste0("tumo", tumourplatekey, "_norm", normalplatekey)
  print(paste("Sample name:", samplename))
  print(paste("Using", subcloneshg38filepath))
  stopifnot(file.exists(vcffilepath))
  
  print(paste("VCF is ", vcffilepath))
  snvs=readVcf(vcffilepath)
  print(paste("vcfloaded", vcffilepath))
  snvs=data.frame(chr=as.character(seqnames(snvs)),
                  pos=start(snvs),
                  ref=ref(snvs),
                  alt=as.character(unlist(alt(snvs))),
                  c_ref=as.numeric(info(snvs)$t_ref_count),
                  c_alt=as.numeric(info(snvs)$t_alt_count),
                  stringsAsFactors=F)
  snvs=snvs[which(snvs$chr %in% 1:22),]
  snvs=snvs[which(snvs$c_ref+snvs$c_alt>=10),]
  
  tumpos=data.frame(Chromosome=snvs$chr,
                    Position=snvs$pos,
                    REF=snvs$ref,
                    ALT=snvs$alt,
                    count_A=0,
                    count_C=0,
                    count_G=0,
                    count_T=0,
                    total_depth=snvs$c_ref+snvs$c_alt)
  
  CONVERT=c(5:8)
  names(CONVERT)=c('A','C','G','T')
  
  tumpos[cbind(1:nrow(snvs),CONVERT[snvs$ref])]=snvs$c_ref
  tumpos[cbind(1:nrow(snvs),CONVERT[snvs$alt])]=snvs$c_alt
  rm(CONVERT,snvs)
  
  if (nrow(tumpos)>50000) {
    set.seed(1234)
    TO_KEEP=sort(sample(1:nrow(tumpos),50000))
    print(paste("Writing all tumour position counts to ", gsub('.txt$','.all.txt',output_file)))
    write.table(tumpos[, c(1, 2, 5:9)], gsub('.txt$','.all.txt',output_file), row.names = F, 
                quote = F, sep = "\t")
    print(paste("Writing subsampled tumour position counts to ", output_file))
    write.table(tumpos[TO_KEEP, c(1, 2, 5:9)], output_file, row.names = F, 
                quote = F, sep = "\t")
    print(paste("Writing all tumour loci to ", gsub('.txt$','.all.txt',output_loci_file)))
    write.table(tumpos[, c(1:4)], gsub('.txt$','.all.txt',output_loci_file), row.names = F, 
                col.names = F, quote = F, sep = "\t")
    print(paste("Writing subsampled tumour loci to ", output_loci_file))
    write.table(tumpos[TO_KEEP, c(1:4)], output_loci_file, row.names = F, 
                col.names = F, quote = F, sep = "\t")
  } else {
    print(paste("Writing tumour position counts to ", output_file))
    write.table(tumpos[, c(1, 2, 5:9)], output_file, row.names = F, 
                quote = F, sep = "\t")
    print(paste("Writing tumour loci to ", output_loci_file))
    write.table(tumpos[, c(1:4)], output_loci_file, row.names = F, 
                col.names = F, quote = F, sep = "\t")
  }
  
  loci_file = output_loci_file
  allele_frequencies_file = output_file
  cellularity_file = rhoandpsifilepath
  subclone_file = subcloneshg38filepath
  output_DPinput_file = output_DPinput_file
  print(paste("Using loci file", loci_file))
  print(paste("Using allele frequencies file", allele_frequencies_file))
  print(paste("Using cellularity file", cellularity_file))
  print(paste("Using subclones file", subclone_file))
  print(paste("Writing to output file", output_DPinput_file))
  if (gender == "FEMALE") {
    gender = "female"
  } else if (gender == "MALE") {
    gender = "male"
  } else {
    print("Gender unspecified, exiting")
    break
  }
  runGetDirichletProcessInfo(loci_file, allele_frequencies_file, 
                             cellularity_file, subclone_file, gender, SNP.phase.file = "NA", 
                             mut.phase.file = "NA", output_file = output_DPinput_file)
}