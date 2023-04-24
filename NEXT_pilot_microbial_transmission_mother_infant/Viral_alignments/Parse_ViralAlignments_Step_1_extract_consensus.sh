### Bash script to extract consensus sequences from mapping of MGS and VLP BAM files 
### To the database of non-redundant viral contigs 
# Date: 21 Apr 2023
# By: Alex Kurilshikov
### Notes:
# requires samtools-1.17

### Loading input data

metadata = read.table("NEXT_Pilot_metadata.txt",as.is = T,header=T)
rpkm = read.table("RPKM_counts_combined_0.95_UPD_final_for_exp.txt",header=T,as.is =T)
mgs_files = dir("MGS_to_viral_decontaminated")
vlp_files = dir("VLP_to_viral_decontaminated")

### Extracting consensus sequences

for (i in 1:nrow(rpkm)){
dir.create(paste0("Consensus3/",rownames(rpkm)[i]))
  samples = colnames(rpkm)[which(rpkm[i,]>0)]
  print(paste(i , length(samples)))
  vlp = grep("V$",samples,value = T)
  mgs = grep("M$",samples,value = T)
  vlp_sam = paste0(metadata[match(vlp,metadata[,"Short_sample_ID"]),1],"_all_vir_alignments.sorted.bam")
  mgs_sam = paste0(metadata[match(mgs,metadata[,"Short_sample_ID"]),1],"_all_vir_alignments.sorted.bam")
  vlp_out = paste0("Consensus3/",
                    rownames(rpkm)[i],
                    "/",
                    metadata[match(vlp,metadata[,"Short_sample_ID"]),"FAM_ID"],
                    "_",
                    metadata[match(vlp,metadata[,"Short_sample_ID"]),"Type"],
                    "_",
                    sub("[0-9]+$","",   metadata[match(vlp,metadata[,"Short_sample_ID"]),"Individual_ID"]),
                    "_",
                    metadata[match(vlp,metadata[,"Short_sample_ID"]),"source"],
                    "_",
                    metadata[match(vlp,metadata[,"Short_sample_ID"]),"Timepoint"],
                    ".fa"
  )
  mgs_out = paste0("Consensus3/",
                   rownames(rpkm)[i],
                   "/",
                   metadata[match(mgs,metadata[,"Short_sample_ID"]),"FAM_ID"],
                   "_",
                   metadata[match(mgs,metadata[,"Short_sample_ID"]),"Type"],
                   "_",
                   sub("[0-9]+$","",   metadata[match(mgs,metadata[,"Short_sample_ID"]),"Individual_ID"]),
                   "_",
                   metadata[match(mgs,metadata[,"Short_sample_ID"]),"source"],
                   "_",
                   metadata[match(mgs,metadata[,"Short_sample_ID"]),"Timepoint"],
                   ".fa"
  )

  for (j in 1:length(vlp_sam)) {

    system(paste0("./samtools-1.17/samtools consensus -m simple -r ", rownames(rpkm)[i]," -o ", vlp_out[j]," VLP_to_viral_decontaminated/",vlp_sam[j]))
  }
  for (j in 1:length(mgs_sam)) {

    system(paste0("./samtools-1.17/samtools consensus -m simple -r ", rownames(rpkm)[i]," -o ", mgs_out[j]," MGS_to_viral_decontaminated/",mgs_sam[j]))
  }

}

### Merging consensus sequences to one file per contig

cd Consensus3

for j in `ls`; do cd ${j}
	for i in `ls *.fa`; do echo ">${i}"|perl -pe 's/[.]fa$//' >> merged.fa;tail -n+2 ${i}|perl -pe 's/N/-/g' >> merged.fa;done
	cd ../
	done
cd ../

