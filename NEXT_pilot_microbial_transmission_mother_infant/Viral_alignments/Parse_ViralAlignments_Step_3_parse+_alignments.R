### R script to make distance matrices from alignments
# Date: 21 Apr 2023
# By: Alex Kurilshikov
### Notes:
# The starting path is within folder Consensus3
# Some files are in aligned fasta format (*.fa extension), some in Clustal (*.aln). 
# The script distingiushes that and parses properly

library(ape)
files = dir("./all_aln/")
fa = grep("fa$",files,value = T)

seqs = list()
for(i in aln){
  print(i)
  seqs[[i]] = read.dna(paste0("./all_aln/",i),format = "clustal")
}

for(i in fa){
  seqs[[i]] = read.dna(paste0("./all_aln/",i),format = "fasta")
}

seqs.trimmed = lapply(seqs,function(x){
  length = dim(x)[2]
  x[,101:(length-100)]
})

dist.matrices = lapply(seqs.trimmed,function(x) dist.dna(x, model = "K80",pairwise.deletion = T,as.matrix = T))

names(dist.matrices) = sub("[.]aln.*","",names(dist.matrices))

lapply(names(dist.matrices),function(x) {
  write.table(dist.matrices[[x]],sep="\t",file = paste0("./dist.matrices/",x,".dist.txt"))
}
)
