############################ Making oordination plots for species of interest ########################################################
library(ggpubr)
library(ggplot2)
library(vegan)
library("gridExtra")

setwd("~/Resilio Sync/NEXT Pilot/04. Scripts /Strain_PCOA /Species_key_names")

# Making new key with descriptive ID's 
metadata = read.table(file = "~/Resilio Sync/NEXT Pilot/03. Phenotypes/Microbiome_metadata_phenos.txt", header = T, row.names = 1, check.names = FALSE, sep = "\t")
key=metadata
key$SAMPLE_ID=row.names(key)
key$Delivery_mode=gsub("Vaginal","VG",key$Delivery_mode)
key$Feeding_type_baby=gsub("Combination","CB",key$Feeding_type_baby)
key$Feeding_type_baby=gsub("Bottlefeeding","FF",key$Feeding_type_baby)
key$Feeding_type_baby=gsub("Breastfeeding","BF",key$Feeding_type_baby)
key$Feeding_type_baby[is.na(key$Feeding_type_baby)]="NA"
key$Place=gsub("Home", "HM", key$Place)
key$Place=gsub("Hospital", "HS", key$Place)
key$id=paste(key$new_sample_id,key$Delivery_mode,key$Feeding_type_baby,key$Place,sep = "_")
key=key[!duplicated(key$id),]
row.names(key)=key$id

distance=read.table("s_Bifidobacterium_bifidum.txt") # here load the distance matrix of species of interest 
data.b.pcoa=cmdscale(distance,k=(nrow(distance)-1),eig=TRUE)
combined<-merge(distance, key, by="row.names")
rownames(combined)<-combined$Row.names
combined$Row.names=NULL

# Adding first 2 PC's to combined 
pcoa = data.frame(PC1 = data.b.pcoa$points[,1], PC2 = data.b.pcoa$points[,2])
combinedPC=merge(pcoa, combined, by="row.names")
rownames(combinedPC)<-combinedPC$Row.names
combinedPC$Row.names=NULL
combinedPC<-combinedPC[!is.na(combinedPC$Type),]#Removing the sample for which type is NA
combinedPC$Type=gsub("Baby","Infant",combinedPC$Type)
combinedPC$Place=gsub("HM","Home",combinedPC$Place)
combinedPC$Place=gsub("HS","Hospital",combinedPC$Place)
#combinedPC<-combinedPC[combinedPC$Timepoint %in% c("B", "P7", "M1"), ] #Only choosing M B or P7 and Baby M1 

names (combinedPC)
key1<-key

#Check significance of place of birth on matrix of species
key1<-key1[rownames(distance),]
names(key1)
ad1<-adonis(distance ~ key1$Type+key1$Pair_ID+key1$Place,permutations=999,parallel=2)
ad1

#Individually looking at phenotypes 
names(combinedPC)

r=ggplot(combinedPC, aes(x=PC1, y=PC2, color=Place, shape=Type)) + geom_point(size=3) + theme_bw()
r
r1=r+ggtitle("Bifidobacterium bifidum")+geom_line(aes(group = Pair_ID), lty = 6, colour = "black")+
  #+geom_text(aes(label=Pair_ID),hjust=1, vjust=2,  size = 2.5)+
  scale_color_manual(values = c("#F62950","#005A9C"))+ #your colors here
  theme(
  plot.title = element_text(color="black", size=24, face="bold"),
  axis.title.x = element_text(color="black", size=24, face="bold"),
  axis.title.y = element_text(color="black", size=24, face="bold")
)
r1
r2=r1+ theme(legend.text = element_text(size=22))+theme(legend.title = element_blank())
r2


#To check CS 
names (combinedPC)
q=ggplot(combinedPC, aes(x=PC1, y=PC2, color=Delivery_mode, shape=Type)) + geom_point(size=3) + theme_bw()
q

q1=q+ggtitle("Bifidobacterium bifidum")+geom_line(aes(group = Pair_ID), lty = 6, colour = "black")+
  #+geom_text(aes(label=Pair_ID),hjust=1, vjust=2,  size = 2.5)+
  scale_color_manual(values = c("#F62950","#005A9C"))+ #your colors here
  theme(
    plot.title = element_text(color="black", size=24, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold")
  )
q1

#Akkeramancia 
distanceAK=read.table("s_Akkermansia_muciniphila.txt") # here load the distance matrix of species of interest 
data.b.pcoaAK=cmdscale(distanceAK,k=(nrow(distanceAK)-1),eig=TRUE)
combinedAK<-merge(distanceAK, key, by="row.names")
rownames(combinedAK)<-combinedAK$Row.names
combinedAK$Row.names=NULL

# Adding first 2 PC's to combined 
pcoaAK = data.frame(PC1 = data.b.pcoaAK$points[,1], PC2 = data.b.pcoaAK$points[,2])
combinedPCAK=merge(pcoaAK, combinedAK, by="row.names")
rownames(combinedPCAK)<-combinedPCAK$Row.names
combinedPCAK$Row.names=NULL
combinedPCAK<-combinedPCAK[!is.na(combinedPCAK$Type),] #Removing the sample for which type is NA
keep=combinedPCAK
keep$Family_ID
keep <-keep[keep$Family_ID %in% c("F05", "F06", "F27"), ]
keep$Type=gsub("Baby","Infant",keep$Type)
#combinedPC<-combinedPC[combinedPC$Timepoint %in% c("B", "P7", "M1"), ] #Only choosing M B or P7 and Baby M1 


#Individually looking at phenotypes 

a=ggplot(keep, aes(x=PC1, y=PC2, color=Type)) + geom_point(size=3) + 
  theme_bw()
a

a1=a+ggtitle("Akkermansia muciniphila")+geom_text(aes(label=Timepoint),hjust=0.5, vjust=2)+geom_line(aes(group = Pair_ID), lty = 6, colour = "black")+
theme(
  plot.title = element_text(color="black", size=24, face="bold"),
  axis.title.x = element_text(color="black", size=24, face="bold"),
  axis.title.y = element_text(color="black", size=24, face="bold")
)
a1
a2=a1+ theme(legend.text = element_text(size=22))+theme(legend.title = element_blank())
a2





grid.arrange(r2, a2,ncol = 1, nrow = 2)

dev.off()

