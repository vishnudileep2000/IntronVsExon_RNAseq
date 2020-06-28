####Created on April 9th, 2020
###author:Vishnu Dileep
library(data.table)
library(parallel)
library(travis) #https://github.com/dvera/travis#


#####make intron exon file###
splitgene=function(x){
  
  exon=NULL;intron=NULL
  exon$start=unlist(strsplit(as.vector(x[10]),","))
  exon$end=unlist(strsplit(as.vector(x[11]),","))
  exon=lapply(exon,as.numeric)
  exon=as.data.frame(exon)
  
  exon=cbind(CHR=rep(x[3], x[9]),exon)
  exon$gname=rep(x[2], x[9])
  
  intron=data.frame(CHR=exon[,1])
  intron$start=exon$end
  intron$end=c(exon$start[-1], x[6])
  intron$gname=exon$gname
  
  return(cbind(exon,intron))
  
}

grn=read.delim("Refseq_mm9_formatted.txt", header=T)
grn=grn[!duplicated(grn$name2),]

ex.int=do.call("rbind",apply(grn,1,splitgene))

for(i in c(1,4,5,8)) ex.int[,i]=as.character(as.vector(ex.int[,i]))
for(i in c(2,3,6,7)) ex.int[,i]=as.integer(as.vector(ex.int[,i]))

exon=ex.int[,1:4]
intron=ex.int[,5:8]

intron$size=intron$end-intron$start
exon$size=exon$end-exon$start

sum(intron$size==0)
sum(exon$size==0)

intron$start2=intron$start+10
intron$end2=intron$end-10
intron$size2=intron$end2-intron$start2
sum(intron$size2<0)

intron[intron$size2<=0,c("start2")]=intron[intron$size2<=0,c("start")]
intron[intron$size2<=0,c("end2")]=intron[intron$size2<=0,c("start")]
intron$size2=intron$end2-intron$start2
sum(intron$size2==0)

intron=intron[,c("CHR","start2","end2","gname","size2")]
names(intron)=c("CHR","start","end","gname","size")


###collapse and aggreagte intron and exon ranges###
collapse.exon=function(x){
  sub=subset(exon,gname==x)
  ranges=paste0(sub$CHR[1],":",paste0(sub$start,"-",sub$end),collapse=",")
  return(ranges)
  }

collapse.intron=function(x){
  sub=subset(intron,gname==x)
  sub=sub[sub$size>0,]
  ranges=paste0(sub$CHR[1],":",paste0(sub$start,"-",sub$end),collapse=",")
  return(ranges)
}


exon.ranges=unlist(mclapply(unique(exon$gname), collapse.exon, mc.cores = 30))
intron.ranges=unlist(mclapply(unique(intron$gname), collapse.intron, mc.cores = 30))

ex.int.ranges=data.frame(name=unique(exon$gname),exon.ranges,intron.ranges)

##aggregate sizes###

aggExon.size=aggregate(size~gname,exon, sum)
names(aggExon.size)=c("name","exonSize")

aggIntron.size=aggregate(size~gname,intron, sum)
names(aggIntron.size)=c("name","intronSize")

exon$size=NULL;intron$size=NULL



write.table(exon,"ExonLocs.bed",sep="\t",row.names = F, col.names = F, quote=F)
write.table(intron,"IntronLocs.bed",sep="\t",row.names = F, col.names = F, quote=F)





#########convert BAM to BED####
setwd("/media/vishnu/Gangsta/Other/marco/intron_analysis/BED")
Bams= files("../BAM/*.bam")
set.cores=30
Beds=bamToBed(Bams, threads=set.cores)


###Format beds###
Beds=files("*.bed")

for(thisbed in Beds){
cat(thisbed,"\n")
bedfile=fread(thisbed, header=T)
bedfile[,c(4:6):=NULL]

newBed=unlist(strsplit(thisbed,split = "_a"))[1]
newBed=paste0(newBed,".bed")
fwrite(bedfile,newBed,quote = F, sep = "\t", row.names = F, col.names = F)
newBed=NULL;bedfile=NULL
}


########split bed file into exon and intron and count###
###deleted old bedfiles###
setwd("/media/vishnu/Gangsta/Other/marco/intron_analysis/BED")
Beds=files("*.bed")
setwd("/media/vishnu/Gangsta/Other/marco/intron_analysis/intersects")


for(thisbed in Beds) {
Ex.outfilename=paste0(unlist(strsplit(thisbed,split = ".bed"))[1],
                   "_exon.intersect.txt")
Int.outfilename=paste0(unlist(strsplit(thisbed,split = ".bed"))[1],
                       "_intron.intersect.txt")

cmd=paste0("bedtools intersect -a ../ExonLocs.bed -b ../BED/",thisbed," -c > ",Ex.outfilename,
           ";bedtools intersect -a ../IntronLocs.bed -b ../BED/",thisbed," -c > ",Int.outfilename)
cat(cmd,"\n")
system(cmd)
}



#####combine counts###
setwd("/media/vishnu/Gangsta/Other/marco/intron_analysis/BED")
Beds=files("*.bed")
setwd("/media/vishnu/Gangsta/Other/marco/intron_analysis/intersects")


finalmerge=NA
for(thisbed in Beds){
cat(thisbed,"\n")
exon=NULL;intron=NULL
sample=unlist(strsplit(thisbed,split = ".bed"))[1]
Ex.filename=paste0(sample,"_exon.intersect.txt")
Int.filename=paste0(sample, "_intron.intersect.txt")

exon=read.delim(Ex.filename, header=F)
intron=read.delim(Int.filename, header=F)

names(exon)=c("CHR","Start","End","name","ex.reads")
names(intron)=c("CHR","Start","End","name","int.reads")


exon.sum=aggregate(ex.reads~name,exon, sum)
intron.sum=aggregate(int.reads~name,intron, sum)

sum.ex.int=merge(exon.sum,intron.sum,by=c("name"))
names(sum.ex.int)=c("name",paste0(sample,".ex.reads"),paste0(sample,".int.reads"))

finalmerge=cbind(finalmerge,sum.ex.int)
}

names(finalmerge)[c(1,seq(5,ncol(finalmerge),3))]
finalmerge[,c(1,seq(5,37,3))]=NULL


####merge other data##

finalmerge=merge(finalmerge,aggExon.size, by="name")
finalmerge=merge(finalmerge,aggIntron.size, by="name")

names(grn)[2]="name"

finalmerge=merge(finalmerge,grn[,1:9], by="name")
finalmerge=merge(finalmerge,ex.int.ranges, by="name")

finalmerge=finalmerge[,c(1,28:35,26,27,2:25,36,37)]

write.table(finalmerge,"Final_exon_intron_rawreads_allSamples.txt",sep="\t",row.names = F, col.names = T, quote=F)


####normalization####

dat=finalmerge

###RPM####
dat[,12:35]=sweep(dat[,12:35],2,colSums(dat[,12:35]),'/') *1e6

###per kil0base##
dat[,seq(12,35,2)]=dat[,seq(12,35,2)]/dat$exonSize *1000
dat[,seq(13,35,2)]=dat[,seq(13,35,2)]/dat$intronSize *1000

write.table(dat,"Final_exon_intron_RPKM_allSamples.txt",sep="\t",row.names = F, col.names = T, quote=F)







###misc###

i=15438;play=0

exon=NULL
exon$start=unlist(strsplit(as.vector(grn$exonStarts[i]),","))
exon$end=unlist(strsplit(as.vector(grn$exonEnds[i]),","))
exon=lapply(exon,as.numeric)
exon=as.data.frame(exon)

exon=cbind(CHR=rep(grn$chrom[i], grn$exonCount[i]),exon)
exon$gname=rep(grn$name2[i], grn$exonCount[i])

intron=data.frame(CHR=exon[,1])
intron$start=exon$end+play
intron$end=c(exon$start[-1]-play, grn$txEnd[i])
intron$gname=exon$gname