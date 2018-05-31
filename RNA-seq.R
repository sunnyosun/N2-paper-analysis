# RNAseq after cuffdiff

# making vocano plot
res <- read.table("gene_exp.diff", header=TRUE,stringsAsFactors = F,sep='\t')
res=res[which(res[,'status']=='OK'),]
res=res[which(res[,'log2.fold_change.']!='-Inf'),]
res=res[which(res[,'log2.fold_change.']!='Inf'),]

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval), pch=20,main="Volcano plot"))

with(subset(res, significant=='yes'), points(log2FoldChange, -log10(pval), pch=20, col="red"))


# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, p_value<.00001 ), points(log2.fold_change., -log10(p_value), pch=20, col="red"))
with(subset(res, abs(log2.fold_change.)>2), points(log2.fold_change., -log10(p_value), pch=20, col="orange"))
with(subset(res, pval<.00001 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="green"))

sig=subset(res, p_value<.01 & abs(log2.fold_change.)>2)
write.table(sig,'pvalue001_logfc2.txt',quote=F,sep='\t',col.names = T,row.names = F)

gene=read.table('gene_name.txt',header=T,sep='\t',stringsAsFactors = F)
for (i in 1:nrow(sig))
{
  sig[i,15]=gene[which(gene[,'id']==sig[i,'gene_id']),'name']
}
colnames(sig)[15]='name'

with(subset(sig), textxy(log2.fold_change., -log10(p_value), labs=name, cex=.8))


#######################
# RNAseq after htseq-count
data1=read.table('BS00295A_S7.txt',sep='\t',stringsAsFactors = F)
data2=read.table('BS00296A_S8.txt',sep='\t',stringsAsFactors = F)
data3=read.table('BS00297A_S9.txt',sep='\t',stringsAsFactors = F)
data4=read.table('BS00298A_S10.txt',sep='\t',stringsAsFactors = F)
data5=read.table('BS00299A_S11.txt',sep='\t',stringsAsFactors = F)
data6=read.table('BS00300A_S12.txt',sep='\t',stringsAsFactors = F)
data7=read.table('BS00301A_S13.txt',sep='\t',stringsAsFactors = F)
data8=read.table('BS00302A_S14.txt',sep='\t',stringsAsFactors = F)
data9=read.table('BS00303A_S15.txt',sep='\t',stringsAsFactors = F)
data10=read.table('BS00304A_S16.txt',sep='\t',stringsAsFactors = F)
data11=read.table('BS00305A_S17.txt',sep='\t',stringsAsFactors = F)
data12=read.table('BS00306A_S18.txt',sep='\t',stringsAsFactors = F)

n16=merge(data1,data2,by='V1')
n16=merge(n16,data3,by='V1')
colnames(n16)=c('gene','n16_1','n16_2','n16_3')
n16=n16[5:nrow(n16),]

n8=merge(data4,data5,by='V1')
n8=merge(n8,data6,by='V1')
colnames(n8)=c('gene','n8_1','n8_2','n8_3')
n8=n8[5:nrow(n8),]

n4=merge(data7,data8,by='V1')
n4=merge(n4,data9,by='V1')
colnames(n4)=c('gene','n4_1','n4_2','n4_3')
n4=n4[5:nrow(n4),]

n2=merge(data10,data11,by='V1')
n2=merge(n2,data12,by='V1')
colnames(n2)=c('gene','n2_1','n2_2','n2_3')
n2=n2[5:nrow(n2),]

data=merge(n16,n8,by='gene')
data=merge(data,n4,by='gene')
data=merge(data,n2,by='gene')

write.table(data,'htseq_count_merged.txt',quote=F,col.names=T,row.names=F,sep='\t')

##
ann=read.table('annotation.txt',sep='\t',stringsAsFactors = F,header=T)
gene=read.table('BY4741_gene_name.txt',sep='\t',stringsAsFactors = F)
c=row.names(Counts)
v=vector()
for (i in 1:length(c))
{
  tmp=gene[which(gene[,9]==c[i]),10]
  v[i]=gene[which(gene[,9]==tmp),10]
}

#####
Counts = read.table( 'htseq_count_merged.txt', header=TRUE, row.names=1,sep='\t',stringsAsFactors = F )

Design = data.frame(row.names = colnames(Counts),condition = c('n16','n16','n16','n8','n8','n8','n4','n4','n4','n2','n2','n2' ),libType = c("single-end", "single-end", "single-end","single-end", "single-end", "single-end", "single-end","single-end","single-end","single-end","single-end","single-end" ))

singleSamples =Design$libType == "single-end"
countTable = Counts[ , singleSamples ]
condition = Design$condition[ singleSamples ]
library( "DESeq" )
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
str( fitInfo(cds) )
plotDispEsts( cds )

### Diffrential expression gene
res = nbinomTest( cds, "n16", "n2" )
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
resSig = res[ res$padj < 0.1, ]
# most significantly differentially expressed genes:
head( resSig[ order(resSig$pval), ] )

# most strongly down-regulated of the significant genes
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )

# most strongly up-regulated ones:
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )

write.csv( res, file="n2n16.csv" )


#info=read.table('n2_info.txt',header=T,sep='\t',stringsAsFactors = F)

# making vocano plot
res <- read.csv("n2n16.csv", header=TRUE,stringsAsFactors = F,row.names = 1)

# filter out deleted genes
info=read.table('n2_info.txt',header=T,sep='\t',stringsAsFactors = F)
name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')
for (i in 1:nrow(info))
{
  chr=name[which(name[,1]==info[i,'chr'])[1],1]
  start=max(name[which(name[,9]==info[i,'start_gene']),4:5])
  end=min(name[which(name[,9]==info[i,'end_gene']),4:5])
  info[i,8]=start
  info[i,9]=end
}

res2=res[1,]
for (i in 2:nrow(res))
{
  mi=min(name[which(name[,9]==res[i,1]),4:5])
  ma=max(name[which(name[,9]==res[i,1]),4:5])
  chr=name[which(name[,9]==res[i,1]),1]
  if (chr=='chrM' | chr=='2micron')
  {
    res2=rbind(res2,res[i,])
  }
  else
  {
    boun=info[which(info[,1]==chr),8:9]
    if (mi>boun[1] & ma<boun[2])
    {
      res2=rbind(res2,res[i,])
    }
  }
}
res=res2
write.table(res,'n2n16_norm_clean.txt',quote = F, sep='\t',col.names=T,row.names=F)

#change 0 to 0.1
n2=read.table('n2n16_norm_clean.txt',sep='\t',header=T,stringsAsFactors = F)
n2=n2[!is.na(n2[,6]),]
n2[which(n2[,'log2FoldChange']=='-Inf'),'log2FoldChange']=log2(0.1/n2[which(n2[,'log2FoldChange']=='-Inf'),'baseMeanA'])
n2[which(n2[,'log2FoldChange']=='Inf'),'log2FoldChange']=log2(n2[which(n2[,'log2FoldChange']=='Inf'),'baseMeanB']/0.1)
res=n2

setEPS()
postscript("n2n16_vp_crop_1e-5.eps",width=6, height=5)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval), ylim=c(0,12),pch=20,main="Volcano plot",col='grey'))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
s=subset(res, pval<.00001 & abs(log2FoldChange)>2)
with(s, points(log2FoldChange, -log10(pval), pch=20, col="red"))
s2=subset(s, baseMeanB==0 | baseMeanA==0)
with(s2, points(log2FoldChange, -log10(pval), pch=20, col="blue"))


na=vector()
for (i in 1:nrow(s))
{
  na[i]=strsplit(s[i,1],'_')[[1]][1]
}
name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')

for (i in 1:length(na))
{
  na[i]=name[which(name[,'V9']==na[i]),'V10']
}

with(s,text(log2FoldChange, -log10(pval)+0.8, labels=na, cex=0.8))
dev.off()



######################## plot 4 section dots
######################
n2=read.table('n8_sig.txt',header=T,sep='\t',stringsAsFactors = F)
n2[which(n2[,16]==-1),17]<-'blue'
n2[which(n2[,16]==0),17]<-'grey'
n2[which(n2[,16]==1),17]<-'red'

setEPS()
postscript("n2_sig.eps",height=4,width=6,paper='special')
par(mfrow=c(1,1),lab=c(x=7,y=5,len=1),las=1,tcl=-0.25,mgp=c(1.75,0.3,0),
    bty='n',mar=c(2.5,2,0.5,0.2),cex.lab=1.5,cex.axis=0.8,font.lab=1,cex=1,
    omi=c(0,0,0,0),xaxt='n')
plot(n2[,c('fused_1','log2.fold_change.')],pch=16,col=n2[,17],ylim=c(-6,6),xlim=c(-1.5,1.5),
     ylab='Log2 (fold change of expression)')
abline(h=0)
abline(v=0)
text(n2[,'fused_1']+0.1, n2[,'log2.fold_change.'], labels=n2[,'name'], cex=0.8)
dev.off()

setEPS()
postscript("n4_sig.eps",height=4,width=6,paper='special')
par(mfrow=c(1,1),lab=c(x=7,y=5,len=1),las=1,tcl=-0.25,mgp=c(1.75,0.3,0),
    bty='n',mar=c(2.5,2,0.5,0.2),cex.lab=1.5,cex.axis=0.8,font.lab=1,cex=1,
    omi=c(0,0,0,0),xaxt='n')
plot(n2[,c('fused_1','log2.fold_change.')],pch=16,col=n2[,17],ylim=c(-6,6),xlim=c(-1.5,1.5),
     ylab='Log2 (fold change of expression)')
abline(h=0)
abline(v=0)
text(n2[,'fused_1']+0.1, n2[,'log2.fold_change.'], labels=n2[,'name'], cex=0.8)
dev.off()


setEPS()
postscript("n8_sig.eps",height=4,width=6,paper='special')
par(mfrow=c(1,1),lab=c(x=7,y=5,len=1),las=1,tcl=-0.25,mgp=c(1.75,0.3,0),
    bty='n',mar=c(2.5,2,0.5,0.2),cex.lab=1.5,cex.axis=0.8,font.lab=1,cex=1,
    omi=c(0,0,0,0),xaxt='n')
plot(n2[,c('fuse_1','log2.fold_change.')],pch=16,col=n2[,17],ylim=c(-6,6),xlim=c(-1.5,1.5),
     ylab='Log2 (fold change of expression)')
abline(h=0)
abline(v=0)
text(n2[,'fuse_1']+0.1, n2[,'log2.fold_change.'], labels=n2[,'name'], cex=0.8)
dev.off()


### plot 4 section dots vs chr locations
len=read.table('chr_length_sacCer2.txt',sep='\t',stringsAsFactors = F)

n2=read.table('n8_sig.txt',header=T,sep='\t',stringsAsFactors = F)
n2[which(n2[,16]==-1),17]<-'blue'
n2[which(n2[,16]==0),17]<-'grey'
n2[which(n2[,16]==1),17]<-'red'
for (i in 1:nrow(n2))
{
  chr=strsplit(n2[i,'locus'],split = ':')[[1]][1]
  tmp=strsplit(n2[i,'locus'],split = ':')[[1]][2]
  start=as.numeric(strsplit(tmp,split = '-')[[1]][1])
  end=as.numeric(strsplit(tmp,split = '-')[[1]][2])
  l=len[which(len[,1]==chr),2]
  if (start>50000)
  {
    loc=l-start
  }
  else
  {
    loc=start
  }
  n2[i,18]=loc
  if (n2[i,17]=='grey')
  {
    n2[i,18]=0
  }
  if (n2[i,17]=='blue')
  {
    n2[i,18]=-(50000-loc)
  }
  if (n2[i,17]=='red')
  {
    n2[i,18]=50000-loc
  }
}

setEPS()
postscript("n8_sig.eps",height=4,width=3,paper='special')
par(mfrow=c(1,1),lab=c(x=5,y=5,len=1),las=1,tcl=-0.25,mgp=c(1.75,0.3,0),
    bty='n',mar=c(2.5,2,0.5,0.2),cex.lab=1.5,cex.axis=0.8,font.lab=1,cex=1,
    omi=c(0,0,0,0))
plot(n2[,'V18'],n2[,'log2.fold_change.'],pch=16,col=n2[,17],ylim=c(-6,6),xlim=c(-50000,50000),
     ylab='',xlab='')
abline(h=0)
abline(v=0)
text(n2[,'V18']+5000, n2[,'log2.fold_change.'], labels=n2[,'name'], cex=0.8)
dev.off()




####################################################
#n=4

### Diffrential expression gene
res = nbinomTest( cds, "n16", "n4" )
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
resSig = res[ res$padj < 0.1, ]
# most significantly differentially expressed genes:
head( resSig[ order(resSig$pval), ] )

# most strongly down-regulated of the significant genes
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )

# most strongly up-regulated ones:
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )

write.csv( res, file="n4n16.csv" )


#info=read.table('n2_info.txt',header=T,sep='\t',stringsAsFactors = F)

# making vocano plot
res <- read.csv("n4n16.csv", header=TRUE,stringsAsFactors = F,row.names = 1)
# filter out deleted genes
info=read.table('n4_info.txt',header=T,sep='\t',stringsAsFactors = F)
name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')
for (i in 1:nrow(info))
{
  chr=name[which(name[,1]==info[i,'chr'])[1],1]
  start=max(name[which(name[,9]==info[i,'start_gene']),4:5])
  end=min(name[which(name[,9]==info[i,'end_gene']),4:5])
  info[i,8]=start
  info[i,9]=end
}

res2=res[1,]
for (i in 2:nrow(res))
{
  mi=min(name[which(name[,9]==res[i,1]),4:5])
  ma=max(name[which(name[,9]==res[i,1]),4:5])
  chr=name[which(name[,9]==res[i,1]),1]
  if (chr=='chrM' | chr=='2micron')
  {
    res2=rbind(res2,res[i,])
  }
  else
  {
    boun=info[which(info[,1]==chr),8:9]
    if (mi>boun[1] & ma<boun[2])
    {
      res2=rbind(res2,res[i,])
    }
  }
}
res=res2
write.table(res,'n4n16_norm_clean.txt',quote = F, sep='\t',col.names=T,row.names=F)

#change 0 to 0.1
n2=read.table('n4n16_norm_clean.txt',sep='\t',header=T,stringsAsFactors = F)
n2=n2[!is.na(n2[,6]),]
n2[which(n2[,'log2FoldChange']=='-Inf'),'log2FoldChange']=log2(0.1/n2[which(n2[,'log2FoldChange']=='-Inf'),'baseMeanA'])
n2[which(n2[,'log2FoldChange']=='Inf'),'log2FoldChange']=log2(n2[which(n2[,'log2FoldChange']=='Inf'),'baseMeanB']/0.1)
res=n2

setEPS()
postscript("n4n16_vp_all.eps",width=6, height=5)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval), pch=20,col='grey',main="Volcano plot"))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
s=subset(res, pval<.00001 & abs(log2FoldChange)>2)
s=s[which(s[,'baseMeanB']!=0),]
with(s, points(log2FoldChange, -log10(pval), pch=20, col="red"))

na=vector()
for (i in 1:nrow(s))
{
  na[i]=strsplit(s[i,1],'_')[[1]][1]
}
name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')

for (i in 1:length(na))
{
  na[i]=name[which(name[,'V9']==na[i]),'V10']
}

with(s,text(log2FoldChange, -log10(pval)+5, labels=na, cex=1))
dev.off()



###################
# n=8
res = nbinomTest( cds, "n16", "n8" )
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
resSig = res[ res$padj < 0.1, ]
# most significantly differentially expressed genes:
head( resSig[ order(resSig$pval), ] )

# most strongly down-regulated of the significant genes
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )

# most strongly up-regulated ones:
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )

write.csv( res, file="n8n16.csv" )

# making vocano plot
res <- read.csv("n8n16.csv", header=TRUE,stringsAsFactors = F,row.names = 1)

# filter out deleted genes
info=read.table('n8_info.txt',header=T,sep='\t',stringsAsFactors = F)
name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')
for (i in 1:nrow(info))
{
  chr=name[which(name[,1]==info[i,'chr'])[1],1]
  start=max(name[which(name[,9]==info[i,'start_gene']),4:5])
  end=min(name[which(name[,9]==info[i,'end_gene']),4:5])
  info[i,8]=start
  info[i,9]=end
}

res2=res[1,]
for (i in 2:nrow(res))
{
  mi=min(name[which(name[,9]==res[i,1]),4:5])
  ma=max(name[which(name[,9]==res[i,1]),4:5])
  chr=name[which(name[,9]==res[i,1]),1]
  if (chr=='chrM' | chr=='2micron')
  {
    res2=rbind(res2,res[i,])
  }
  else
  {
    boun=info[which(info[,1]==chr),8:9]
    if (mi>boun[1] & ma<boun[2])
    {
      res2=rbind(res2,res[i,])
    }
  }
}
res=res2
write.table(res,'n8n16_norm_clean.txt',quote = F, sep='\t',col.names=T,row.names=F)


setEPS()
postscript("n8n16_vp_all.eps",width=10, height=5)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval),ylim=c(0,150), pch=20,main="Volcano plot",col='grey'))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
s=subset(res, pval<.00001 & abs(log2FoldChange)>2)
s=s[which(s[,'baseMeanB']!=0),]
with(s, points(log2FoldChange, -log10(pval), pch=20, col="red"))

na=vector()
for (i in 1:nrow(s))
{
  na[i]=strsplit(s[i,1],'_')[[1]][1]
}
name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')

for (i in 1:length(na))
{
  na[i]=name[which(name[,'V9']==na[i]),'V10']
}

with(s,text(log2FoldChange, -log10(pval)+5, labels=na, cex=1))
dev.off()


#################################
# use cuffdiff output
setwd("/Volumes/My Book for Mac/JL_20161221")
difffile='Jingchuan/N16N8/gene_exp.diff'
res=read.table(difffile,header=T,sep='\t',stringsAsFactors = F)

# remove the deleted genes due to fusion
# filter out deleted genes
infoname='Jingchuan/n8_info.txt'
outfile='Jingchuan/N16N8/n8gene_exp_clean.diff'

info=read.table(infoname,header=T,sep='\t',stringsAsFactors = F)
name=read.table('Jingchuan/BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')
for (i in 1:nrow(info))
{
  chr=name[which(name[,1]==info[i,'chr'])[1],1]
  start=min(name[which(name[,9]==info[i,'start_gene']),4:5])
  end=max(name[which(name[,9]==info[i,'end_gene']),4:5])
  info[i,8]=start
  info[i,9]=end
}

res2=res[1,]
for (i in 2:nrow(res))
{
  mi=min(name[which(name[,9]==res[i,1]),4:5])
  ma=max(name[which(name[,9]==res[i,1]),4:5])
  chr=name[which(name[,9]==res[i,1]),1]
  if (chr=='chrM' | chr=='2micron')
  {
    res2=rbind(res2,res[i,])
  }
  else
  {
    boun=info[which(info[,1]==chr),8:9]
    if (mi>=boun[1] & ma<=boun[2])
    {
      res2=rbind(res2,res[i,])
    }
  }
}
res=res2
write.table(res,outfile,quote = F, sep='\t',col.names=T,row.names=F)

####### volcano plot
# read clean data
setwd("/Volumes/My Book for Mac/JL_20161221")
cleanfile='Jingchuan/N16N8/n8gene_exp_clean.diff'
res=read.table(cleanfile,header=T,sep='\t',stringsAsFactors = F)
name=read.table('Jingchuan/BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')

# remove 2micron genes
genes_2m=subset(name,V1=='2micron')[,10]
res2=subset(res,!(res$gene %in% genes_2m))
res=res2

# change the pvalue=0 to p
p=unique(res[order(res$p_value),'p_value'])[2]
res[which(res$p_value==0),'p_value']=p
res2=subset(res,status=='OK')
res=res2
res$value_1=res$value_1+0.1
res$value_2=res$value_2+0.1

# re-calculate the fold change
res$log2.fold_change.=log2(res$value_2/res$value_1)

# plot
with(res, plot(log2.fold_change., -log10(p_value), pch=20,main="Volcano plot"), col='grey')
sig=subset(res, p_value<.00001 & abs(log2.fold_change.)>1)
with(sig, points(log2.fold_change., -log10(p_value), pch=20, col="red"))
abline(h=5,lty=2)
abline(v=1,lty=2)
abline(v=-1,lty=2)
#write.table(sig,'pvalue001_logfc2.txt',quote=F,sep='\t',col.names = T,row.names = F)

gene=read.table('Jingchuan/gene_name.txt',header=T,sep='\t',stringsAsFactors = F)
for (i in 1:nrow(sig))
{
  sig[i,15]=gene[which(gene[,'id']==sig[i,'gene_id']),'name']
}
colnames(sig)[15]='name'

with(sig,text(log2.fold_change., -log10(p_value)+0.2, labels=name, cex=0.8))

write.table(sig,'Jingchuan/N16N8/n8_sig.txt',quote=F,sep='\t',col.names = T,row.names=F)


#### calculate the significance for genes enriched at telomeres (20kb from telomeres)
#### using fisher's exact test
len=read.table('chr_length_sacCer2.txt',sep='\t',stringsAsFactors = F)
sig=read.table('n2_sig.txt',sep='\t',stringsAsFactors = F,header=T)
out=matrix(NA,nrow=1,ncol=ncol(sig))
colnames(out)=colnames(sig)
out2=matrix(NA,nrow=1,ncol=ncol(sig))
colnames(out2)=colnames(sig)
n=0
for (i in 1:nrow(sig))
{
  chr=strsplit(sig[i,'locus'],':')[[1]][1]
  l=len[which(len[,1]==chr),2]
  if (as.numeric(strsplit(sig[i,'locus'],'-')[[1]][2]) <=20000 | as.numeric(strsplit(sig[i,'locus'],'-')[[1]][2]) >=(l-20000))
  {
    n=n+1
    out=rbind(out,sig[i,])
  }
  else
  {
    out2=rbind(out2,sig[i,])
  }
}

allgenes=read.table('BY4741_gene_name_noCHRM.txt',sep='\t',stringsAsFactors = F)
n2_sig=read.table('n8_sig.txt',sep='\t',stringsAsFactors = F,header=T)[,2]
n2_sig=as.data.frame(n2_sig)
colnames(n2_sig)='gene_id'
colnames(allgenes)=c('chr','sgd','class','start','end','one','strand','two','gene_id','gene_name')
allgenes=merge(n2_sig,allgenes,by.x='gene_id',by.y='gene_id')
n=0
for (i in 1:nrow(allgenes))
{
  chr=allgenes[i,2]
  l=len[which(len[,1]==chr),2]
  if (as.numeric(allgenes[i,6]) <=20000 | as.numeric(allgenes[i,6]) >=(l-20000))
  {
    n=n+1
  }
}

n2=read.table('n8gene_exp_clean.diff',sep='\t',stringsAsFactors = F,header=T)[,2]
n2=as.data.frame(n2)
colnames(n2)='gene_id'
allgenes=merge(n2,allgenes,by.x='gene_id',by.y='gene_id')



# n=2 total genes 6572
test=matrix(c(15,24,213,6359),nrow=2,byrow = F,dimnames =list(c("Subtelomeric_genes", "Non-subtelomeric_genes"),c("Significantly_changed_in_n=2", "All_genes")))
test
fisher.test(test, alternative = "greater")

# n=4 total genes 6586
test=matrix(c(16,0,227,6359),nrow=2,byrow = F,dimnames =list(c("Subtelomeric_genes", "Non-subtelomeric_genes"),c("Significantly_changed_in_n=4", "All_genes")))
test
fisher.test(test, alternative = "greater")

# n=8 total genes 6609
test=matrix(c(9,0,250,6359),nrow=2,byrow = F,dimnames =list(c("Subtelomeric_genes", "Non-subtelomeric_genes"),c("Significantly_changed_in_n=8", "All_genes")))
test
fisher.test(test, alternative = "greater")



############################### n=2 subtelomere 213 genes
test=matrix(c(9,1,5,198),nrow=2,byrow = T,dimnames =list(c("Significant_in_Ellahi_paper", "not_significant_in_Ellahi_paper"),c("Significantly_changed_in_n=2", "not_significantly_changed_in_n=2")))
test
fisher.test(test)

# n=4 subtelomere 227 genes
test=matrix(c(6,2,10,209),nrow=2,byrow = T,dimnames =list(c("Significant_in_Ellahi_paper", "not_significant_in_Ellahi_paper"),c("Significantly_changed_in_n=4", "not_significantly_changed_in_n=4")))
test
fisher.test(test)

# n=8 subtelomere 250 genes
test=matrix(c(5,3,4,238),nrow=2,byrow = T,dimnames =list(c("Significant_in_Ellahi_paper", "not_significant_in_Ellahi_paper"),c("Significantly_changed_in_n=8", "not_significantly_changed_in_n=8")))
test
fisher.test(test)



