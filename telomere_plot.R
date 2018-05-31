### plot gene expression vs distance to telomeres

#############
setwd("/Volumes/My Book for Mac/JL_20161221")
n2=read.table('Jingchuan/N16N2/n2gene_exp_clean.diff',sep='\t',header=T,stringsAsFactors = F)
name=read.table('Jingchuan/BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')
# filter out chrM and 2micron
for (i in 1:nrow(n2))
{
  chr=name[which(name[,9]==n2[i,1]),1]
  start=name[which(name[,9]==n2[i,1]),4]
  end=name[which(name[,9]==n2[i,1]),5]
  mid=mean(start,end)
  gene=name[which(name[,9]==n2[i,1]),10]
  n2[i,15]=gene
  n2[i,16]=chr
  n2[i,17]=mid
}
#remove na
n2=n2[!is.na(n2[,6]),]
#change 0 to 0.1
n2[which(n2[,'log2FoldChange']=='-Inf'),'log2FoldChange']=log2(0.1/n2[which(n2[,'log2FoldChange']=='-Inf'),'baseMeanA'])
n2[which(n2[,'log2FoldChange']=='Inf'),'log2FoldChange']=log2(n2[which(n2[,'log2FoldChange']=='Inf'),'baseMeanB']/0.1)

info=read.table('Jingchuan/n2_info.txt',header=T,sep='\t',stringsAsFactors = F)
for (i in 1:16)
{
  chr=info[i,1]
  setEPS()
  postscript(paste(chr,'_telomere_plot.eps'),width=8, height=5)
  x=n2[which(n2[,10]==chr),11]
  y=n2[which(n2[,10]==chr),'log2FoldChange']
  plot(x,y,pch=16,col='dodgerblue',xlab=chr,ylab='fold change',ylim=c(min(y)-0.3,max(y)+0.3))
  abline(h=0)
  dev.off()
}

## only plot 100kb each side
for (i in 1:16)
{
  chr=info[i,1]
  setEPS()
  postscript(paste(chr,'_telomere_plot_100kb_each_side.eps'),width=8, height=5)
  tmp=n2[which(n2[,10]==chr),]
  tmp2=tmp[(tmp[,11]<=100000 | tmp[,11]>=(max(tmp[,11]-100000))),]
  a=which(tmp2[,11]>=(max(tmp2[,11]-100000)))
  tmp2[a,11]=tmp2[a,11]-min(tmp2[a,11])+100000
  x=tmp2[,11]
  y=tmp2[,'log2FoldChange']
  plot(x,y,pch=16,col='dodgerblue',xlab=chr,ylab='log2 fold change',ylim=c(min(y)-0.3,max(y)+0.3),xlim=c(0,220000))
  abline(h=0)
  abline(v=100000,col='red',lty=3,lwd=3)
  dev.off()
}

### metachromosome

name=read.table('BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')
info=read.table('n8_info.txt',header=T,sep='\t',stringsAsFactors = F)

setEPS()
postscript('n8n16_telomere_plot_metachromosome.eps',width=6, height=5)
i=1
chr=info[i,1]
tmp=n2[which(n2[,10]==chr),]
len=max(name[which(name[,1]==chr),4:5])+300
mid =len/2
tmp[which(tmp[,11]>mid),11]=len-tmp[which(tmp[,11]>mid),11]
plot(tmp[,11]/1000, tmp[,'log2FoldChange'], ylim=c(-10,10),xlim=c(0,800),pch=16,,col='dodgerblue',xlab='metachromosome (kb)', ylab='log2FoldChange')
for(i in 2:16)
{
  chr=info[i,1]
  tmp=n2[which(n2[,10]==chr),]
  len=max(name[which(name[,1]==chr),4:5])+300
  mid =len/2
  tmp[which(tmp[,11]>mid),11]=len-tmp[which(tmp[,11]>mid),11]
  points(tmp[,11]/1000, tmp[,'log2FoldChange'],pch=16,,col='dodgerblue')
}
dev.off()


########
#
# read clean data
name=read.table('Jingchuan/BY4741_gene_name.txt',stringsAsFactors = F,sep='\t')

cleanfile='Jingchuan/N16N2/n2gene_exp_clean.diff'
res=read.table(cleanfile,header=T,sep='\t',stringsAsFactors = F)
# remove 2micron and chrM genes
genes_2m=subset(name,V1=='2micron')[,10]
res2=subset(res,!(res$gene %in% genes_2m))
res=res2
genes_mit=subset(name,V1=='chrM')[,10]
res2=subset(res,!(res$gene %in% genes_mit))
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

# filter out chrM and 2micron
n2=res
for (i in 1:nrow(n2))
{
  chr=name[which(name[,9]==n2[i,1]),1]
  start=name[which(name[,9]==n2[i,1]),4]
  end=name[which(name[,9]==n2[i,1]),5]
  mid=mean(start,end)
  gene=name[which(name[,9]==n2[i,1]),10]
  n2[i,15]=gene
  n2[i,16]=chr
  n2[i,17]=mid
}

telo=read.table('Jingchuan/n2_telo.txt',fill=T,sep='\t',stringsAsFactors = F)
telo_fused=data.frame(matrix(0,nrow=1,ncol=17))
colnames(telo_fused)=colnames(n2)
telo_remaining=data.frame(matrix(0,nrow=1,ncol=17))
colnames(telo_remaining)=colnames(n2)

for (i in 1:nrow(telo))
{
  tmp=n2[which(n2$V16==telo[i,1]),]
  if (telo[i,2]=='left')
  {
    telo_remaining=rbind(telo_remaining,tmp[which(tmp$V17<=100000),])
  }
  if (telo[i,2]=='right')
  {
    telo_remaining=rbind(telo_remaining,tmp[which(tmp$V17>=(max(tmp$V17)-100000)),])
  }
}
telo_remaining=telo_remaining[2:nrow(telo_remaining),]

telo_fused=data.frame(setdiff(n2[,1],telo_remaining[,1]))
colnames(telo_fused)='test_id'
telo_fused_data=merge(telo_fused,n2)


info=read.table('Jingchuan/n2_info.txt',header=T,sep='\t',stringsAsFactors = F)
for (i in 1:16)
{
  chr=info[i,1]
  setEPS()
  postscript(paste(chr,'_telomere_plot.eps'),width=8, height=5)
  x=n2[which(n2[,16]==chr),17]
  y=n2[which(n2[,16]==chr),'log2.fold_change.']
  plot(x/1000,y,pch=16,col='dodgerblue',xlab=chr,ylab='fold change',ylim=c(-5,5))
  abline(h=0)
  dev.off()
}


## only plot 100kb each side
for (i in 1:16)
{
  chr=info[i,1]
  setEPS()
  postscript(paste(chr,'_telomere_plot_100kb_each_side.eps'),width=8, height=5)
  tmp=n2[which(n2[,16]==chr),]
  tmp2=tmp[(tmp[,17]<=100000 | tmp[,17]>=(max(tmp[,17]-100000))),]
  a=which(tmp2[,17]>=(max(tmp2[,17]-100000)))
  tmp2[a,17]=tmp2[a,17]-min(tmp2[a,17])+100000
  x=tmp2[,17]
  y=tmp2[,'log2.fold_change.']
  plot(x/1000,y,pch=16,col='dodgerblue',xlab=chr,ylab='log2 fold change',ylim=c(-5,5),xlim=c(0,220))
  abline(h=0)
  abline(v=100,col='red',lty=3,lwd=3)
  dev.off()
}


### meta
n2=telo_fused_data
n2=telo_remaining

setEPS()
postscript('n2_remaining_telomere_plot_metachromosome.eps',width=6, height=5)
i=1
chr=info[i,1]
tmp=n2[which(n2[,16]==chr),]
len=max(name[which(name[,1]==chr),4:5])+300
mid =len/2
tmp[which(tmp[,17]>mid),17]=len-tmp[which(tmp[,17]>mid),17]
plot(tmp[,17]/1000, tmp[,'log2.fold_change.'], ylim=c(-10,10),xlim=c(0,100),pch=16,,col='dodgerblue',xlab='metachromosome (kb)', ylab='log2FoldChange')
for(i in 2:16)
{
  chr=info[i,1]
  tmp=n2[which(n2[,16]==chr),]
  len=max(name[which(name[,1]==chr),4:5])+300
  mid =len/2
  tmp[which(tmp[,17]>mid),17]=len-tmp[which(tmp[,17]>mid),17]
  points(tmp[,17]/1000, tmp[,'log2.fold_change.'],pch=16,,col='dodgerblue')
}
dev.off()

#################
# extract boundary genes
cleanfile='Jingchuan/N16N8/n8gene_exp_clean.diff'
res=read.table(cleanfile,header=T,sep='\t',stringsAsFactors = F)
info=read.table('Jingchuan/n8_info.txt',header=T,sep='\t',stringsAsFactors = F)


out=as.data.frame(matrix(NA,nrow=32,ncol=8))
colnames(out)=c('gene','chr','start','end','n=16','n=2','log2.fold_change.','p_value')
for (i in 1:nrow(info))
{
  out[(2*i-1),2]=info[i,'chr']
  out[(2*i-1),1]=info[i,'start_gene']
  out[(2*i-1),3]=info[i,'start1']
  out[(2*i-1),4]=info[i,'end1']
  out[(2*i-1),5]=res[which(res[,'gene']==info[i,'start_gene']),'value_1']
  out[(2*i-1),6]=res[which(res[,'gene']==info[i,'start_gene']),'value_2']
  out[(2*i-1),7]=res[which(res[,'gene']==info[i,'start_gene']),'log2.fold_change.']
  out[(2*i-1),8]=res[which(res[,'gene']==info[i,'start_gene']),'p_value']
  out[i*2,2]=info[i,'chr']
  out[i*2,1]=info[i,'end_gene']
  out[i*2,3]=info[i,'start2']
  out[i*2,4]=info[i,'end2']
  out[i*2,5]=res[which(res[,'gene']==info[i,'end_gene']),'value_1']
  out[i*2,6]=res[which(res[,'gene']==info[i,'end_gene']),'value_2']
  out[i*2,7]=res[which(res[,'gene']==info[i,'end_gene']),'log2.fold_change.']
  out[i*2,8]=res[which(res[,'gene']==info[i,'end_gene']),'p_value']
}
write.table(out,'n=8_boundary_genes.txt',sep='\t',quote=F,col.names=T,row.names=F)



##############
# calculating percentage of deleted length
setwd("/Volumes/My Book for Mac/JL_20161221")
cleanfile='Jingchuan/N16N2/n2gene_exp_clean.diff'
res=read.table(cleanfile,header=T,sep='\t',stringsAsFactors = F)

file='Jingchuan/N16N2/gene_exp.diff'
res2=read.table(file,header=T,sep='\t',stringsAsFactors = F)

d=setdiff(res2[,'gene'],res[,'gene'])
out=as.data.frame(matrix(NA,ncol=ncol(res2),nrow=length(d)))
colnames(out)=colnames(res2)
for (i in 1:length(d))
{
  out[i,]=res2[which(res2[,'gene']==d[i]),]
}

write.table(out,'n=2_deleted_genes.txt',quote=F,sep='\t',col.names=T,row.names=F)


####
################## verified ORFs
v=read.table('sacCer_verified_ORFs.txt',stringsAsFactors = F,sep='\t',quote = '\"',fill = T)
setwd("/Volumes/My Book for Mac/JL_20161221")
n=read.table('Jingchuan/n=8_deleted_genes.txt',stringsAsFactors = F,header = T)
x=intersect(n[,'gene_id'],v[,2])
out=as.data.frame(matrix(NA,nrow=length(x),ncol=ncol(v)))
for (i in 1:length(x))
{
  out[i,]=v[which(v[,2]==x[i]),]
}
write.table(out,'n=8_deleted_genes_verified_ORF.txt',quote=F,sep='\t',col.names = F,row.names = F)







