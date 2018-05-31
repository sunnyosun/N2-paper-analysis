# processing vcf files after variant calling
# Jingchuan fusion chromosome paper
# 20170711

# SNPs
n16_snp=read.table('BS00557A_S35_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)
n15_snp=read.table('BS00259A_S3_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)
n12_snp=read.table('BS00261A_S5_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)
n8_snp=read.table('BS00262A_S6_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)
n4_snp=read.table('BS00263A_S7_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)
n2_snp=read.table('BS00264A_S8_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)
wt_snp=read.table('BY4741_TGGAACAA_L007_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)

n13_snp=read.table('BS00260A_S4_filtered_snps_final.vcf',skip=50,stringsAsFactors = F)


# compare to wt BY4741
data1=n13_snp
colnames(data1)[1:2]=c('chr','pos')
data2=n16_snp
colnames(data2)[1:2]=c('chr','pos')
chr=unique(n13_snp[,1])
diff=data.frame(matrix(0),nrow=1)
colnames(diff)=c('chr','pos')
for (i in 1:17)
{
  c1=data1[which(data1[,1]==chr[i]),1:2]
  c2=data2[which(data2[,1]==chr[i]),1:2]
  if (length(setdiff(c1[,2],c2[,2])) != 0)
  {
    tmp=data.frame(matrix(chr[i],nrow=length(setdiff(c1[,2],c2[,2]))))
    tmp=cbind(tmp,setdiff(c1[,2],c2[,2]))
    colnames(tmp)=c('chr','pos')
    diff=rbind(diff,tmp)
  }
}
diff=diff[2:nrow(diff),]
out=merge(diff,data1)
# DP>=10, V6>1000, AF=1
out=out[which(out$V6>=800),]
out1=data.frame(matrix(0,ncol=10))
colnames(out1)=colnames(out)
for (i in 1:nrow(out))
{
  line8=out[i,8]
  line10=out[i,10]
  if (strsplit(line10,':')[[1]][3]>10 && strsplit(line8,';')[[1]][2]=="AF=1.00")
  {
    out1=rbind(out1,out[i,])
  }
}
out1=out1[2:nrow(out1),]
write.table(out1,'n13_wt_snp.txt',col.names=T,row.names=F,quote=F,sep='\t')


# INDELs
n16_indel=read.table('BS00557A_S35_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
n15_indel=read.table('BS00259A_S3_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
n12_indel=read.table('BS00261A_S5_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
n8_indel=read.table('BS00262A_S6_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
n4_indel=read.table('BS00263A_S7_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
n2_indel=read.table('BS00264A_S8_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
wt_indel=read.table('BY4741_TGGAACAA_L007_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)
n13_indel=read.table('BS00260A_S4_filtered_indels_recal.vcf',skip=50,stringsAsFactors = F)

data1=n13_indel
colnames(data1)[1:2]=c('chr','pos')
data2=n16_indel
colnames(data2)[1:2]=c('chr','pos')
chr=unique(n13_indel[,1])
diff=data.frame(matrix(0),nrow=1)
colnames(diff)=c('chr','pos')
for (i in 1:17)
{
  c1=data1[which(data1[,1]==chr[i]),1:2]
  c2=data2[which(data2[,1]==chr[i]),1:2]
  if (length(setdiff(c1[,2],c2[,2])) != 0)
  {
    tmp=data.frame(matrix(chr[i],nrow=length(setdiff(c1[,2],c2[,2]))))
    tmp=cbind(tmp,setdiff(c1[,2],c2[,2]))
    colnames(tmp)=c('chr','pos')
    diff=rbind(diff,tmp)
  }
}
diff=diff[2:nrow(diff),]
out=merge(diff,data1)
# DP>=10, V6>1000, AF=1
out=out[which(out$V6>=800),]
out1=data.frame(matrix(0,ncol=10))
colnames(out1)=colnames(out)
for (i in 1:nrow(out))
{
  line8=out[i,8]
  line10=out[i,10]
  if (strsplit(line10,':')[[1]][3]>10 && strsplit(line8,';')[[1]][2]=="AF=1.00")
  {
    out1=rbind(out1,out[i,])
  }
}
out1=out1[2:nrow(out1),]
write.table(out1,'n13_wt_indel.txt',col.names=T,row.names=F,quote=F,sep='\t')


######
# compare to n15
n15_snp=read.table('n15_wt_snp.txt',stringsAsFactors = F,header=T)
n12_snp=read.table('n12_wt_snp.txt',stringsAsFactors = F,header=T)
n8_snp=read.table('n8_wt_snp.txt',stringsAsFactors = F,header=T)
n4_snp=read.table('n4_wt_snp.txt',stringsAsFactors = F,header=T)
n2_snp=read.table('n2_wt_snp.txt',stringsAsFactors = F,header=T)
n13_snp=read.table('n13_wt_snp.txt',stringsAsFactors = F,header=T)


n15_indel=read.table('n15_wt_indel.txt',stringsAsFactors = F,header=T)
n12_indel=read.table('n12_wt_indel.txt',stringsAsFactors = F,header=T)
n8_indel=read.table('n8_wt_indel.txt',stringsAsFactors = F,header=T)
n4_indel=read.table('n4_wt_indel.txt',stringsAsFactors = F,header=T)
n2_indel=read.table('n2_wt_indel.txt',stringsAsFactors = F,header=T)
n13_indel=read.table('n13_wt_indel.txt',stringsAsFactors = F,header=T)


data1=n13_snp
data2=n15_snp
chr=unique(c(data1[,1],data2[,1]))
diff=data.frame(matrix(0),nrow=1)
colnames(diff)=c('chr','pos')
for (i in 1:17)
{
  c1=data1[which(data1[,1]==chr[i]),1:2]
  c2=data2[which(data2[,1]==chr[i]),1:2]
  if (length(setdiff(c1[,2],c2[,2])) != 0)
  {
    tmp=data.frame(matrix(chr[i],nrow=length(setdiff(c1[,2],c2[,2]))))
    tmp=cbind(tmp,setdiff(c1[,2],c2[,2]))
    colnames(tmp)=c('chr','pos')
    diff=rbind(diff,tmp)
  }
}
diff=diff[2:nrow(diff),]
out=merge(diff,data1)
write.table(out,'n13_n15_snp.txt',col.names=T,row.names=F,quote=F,sep='\t')


###
# filtering out the fused part
n12_snp=read.table('n12_n15_snp.txt',stringsAsFactors = F,header=T)
n8_snp=read.table('n8_n15_snp.txt',stringsAsFactors = F,header=T)
n4_snp=read.table('n4_n15_snp.txt',stringsAsFactors = F,header=T)
n2_snp=read.table('n2_n15_snp.txt',stringsAsFactors = F,header=T)
n13_snp=read.table('n13_n15_snp.txt',stringsAsFactors = F,header=T)

n12_indel=read.table('n12_n15_indel.txt',stringsAsFactors = F,header=T)
n8_indel=read.table('n8_n15_indel.txt',stringsAsFactors = F,header=T)
n4_indel=read.table('n4_n15_indel.txt',stringsAsFactors = F,header=T)
n2_indel=read.table('n2_n15_indel.txt',stringsAsFactors = F,header=T)
n13_indel=read.table('n13_n15_indel.txt',stringsAsFactors = F,header=T)

library(xlsx)
info=read.xlsx('~/Desktop/Jingchuan/variant_calling/n=4_remaining coordinate.xlsx', 1, header=T)
data=n4_indel
out=data.frame(matrix(0,ncol=10))
colnames(out)=colnames(data)
for (i in 1:nrow(data))
{
  chr=data[i,1]
  line=info[which(info[,1]==chr),]
  start=line[1,2]
  end=line[1,3]
  cen1=line[1,4]
  cen2=line[1,5]
  if (chr=='chrM')
  {
    out=rbind(out,data[i,])
  }
  else if (is.na(cen1)==T)
  {
    if (data[i,2]>start && data[i,2]<end)
    {
      out=rbind(out,data[i,])
    }
  }
  else if (is.na(cen1)==F)
  {
    if (((data[i,2]>start && data[i,2]<cen1) | (data[i,2]>cen2 && data[i,2]<end))==T)
    {
      out=rbind(out,data[i,])
    }
  }
}
out=out[2:nrow(out),]

out1=data.frame(matrix(0,ncol=10))
colnames(out1)=colnames(data)
for (i in 1:nrow(out))
{
  chr=out[i,1]
  if (chr=='chrM')
  {
    out1=rbind(out1,out[i,])
  }
  else
  {
    line=info[which(info[,1]==chr),]
    snp=strsplit(as.character(as.matrix(line[6])),split = '; ')[[1]]
    s=vector()
    for (j in 1:length(snp))
    {
      tmp=strsplit(as.character(snp[j]),split = '; ')[[1]]
      s=c(s,strsplit(as.character(tmp),split = '[ ]')[[1]][1])
    }
    if ((out[i,2] %in% s)==F)
    {
      out1=rbind(out1,out[i,])
    }
  }
}
out1=out1[2:nrow(out1),]
write.table(out1,'n4_n15_indel_filtered.txt',col.names=T,row.names=F,quote=F,sep='\t')

