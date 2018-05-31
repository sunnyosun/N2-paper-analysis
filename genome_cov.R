# plotting the genome coverage by chromosomes
# Jingchuan's fusion chromosome genome sequencing data

# load in the genomecov data
name='wt'
data=read.table('BS00259A_S3_genomecov.bedgraph',sep='\t',stringsAsFactors = F)
chrs_sort=c('chrI','chrVI','chrIII','chrIX','chrVIII','chrV','chrXI','chrX','chrXIV','chrII','chrXIII','chrXVI','chrXII','chrVII','chrXV','chrIV')
d=list()
for (i in 1:16)
{
  d[[i]]=data[which(data[,1]==chrs_sort[i]),]
}
rm(data)

d_sm=list()
for (n in 1:16)
{
  tmp=d[[n]]
  sm=as.data.frame(matrix(0,ncol=3,nrow=(floor(nrow(tmp)/100)-1)))
  for (i in 1:(floor(nrow(tmp)/100)-1))
  {
    sm[i,]=colMeans(tmp[(100*(i-1)+1):(100*i),2:4])
  }
  sm[,4]=rowMeans(sm[,1:2])
  d_sm[[n]]=sm[,c(4,3)]
}
d=d_sm


setEPS()
postscript("wt_genomecov.eps",height=30,width=5,paper='special')
par(mfrow=c(16,1),lab=c(x=4,y=3,len=1),las=1,tcl=-0.25,mgp=c(0.3,0.2,0),
    bty='n',mar=c(1.2,1.6,1,1),cex.lab=0.7,cex.axis=0.7,font.lab=1,cex=1,
    omi=c(0.3,0.3,0.3,0.3))
for (n in 1:16)
{
  plot(d[[n]][,1]/1000,d[[n]][,2],pch=16,cex=0.7,xlab='',xlim=c(0,1531),ylim=c(0,400),ylab='',xaxt='n',main=chrs_sort[n],col='blue')
  t=floor(max(d[[n]][,1]/1000)/100)
  axis(1,c(100*(0:t),max(d[[n]][,1]/1000)),c(100*(0:t),max(d[[n]][,1]/1000)))
}
dev.off()


