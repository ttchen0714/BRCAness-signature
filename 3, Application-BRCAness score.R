#####################    Input
rev_pairs=as.matrix(read.table("Result/BRCAness_signature.txt",header=TRUE,sep="\t")); 
gene1_list<-as.numeric(rev_pairs[,1]);
gene2_list<-as.numeric(rev_pairs[,2]);
#####################    ICGC-OV
exp=read.table("Data/expProfie_ICGC.txt",header=TRUE,row.names=1,sep="\t");
ind1<-match(gene1_list,rownames(exp));
ind2<-match(gene2_list,rownames(exp));
naomit_num=length(na.omit(ind1+ind2));naomit_num
myrank<-as.matrix(unlist(apply(exp,2,rank)))
g1g2_pairnum<-unlist(apply(myrank,2,function(x) length(which(x[ind1]-x[ind2]>0))))
sample_class<-1-g1g2_pairnum/naomit_num;
sample_class<-cbind("g1g2_pairnum",t(sample_class));
colnames(sample_class)<-c("sample",colnames(exp));
#####################    Output
write.table(sample_class,'Result/ICGC-score.txt',sep='\t',row.names=FALSE,quote=FALSE);
