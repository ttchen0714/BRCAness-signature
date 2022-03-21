setwd("ttchen0714/BRCAness-signature/2-gene pairs"); 
dyn.load("FisherExactTest64.dll")
#####################    Input
TCGA_trainset = as.matrix(read.table("Data/expProfie_TCGA.txt",header=TRUE,row.names=1,sep="\t",quote=""))
colnames(TCGA_trainset)<-gsub("\\.","-",colnames(TCGA_trainset));
rowname=rownames(TCGA_trainset)
TCGA_trainset=apply(TCGA_trainset,2,as.numeric)
rownames(TCGA_trainset)=rowname
# duplicated samples after remove other information? 
if(length(unique(substr(colnames(TCGA_trainset),1,12)))==length(substr(colnames(TCGA_trainset),1,12)))
{
    print("No duplicated samples after remove other information in sample names!\n")
}else{
    print("Duplicated samples after remove other information in sample names!\n")
}
#####################    Input
TCGA_class = as.matrix(read.table("Data/brca-mutation-TCGA_platin.txt",header=TRUE,row.names=1,sep="\t",quote=""))
colnames(TCGA_class)<-as.matrix(unlist(apply(as.matrix(colnames(TCGA_class)),1,function(x) substr(x,1,12))))
colnames(TCGA_class)<-gsub("\\.","-",colnames(TCGA_class))
tcga_BRCAmut = TCGA_trainset[,intersect(colnames(TCGA_trainset),names(which(TCGA_class[3,]=="1")))]
tcga_BRCAwild = TCGA_trainset[,intersect(colnames(TCGA_trainset),names(which(TCGA_class[3,]=="0")))]
#####################   Stable gene pairs in mutated samples
cutoff =0.99
stable_matrix_mut = matrix(0,nrow=nrow(tcga_BRCAmut),ncol=nrow(tcga_BRCAmut),dimnames=list(rownames(tcga_BRCAmut),rownames(tcga_BRCAmut)))
# create a progress bar
pb <- txtProgressBar(min = 0, max = nrow(tcga_BRCAmut), style = 3);
for(i in 1:nrow(tcga_BRCAmut))
{ 
    tmp = unlist(apply(matrix(tcga_BRCAmut[i,],nrow=nrow(tcga_BRCAmut),ncol=ncol(tcga_BRCAmut),byrow=T)-tcga_BRCAmut,1,function(x){return(length(which(x>0)))}))
    tmp[!(tmp>=cutoff*ncol(tcga_BRCAmut))] = 0
    stable_matrix_mut[i,] = tmp
    setTxtProgressBar(pb,i)
}
close(pb)
colnames(stable_matrix_mut)<-rowname;
rownames(stable_matrix_mut)<-rowname;
#####################   Reversed gene pairs in wild samples
fisher_pvalue = matrix(NA,nrow=sum(stable_matrix_mut>0),ncol=5)
pb <- txtProgressBar(min = 0, max = nrow(stable_matrix_mut), style = 3)
index = 1
for(i in 1:nrow(stable_matrix_mut))
{
    # skip the gene without stable pairs
    if(all(stable_matrix_mut[i,]==0))
    {
        next
    }else{
    tmp_a = stable_matrix_mut[i,stable_matrix_mut[i,]>0]
    tmp_c = ncol(tcga_BRCAmut)-tmp_a
    tmp_b = unlist(apply(matrix(rep(tcga_BRCAwild[i,],length(which(stable_matrix_mut[i,]>0))),nrow=length(which(stable_matrix_mut[i,]>0)),byrow=TRUE)-tcga_BRCAwild[which(stable_matrix_mut[i,]>0),],1,function(x){length(which(x>0))}))
    tmp_d = ncol(tcga_BRCAwild)-tmp_b
    tmp = c(rbind(tmp_a,tmp_c,tmp_b,tmp_d))
    tmp_pvalue = .C("FisherExactTest",as.integer(tmp),as.integer(length(tmp)),as.numeric(c(0,lfactorial(1:ncol(TCGA_trainset)))),p.value=as.numeric(rep(0,length(tmp)/4)))$p.value
    fisher_pvalue[index:(index+length(tmp_a)-1),] = cbind(rep(rownames(stable_matrix_mut)[i],length(tmp_a)),rownames(stable_matrix_mut)[stable_matrix_mut[i,]>0],tmp_a,tmp_b,tmp_pvalue)
    index = index+length(tmp_a)
    setTxtProgressBar(pb,i)
    }
}
close(pb)
#####################    Significance
result_filter_cutoff = 0.05 
fisher_filter_result=fisher_pvalue[which(as.numeric(fisher_pvalue[,5])<result_filter_cutoff),];
fisher_filter_result=fisher_filter_result[order(as.numeric(fisher_filter_result[,5])),];
colnames(fisher_filter_result) = c("Gene1_ID","Gene2_ID","G1>G2_mut","G1>G2_wild","fisher_pvalue")
#####################    Id to symbol
gene_info = as.matrix(read.table("Data/Homo_sapiens.gene_info.txt",comment.char="",header=F,skip=1,sep="\t",quote=""))[,c(2,3,5)];
gene_info[,1] = as.character(as.numeric(gene_info[,1]))
id1<-fisher_filter_result[,1];
id2<-fisher_filter_result[,2];
sym1<-gene_info[match(id1,gene_info[,1]),2];
sym2<-gene_info[match(id2,gene_info[,1]),2];
#####################    Output
fisher_filter_result<-cbind(fisher_filter_result[,1:2],sym1,sym2,fisher_filter_result[,-c(1,2)]);
colnames(fisher_filter_result) = c("Gene1_ID","Gene2_ID","Gene1_Symbol","Gene2_Symbol","G1>G2_mut","G1>G2_wild","fisher_pvalue")
write.table(fisher_filter_result,"Result/reverse_pairs_99%fisher0.05_sym.txt",col.names=TRUE,row.names=F,sep="\t",quote=F)
