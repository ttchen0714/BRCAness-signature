#####################   Input
reverse<-as.matrix(read.table("Result/reverse_pairs_99%fisher0.05_sym.txt",header=TRUE,sep="\t",quote="")) 
reverse[,1:2]<-as.numeric(reverse[,1:2])
g1=reverse[,1];
g2=reverse[,2];
#####################   Input
TCGA_class<-as.matrix(read.table('Data/brca-mutation-TCGA_platin.txt',header=TRUE,row.names=1,sep="\t",quote=""));
colnames(TCGA_class)<-gsub("\\.","-",colnames(TCGA_class));
BRCAmut<-colnames(TCGA_class)[which(TCGA_class[3,]==1)];
BRCAwild<-colnames(TCGA_class)[which(TCGA_class[3,]==0)];
TCGA_trainset=as.matrix(read.table("Data/expProfie_TCGA.txt",header=TRUE,row.names=1,sep="\t",quote=""))
colnames(TCGA_trainset)<-gsub("\\.","-",colnames(TCGA_trainset));
myrank<-as.matrix(unlist(apply(TCGA_trainset,2,rank))) 
gene1_list<-reverse[,1];
gene2_list<-reverse[,2];
#####################   REO
result=matrix(0,ncol=dim(myrank)[2],nrow=dim(reverse)[1]);
diff=myrank[gene1_list,]-myrank[gene2_list,];
result[which(diff>0)]=1
pair<-paste(reverse[,1],reverse[,2],sep=">");
rownames(result)=pair;
class_all=cbind(pair,reverse[,1:2],result);
####################
BRCAmut<-colnames(TCGA_class)[which(TCGA_class[3,]==1)];
BRCAwild<-colnames(TCGA_class)[which(TCGA_class[3,]==0)];#
vary<-result;
colnames(vary)<-colnames(myrank);
#####################   Coverage max
sample_com=intersect(colnames(vary),BRCAwild);
vary=vary[,sample_com];
label=TCGA_class[3,sample_com]
#####################   Greedy Algorithm
mode(vary)="numeric";
coverage_all=unlist(apply(vary,1,function(x) table(x)[1]/length(x)));
sorted_coverage=sort(coverage_all,T)
current_coverage= as.numeric(sorted_coverage[1])  
current_pair=names(sorted_coverage[1])    
current_class =matrix(vary[names(sorted_coverage[1]),],nrow=1)
rownames(current_class)=names(sorted_coverage[1]);   
toAdd_class=vary[row.names(vary)!=names(sorted_coverage[1]),]

for(i in 1:dim(toAdd_class)[1]){
	class_matrix=t(t(toAdd_class)+as.vector(current_class));
	class_matrix[which(class_matrix=="1")]=0;
	class_matrix[which(class_matrix=="2")]=1;
	print(i)
	coverage=unlist(apply(class_matrix,1,function(x) table(x)[1]/length(x)));
	ind_max=which(coverage==max(coverage),arr.ind=TRUE)[1];
	if(coverage[ind_max]>(current_coverage+0.01)){
		current_coverage=as.numeric(coverage[ind_max]);
		current_pair=append(current_pair,rownames(toAdd_class)[ind_max]);
		current_class=class_matrix[ind_max];
		toAdd_class=toAdd_class[-ind_max,];
      	print("add one pair")
	}else break
}
#####################   Output
pair<-paste(reverse[,1],reverse[,2],sep=">");
reverse_gene=reverse[match(current_pair,pair),];
write.table(reverse_gene,'Result/BRCAness_signature.txt',sep='\t',row.names=FALSE,quote=FALSE);









