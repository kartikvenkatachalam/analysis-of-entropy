setwd("~/") # establishes working directory
library(tidyverse)
library(readr)
library(combinat)
library(varhandle)


disease <- c("BLCA") # defines disease based on TCGA abbreviations
gene_examined <- c("ENSG00000090674") # defines gene of interest, default is MCOLN1


# following block downloads mutations in chosen disease from Xenahubs, and removes all synonymous mutations
url <- paste0("https://gdc.xenahubs.net/download/TCGA-",disease,".mutect2_snv.tsv.gz")
tmp <- tempfile()
download.file(url,tmp)
all_mutations <- read.delim(gzfile(tmp),stringsAsFactors=FALSE)
all_mutations <- filter (all_mutations, all_mutations$effect!="synonymous_variant")
tumors <- as.data.frame(unique(all_mutations$Sample_ID))
colnames(tumors) = "tumors"
mutation_types <- as.data.frame(unique(all_mutations$effect))
all_mutations[,c(9)]<-as.numeric(1)


# following block generates a matrix with CGC list
consensus <-  read_csv("Census_allFri Aug 30 15_11_07 2019.csv") # can be downloaded from https://cancer.sanger.ac.uk/census; registration will be required; default version was from 08/2019
consensus2 <- as.data.frame(unique(consensus$`Gene Symbol`))
colnames(consensus2)<-"cancer_genes"
mut_matrix <- as.data.frame(matrix(0,ncol=nrow(consensus2),nrow=nrow(tumors)))
rownames(mut_matrix) <- tumors$tumors
colnames(mut_matrix) <- consensus2$cancer_genes
mut_matrix<-as.data.frame(cbind(tumors$tumors,mut_matrix))
mut_matrix <- unfactor(mut_matrix)
for (i in 1:nrow(consensus2)){
  mut <- as.character(consensus2[i,1])
  pull1a <- which(all_mutations$gene==mut)
  pull1 <- as.data.frame(all_mutations[c(pull1a),])
  pull1[,c(3)] <- duplicated(pull1[,c(1)])
  pull2 <- filter(pull1, pull1$chrom=="TRUE")
  pull3 <- as.data.frame(unique(pull2$Sample_ID))
  if (nrow(pull3)!=0){
    for ( j in 1:nrow(pull3)){
      pull4 <- as.character(pull3[j,1])
      pull1[c(which(pull1$Sample_ID==pull4)),9] <- as.numeric(nrow(pull3))
    }
    pull1 <- filter(pull1,pull1$chrom!="TRUE")
  }
  for (k in 1:nrow(pull1)){
    pull8 <- as.character(pull1[c(k),1])
    pull10 <- which(mut_matrix[,c(1)] == pull8)
    mut_matrix[c(pull10),c(i+1)]=1
  }
}
mut_matrix[,c(1)]<-NULL


# following block downloads gene expression FPKM values in chosen disease from Xenahubs and median-normalizes the expression of the chosen gene (default, MCOLN1)
url <- paste0("https://gdc.xenahubs.net/download/TCGA-",disease,".htseq_fpkm.tsv.gz")
tmp <- tempfile()
download.file(url,tmp)
fpkm <- read.delim(gzfile(tmp), stringsAsFactors=FALSE)
gene_exam <- as.data.frame(colnames(fpkm))
gene_exam[,c(1)]<-gsub("[.]","-",gene_exam[,c(1)])
gene_examb <- as.data.frame(t(filter(fpkm,grepl(gene_examined,fpkm$Ensembl_ID))))
gene_exam <- as.data.frame(cbind(gene_exam,gene_examb))
gene_exam <- gene_exam[-c(1),]
gene_exam <- unfactor(gene_exam)
tumors2 <- tumors
for (k in 1:nrow(tumors)){
  tumor <- as.character(tumors[c(k),1])
  val<- which(gene_exam[,c(1)]==tumor)
  if (is_empty(val)==F){
    tumors2[c(k),2] <- as.numeric(gene_exam[c(val),2])
  } else {tumors2[c(k),2] <-NA}
}
colnames(tumors2)<-c("gene_exam_expression","gene_exam_expression")
mut_matrix <- as.data.frame(cbind(tumors2,mut_matrix))
mut_matrix[,c(1)]<-NULL

# following block removes tumors in which gene expression returned NA
gene_exam_rem1<-as.data.frame(mut_matrix$gene_exam_expression)
rem_list <- which(is.na(gene_exam_rem1$`mut_matrix$gene_exam_expression`==T))
if (is_empty(rem_list)==FALSE){
  mut_matrix2 <- as.data.frame(mut_matrix[-c(rem_list),])
} else{mut_matrix2<-mut_matrix}

#following block removes genes that do not exhibit any mutations in the chosen disease
count_matrix <- mut_matrix2
remove1 <- which(colSums(count_matrix)==0)
count_matrix2 <- as.data.frame(count_matrix[,-c(remove1)])



# following block calculates entropy 
count_matrix2 <- count_matrix2 %>% arrange(desc(gene_exam_expression))
n1=round(nrow(count_matrix2)*0.33) # determines the split in the dataset, default is top and bottom third
n2<-nrow(count_matrix2)
cluster1 <- as.data.frame(count_matrix2[c(1:n1),])
cluster1b <- as.data.frame(count_matrix2[-c(1:n1),])
cluster2 <- as.data.frame(count_matrix2[c((n2-n1+1):n2),])
cluster2b <- as.data.frame(count_matrix2[-c((n2-n1+1):n2),])

#cluster1-highest expressors
entropy <- as.data.frame(colnames(count_matrix2))
entropy <- as.data.frame(entropy[-c(1:3),])
colnames(entropy) <- c("genes")
entropy <- entropy %>% mutate(cluster1_information=NA,cluster1_p_value=NA,cluster2_information=NA,cluster2_p_value=NA)
for(j in 1:nrow(entropy)){
  gene1_1<-as.character(entropy[j,1])
  gene2_1<-which(colnames(cluster1)==gene1_1)
  gene3_1<- as.data.frame(cluster1[,gene2_1])
  gene4_1<- as.data.frame(cluster1b[,gene2_1])
  gene10_1<- as.data.frame(count_matrix2[,gene2_1])
  P_A <- as.double(nrow(cluster1)/nrow(count_matrix2))
  P_BA <- as.double(colSums(gene3_1)/nrow(gene3_1))
  P_C <- as.double(colSums(gene4_1)/nrow(gene4_1))
  P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
  prob1_1 = as.double(P_A*P_BA/P_B)
  prob1_1b <- 1-(prob1_1)
  if (prob1_1==0 || prob1_1==1 || prob1_1b<0){
    entrop1_1=0
  } else {entrop1_1 <- as.double((-prob1_1)*(log2(prob1_1))+(-prob1_1b)*(log2(prob1_1b)))}
  entrop1=entrop1_1
  P_A<-1-P_A
  P_BA <- as.double(colSums(gene4_1)/nrow(gene4_1))
  P_C <- as.double(colSums(gene3_1)/nrow(gene3_1))
  P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
  prob2_1 <- as.double((P_A*P_BA)/P_B)
  prob2_1b =1-prob2_1
  if (prob2_1==0 || prob2_1==1|| prob2_1b<0){
    entrop2_1=0
  } else {entrop2_1 <- as.double((-prob2_1)*(log2(prob2_1))+(-prob2_1b)*(log2(prob2_1b)))}
  entrop2=entrop2_1
  prob3_1 <- as.double(colSums(gene10_1)/nrow(gene10_1))
  prob3_1b <- 1-prob3_1
  if (prob3_1==0 || prob3_1==1){
    entrop3_1=0
  } else {entrop3_1 <- as.double((-prob3_1)*(log2(prob3_1))+(-prob3_1b)*(log2(prob3_1b)))}
  entrop3=entrop3_1
  entrop_final<- as.double(entrop3-(((nrow(gene3_1)/nrow(gene10_1))*entrop1)+((nrow(gene4_1)/nrow(gene10_1))*entrop2)))
  entropy[j,2]=entrop_final
}

#cluster2-lowest expressors
for(j in 1:nrow(entropy)){
  gene1_1<-as.character(entropy[j,1])
  gene2_1<-which(colnames(cluster2)==gene1_1)
  gene3_1<- as.data.frame(cluster2[,gene2_1])
  gene4_1<- as.data.frame(cluster2b[,gene2_1])
  gene10_1<- as.data.frame(count_matrix2[,gene2_1])
  P_A <- as.double(nrow(cluster2)/nrow(count_matrix2))
  P_BA <- as.double(colSums(gene3_1)/nrow(gene3_1))
  P_C <- as.double(colSums(gene4_1)/nrow(gene4_1))
  P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
  prob1_1 = as.double(P_A*P_BA/P_B)
  prob1_1b <- 1-(prob1_1)
  if (prob1_1==0 || prob1_1==1 || prob1_1b<0){
    entrop1_1=0
  } else {entrop1_1 <- as.double((-prob1_1)*(log2(prob1_1))+(-prob1_1b)*(log2(prob1_1b)))}
  entrop1=entrop1_1
  P_A<-1-P_A
  P_BA <- as.double(colSums(gene4_1)/nrow(gene4_1))
  P_C <- as.double(colSums(gene3_1)/nrow(gene3_1))
  P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
  prob2_1 <- as.double((P_A*P_BA)/P_B)
  prob2_1b =1-prob2_1
  if (prob2_1==0 || prob2_1==1|| prob2_1b<0){
    entrop2_1=0
  } else {entrop2_1 <- as.double((-prob2_1)*(log2(prob2_1))+(-prob2_1b)*(log2(prob2_1b)))}
  entrop2=entrop2_1
  prob3_1 <- as.double(colSums(gene10_1)/nrow(gene10_1))
  prob3_1b <- 1-prob3_1
  if (prob3_1==0 || prob3_1==1){
    entrop3_1=0
  } else {entrop3_1 <- as.double((-prob3_1)*(log2(prob3_1))+(-prob3_1b)*(log2(prob3_1b)))}
  entrop3=entrop3_1
  entrop_final<- as.double(entrop3-(((nrow(gene3_1)/nrow(gene10_1))*entrop1)+((nrow(gene4_1)/nrow(gene10_1))*entrop2)))
  entropy[j,4]=entrop_final
}

#following block performes permuation test, default number of permutations is 1000
permutations = 1000 # can be changed
entropy_random <- as.data.frame(colnames(count_matrix2))
entropy_random <- as.data.frame(entropy_random[-c(1:3),])
colnames(entropy_random) <- c("genes")
entropy_random2 <- entropy_random
for (x in 1:permutations){
  matrix_random1 <- as.data.frame(count_matrix2[,c(1:3)])
  matrix_random2 <- as.data.frame(count_matrix2[,-c(1:3)])
  matrix_random1 <- matrix_random1[sample(nrow(matrix_random1)),]
  matrix_random <- as.data.frame(cbind(matrix_random1,matrix_random2))
  matrix_random <-  matrix_random %>% arrange(desc(gene_exam_expression))
  cluster1 <- as.data.frame(matrix_random[c(1:n1),])
  cluster1b <- as.data.frame(matrix_random[-c(1:n1),])
  cluster2 <- as.data.frame(matrix_random[c((n2-n1+1):n2),])
  cluster2b <- as.data.frame(matrix_random[-c((n2-n1+1):n2),])
  for(j in 1:nrow(entropy)){
    gene1_1<-as.character(entropy[j,1])
    gene2_1<-which(colnames(cluster1)==gene1_1)
    gene3_1<- as.data.frame(cluster1[,gene2_1])
    gene4_1<- as.data.frame(cluster1b[,gene2_1])
    gene10_1<- as.data.frame(count_matrix2[,gene2_1])
    P_A <- as.double(nrow(cluster1)/nrow(count_matrix2))
    P_BA <- as.double(colSums(gene3_1)/nrow(gene3_1))
    P_C <- as.double(colSums(gene4_1)/nrow(gene4_1))
    P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
    prob1_1 = as.double(P_A*P_BA/P_B)
    prob1_1b <- 1-(prob1_1)
    if (prob1_1==0 || prob1_1==1 || prob1_1b<0){
      entrop1_1=0
    } else {entrop1_1 <- as.double((-prob1_1)*(log2(prob1_1))+(-prob1_1b)*(log2(prob1_1b)))}
    entrop1=entrop1_1
    P_A<-1-P_A
    P_BA <- as.double(colSums(gene4_1)/nrow(gene4_1))
    P_C <- as.double(colSums(gene3_1)/nrow(gene3_1))
    P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
    prob2_1 <- as.double((P_A*P_BA)/P_B)
    prob2_1b =1-prob2_1
    if (prob2_1==0 || prob2_1==1|| prob2_1b<0){
      entrop2_1=0
    } else {entrop2_1 <- as.double((-prob2_1)*(log2(prob2_1))+(-prob2_1b)*(log2(prob2_1b)))}
    entrop2=entrop2_1
    prob3_1 <-as.double(colSums(gene10_1)/nrow(gene10_1))
    prob3_1b <- 1-prob3_1
    if (prob3_1==0 || prob3_1==1){
      entrop3_1=0
    } else {entrop3_1 <- as.double((-prob3_1)*(log2(prob3_1))+(-prob3_1b)*(log2(prob3_1b)))}
    entrop3=entrop3_1
    entrop_final<- as.double(entrop3-(((nrow(gene3_1)/nrow(gene10_1))*entrop1)+((nrow(gene4_1)/nrow(gene10_1))*entrop2)))
    entropy_random[j,(x+1)]=entrop_final
    entropy_random <- as.data.frame(entropy_random)
  }
  for(j in 1:nrow(entropy)){
    gene1_1<-as.character(entropy[j,1])
    gene2_1<-which(colnames(cluster2)==gene1_1)
    gene3_1<- as.data.frame(cluster2[,gene2_1])
    gene4_1<- as.data.frame(cluster2b[,gene2_1])
    gene10_1<- as.data.frame(count_matrix2[,gene2_1])
    P_A <- as.double(nrow(cluster2)/nrow(count_matrix2))
    P_BA <- as.double(colSums(gene3_1)/nrow(gene3_1))
    P_C <- as.double(colSums(gene4_1)/nrow(gene4_1))
    P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
    prob1_1 = as.double(P_A*P_BA/P_B)
    prob1_1b <- 1-(prob1_1)
    if (prob1_1==0 || prob1_1==1 || prob1_1b<0){
      entrop1_1=0
    } else {entrop1_1 <- as.double((-prob1_1)*(log2(prob1_1))+(-prob1_1b)*(log2(prob1_1b)))}
    entrop1=entrop1_1
    P_A<-1-P_A
    P_BA <- as.double(colSums(gene4_1)/nrow(gene4_1))
    P_C <- as.double(colSums(gene3_1)/nrow(gene3_1))
    P_B <- as.double((P_A*P_BA)+((1-P_A)*(P_C)))
    prob2_1 <- as.double((P_A*P_BA)/P_B)
    prob2_1b =1-prob2_1
    if (prob2_1==0 || prob2_1==1|| prob2_1b<0){
      entrop2_1=0
    } else {entrop2_1 <- as.double((-prob2_1)*(log2(prob2_1))+(-prob2_1b)*(log2(prob2_1b)))}
    entrop2=entrop2_1
    prob3_1 <-as.double(colSums(gene10_1)/nrow(gene10_1))
    prob3_1b <- 1-prob3_1
    if (prob3_1==0 || prob3_1==1){
      entrop3_1=0
    } else {entrop3_1 <- as.double((-prob3_1)*(log2(prob3_1))+(-prob3_1b)*(log2(prob3_1b)))}
    entrop3=entrop3_1
    entrop_final<- as.double(entrop3-(((nrow(gene3_1)/nrow(gene10_1))*entrop1)+((nrow(gene4_1)/nrow(gene10_1))*entrop2)))
    entropy_random2[j,(x+1)]=entrop_final
    entropy_random2 <- as.data.frame(entropy_random2)
  }
}

#following block provides nominal P-values
for (y in 1:nrow(entropy)){
  a1 <- entropy[y,2]
  b1<-entropy_random[y,]
  b1<-b1[,-c(1)]
  b1<- t(b1)
  test1 <- mean(b1>a1)
  if (var(b1)!=0){
    test2 <- t.test(b1,mu=a1,alternative = "two.sided")
    entropy[y,3]=as.numeric(test2$p.value)
  } else{entropy[y,3]=NA}
  a2 <- entropy[y,4]
  b2<-entropy_random2[y,]
  b2<-b2[,-c(1)]
  b2<- t(b2)
  test1b <- mean(b2>a2)
  if (var(b2)!=0){
    test2b <- t.test(b2,mu=a2,alternative = "two.sided")
    entropy[y,5]=as.numeric(test2b$p.value)
  } else{entropy[y,5]=NA}
}

#following block determines FDR and identified genes with FDR < 0.05 and information >0
names1 <- as.data.frame(entropy$genes)
names2 <- as.data.frame(rbind(names1,names1))
vals1 <- as.data.frame(entropy[,c(2:3)])
vals2 <- as.data.frame(entropy[,c(4:5)])
colnames(vals1)<-colnames(vals2) 
vals3 <- as.data.frame(rbind(vals1,vals2))
entropy2 <- as.data.frame(cbind(names2,vals3))
entropy2 <- entropy2 %>% mutate(p_adj = p.adjust(entropy2$cluster2_p_value,method = "fdr"))
entropy3 <- filter(entropy2, entropy2$p_adj<0.05)
colnames(entropy3) <- c("genes","information","nominal_p-val","FDR")
entropy3 <- entropy3 %>% arrange(desc(information))
entropy3 <- filter(entropy3,entropy3$information>0)

# following block removes duplicates
remove_dups <- as.numeric(which(duplicated(entropy3$genes)==T))
if (is_empty(remove_dups)==F){
  entropy3 <- as.data.frame(entropy3[-c(remove_dups),])
} else {entropy3 <- entropy3}


# following block determines the number of tumors with mutations in the genes
for (i in 1:nrow(entropy3)){
  gene <- as.character(entropy3[c(i),1])
  gene_column <- as.data.frame(count_matrix2[,c(which(colnames(count_matrix2)==gene))])
  entropy3[i,5] <- as.numeric(colSums(gene_column))
}
colnames(entropy3) <- c("genes","information","nominal_p-val","FDR","number of tumors with mutations")
entropy3 <- entropy3 %>% arrange (desc(`number of tumors with mutations`))
write_csv(entropy3,paste0(gene_examined," entropy analysis in ",disease,"_final results.csv")) #results are saved to working directory
entropy4 <- filter(entropy3,!entropy3$`number of tumors with mutations`<15)
write_csv(entropy4,paste0(gene_examined," entropy analysis in ",disease,"_final results, >15 mutations.csv")) #results are saved to working directory
