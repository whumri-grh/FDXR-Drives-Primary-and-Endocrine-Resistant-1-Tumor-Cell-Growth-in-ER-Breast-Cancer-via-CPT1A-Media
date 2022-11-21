setwd("dir")
getwd()
rm(list = ls())
#Import data and filter
mycounts <- read.csv("FDXRshRNA-Metabolomics from Xiaojing Liu.csv",header=T, row.names = 1)
mycounts1 <- subset(mycounts,mycounts$II_434_1 > 1000)
mycounts2 <- subset(mycounts1,mycounts1$II_434_2 > 1000)
mycounts3 <- subset(mycounts2,mycounts2$II_434_3 > 1000)
mycounts4 <- na.omit(mycounts3)
rm(mycounts1,mycounts2,mycounts3)
group <- rep(c('II_Scr', 'II_434'), each = 3)
#EdgeR analysis of differentially expressed genes
library(edgeR)
dgelist <- DGEList(counts = mycounts4, group = group)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~rev(group))
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
#Export Analysis Results
resdata <- as.data.frame(lrt)
resdata$id <- rownames(resdata)
write.table(resdata, file="FDXRshRNA.edgeR_output.csv",quote = F,col.names = T,row.names = F,sep = ",")
#Differential gene sorting and threshold setting
resdata <- resdata[order(resdata$FDR, resdata$logFC, decreasing = c(FALSE, TRUE)), ]
resdata$threshold = factor(ifelse(resdata$PValue < 0.05 & resdata$logFC >= log(1.2,2),'Sig.Up',ifelse(resdata$PValue < 0.05 & resdata$logFC<=log(0.8,2) ,'Sig.Down','Unsig')),levels=c('Sig.Up','Sig.Down','Unsig'))
summary(resdata)
#Draw heat map
mycounts$id <- rownames(mycounts)
A1 <- subset(resdata,resdata$threshold=="Sig.Up")
A1 <- merge(A1,mycounts,by="id")
A2 <- subset(resdata,resdata$threshold=="Sig.Down")
A2 <- merge(A2,mycounts,by="id")
A3 <- rbind(A1,A2)
rownames(A3) <-A3$id
B <- A3[,c(8:13)]
library(ComplexHeatmap)
for (i in 1:nrow(B)) B[i, ] <- scale(log(unlist(B[i, ] + 1), 2))
B <- as.matrix(B)
Heatmap(B,
        name ="Group",
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        clustering_distance_rows = "manhattan",
        border_gp = gpar(col="black"),
        #row_km = 2,
        #column_km = 2,
        #row_gap = unit(4, "mm"),
        #column_gap = unit(4,"mm"), 
        cluster_columns = F,
        cluster_rows = F,
        show_row_names = F,
        show_column_names = F,
        row_title =NULL,
        column_title = NULL,
        heatmap_width = unit(4,"cm"),
        heatmap_height = unit(8,"cm"),
        show_heatmap_legend = T)
#Limma analysis microarray data
mycounts <- read.csv("FDXR shRNA microarray for heatmap_1.csv",header=T, row.names = 1)
mycounts1 <- mycounts[rowSums(mycounts > 1) >= 3, ]
group <- factor(c(rep("Scr",2),rep("FDXR_sh434",2)), levels = c("Scr","FDXR_sh434"))
library(edgeR)
library(limma)
design <- model.matrix(~group)
rownames(design)=colnames(mycounts1)
dge <- DGEList(counts = mycounts1)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)
fit2=eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom = na.omit(tempOutput)
DEG_voom <- DEG_voom[order(DEG_voom$P.Value, DEG_voom$logFC, decreasing = c(TRUE, FALSE)), ]
DEG_voom <- DEG_voom[DEG_voom$P.Value < 0.05 & abs(DEG_voom$logFC) > log(1.05,2),]
DEG_voom$threshold = factor(ifelse(DEG_voom$P.Value < 0.05 & abs(DEG_voom$logFC) >= log(1.07,2), ifelse(DEG_voom$logFC>= log(1.07,2) ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
summary(DEG_voom)
A1 <- subset(DEG_voom,DEG_voom$threshold=="Up")
A2 <- subset(DEG_voom,DEG_voom$threshold=="Down")
A3 <- rbind(A1,A2)
A3$ID <- row.names(A3)
gene_trans<-read.table("gene_trans.txt",fileEncoding = "UTF16",row.names = NULL,header = TRUE)
colnames(gene_trans)<-c("ID","gene_name")
DEG_voom <- merge(A1,gene_trans,by = "ID")
mycounts1$ID <- rownames(mycounts1)
DEG <- merge(DEG_voom,mycounts1,by="ID")
#Export FDXR Downstream Regulated Genes
A2$ID <- row.names(A2)
FDXR_Downstream <- merge(A2,gene_trans,by = "ID")
mycounts1$ID <- rownames(mycounts1)
FDXR_Downstream <- merge(FDXR_Downstream,mycounts1,by="ID")
write.table(FDXR_Downstream, file="FDXR_Downstream_gene.csv",quote = F,col.names = T,row.names = F,sep = ",")
#Export FDXR Upstream Regulated Genes
A1$ID <- row.names(A1)
FDXR_Upstream <- merge(A1,gene_trans,by = "ID")
mycounts1$ID <- rownames(mycounts1)
FDXR_Upstream <- merge(FDXR_Upstream,mycounts1,by="ID")
write.table(FDXR_Downstream, file="FDXR_Upstream_gene.csv",quote = F,col.names = T,row.names = F,sep = ",")
#Draw heat map
library(ComplexHeatmap)
mycounts2 <- DEG_DOWN[,9:13]
rownames(mycounts2) <- mycounts2$gene_name
mycounts2 <- mycounts2[,-1]
for (i in 1:nrow(mycounts2)) mycounts2[i, ] <- scale(log(unlist(mycounts2[i, ] + 1), 2))
mycounts2 <- as.matrix(mycounts2)
Heatmap(mycounts2,
        name ="Group",
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        clustering_distance_rows = "manhattan",
        border_gp = gpar(col="black"),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = T,
        show_column_names = F,
        row_names_gp = gpar(fontsize=8),
        row_title =NULL,
        column_title = NULL,
        heatmap_width = unit(4,"cm"),
        heatmap_height = unit(8,"cm"),
        show_heatmap_legend = T,
        heatmap_legend_param = list(title_position="topleft"))
#TCGA Data Processing
install.packages("rjson")
library("rjson")
result <- fromJSON(file = "metadata.cart.2022-09-19.json")
metadata <- data.frame(t(sapply(result,function(x){
  id <- x$associated_entities[[1]]$entity_submitter_id
  file_name <- x$file_name
  all <- cbind(id,file_name)
})))
rownames(metadata) <- metadata[,2]
dir <- './all/'
samples <- list.files(dir)
sampledir <- paste0(dir,samples)
mat <- do.call(cbind,lapply(sampledir,function(x){
  rt <- data.table::fread(x,data.table = F)
  rownames(rt) <- rt[,1]
  rt <- rt[,8]
}))
rt <- data.table::fread('./all/00063026-2a67-41de-b105-f80dca277978.rna_seq.augmented_star_gene_counts.tsv',data.table = F)
colnames(mat) <- sapply(strsplit(sampledir,'/'),'[',3)
rownames(mat) <- rt$gene_id
mat1 <- t(mat)
same <- intersect(row.names(metadata),row.names(mat1))
data <- cbind(metadata[same,],mat1[same,])
rownames(data) <- data[,1]
tcga_stad <- t(data)
tcga_stad <- tcga_stad[-c(1:6),]
rownames(rt) <- rt[,1]
same2 <- intersect(row.names(rt),row.names(tcga_stad))
tcga <- cbind(rt[same2,],tcga_stad[same2,])
tcga <- tcga[-c(1,4:9)]
library(limma)
rt <- tcga[,-2]
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data)>0,]
write.table(data,"mRNAmatrix.txt",sep = "\t",quote = F)
write.table(data, file="mRNAmatrix.csv",quote = F,col.names = T,row.names = F,sep = ",")
#Combined TCGA gene expression data with clinical data
library("XML")
library("methods")
dir <- "F:/数据分析结果/R数据分析/单基因单癌分析/clinic_data/"
all_fiels=list.files(path = dir,pattern='*.xml$',recursive=T)
cl = lapply(all_fiels, function(x){
  result <- xmlParse(file = file.path(dir,x)) 
  rootnode <- xmlRoot(result)  
  xmldataframe <- xmlToDataFrame( rootnode[2] ) 
  return(t(xmldataframe)) })
clinical <- t(do.call(cbind,cl))
write.table(clinical, file="clinic.csv",quote = F,col.names = T,row.names = F,sep = ",")
library(limma)
expFile <- "mRNAmatrix.txt"
cliFile <- "time.txt"
rt <- read.table(expFile,header = T,sep = "\t",check.names = F)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
data <- data[rowMeans(data)>0,]
data <- t(data)
rownames(data) <- substr(rownames(data),1,12)
cli <- read.table(cliFile,header = T,sep = "\t",check.names = F,row.names = 1)
cli$status <- gsub("Alive","0",cli$status)
cli$status <- gsub("Dead","1",cli$status)
sameSample <- intersect(row.names(data),row.names(cli))
data <- data[sameSample,]
cli <- cli[sameSample,]
out <- cbind(cli,data)
out <- cbind(id=row.names(out),out)
write.table(out,file = "expTime06.txt",sep = "\t",row.names = F,quote = F)
write.table(rt, file="expTime06.csv",quote = F,col.names = T,row.names = F,sep = ",")
#Plotting survival curves
#TCGA
library(limma)
library(survival)
library(survminer)
gene <- "CPT1A"
expFile <- "mRNAmatrix.txt"
surFile <- "time.txt"
geneexp <- read.table(expFile,header = T,sep = "\t",check.names = F)
geneexp <- as.matrix(geneexp)
rownames(geneexp) <- geneexp[,1]
exp <- geneexp[,2:ncol(geneexp)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
data <- t(data[gene,,drop=F])
group <- sapply(strsplit(rownames(data),"\\-"),"[",4)
group <- sapply(strsplit(group,""),"[",1)
group <- gsub("2","1",group)
data <- data[group==0,,drop=F]
rownames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
data <- avereps(data)
surTime <- read.table(surFile,header = T,sep = "\t",check.names = F,row.names = 1,quote = "")
surTime <- surTime[which(surTime$breast_carcinoma_estrogen_receptor_status=="Positive1"),]
sameSample <- intersect(row.names(data),row.names(surTime))
data <- data[sameSample,,drop=F]
surTime <- surTime[sameSample,,drop=F]
rt <- cbind(surTime,data)
rt$survival_time=rt$survival_time/365
colnames(rt)[ncol(rt)] <- "gene"
rt$status <- gsub("Alive","0",rt$status)
rt$status <- gsub("Dead","1",rt$status)
rt$status <- as.numeric(rt$status)
rt <- as.data.frame(rt)
rt$group <- ifelse(rt$gene>median(rt$gene),"High","Low")
rt$ID <- rownames(rt)
write.table(rt, file="CPT1A-survival.csv",quote = F,col.names = T,row.names = F,sep = ",")
fit <- survfit(Surv(survival_time,status) ~group, data = rt)
diff <- survdiff(Surv(survival_time,status) ~group, data = rt)
pValue <- 1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{pValue=paste0("p=",sprintf("%.03f",pValue))}
ggsurvplot(fit,
           #data=rt,
           conf.int = F,
           pval = pValue,
           pval.size=6,
           legend.title=gene,
           legend.labs=c("High","Low"),
           xlab="Time(years)",
           break.time.by=2,
           xlim=c(0,10),
           palette = c("red","blue"),
           risk.table = F,
           risk.table.title="",
           risk.table.height=.25)
#METABRIC
gene_exp <- read.table("mRNA expression z-scores relative to all samples (log microarray).txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(gene_exp)[2] <- c("ID")
clinical_data <- read.table("brca_metabric_clinical_data.txt",header = T,sep = "\t",check.names = F,quote = "")
clinical_data <- na.omit(clinical_data)
clinical_data <- subset(clinical_data,clinical_data$`ER Status`=="Positive")
clinical_data$`Overall Survival (Months)` <- clinical_data$`Overall Survival (Months)`/12
clinical_data$`Overall Survival Status` <- gsub("0:LIVING","0",clinical_data$`Overall Survival Status`)
clinical_data$`Overall Survival Status` <- gsub("1:DECEASED","1",clinical_data$`Overall Survival Status`)
clinical_data_Cox <- clinical_data[,c(3,26,27)]
colnames(clinical_data_Cox)[1] <- c("ID")
CPT1A_COX <- merge(gene_exp,clinical_data_Cox,by ="ID" )
colnames(CPT1A_COX)[4] <- c("survival_time")
colnames(CPT1A_COX)[5] <- c("status")
CPT1A_COX$CPT1A <- as.numeric(CPT1A_COX$CPT1A)
CPT1A_COX$survival_time <- as.numeric(CPT1A_COX$survival_time)
CPT1A_COX$status <- as.numeric(CPT1A_COX$status)
write.table(CPT1A_COX, file="CPT1A-survival.csv",quote = F,col.names = T,row.names = F,sep = ",")
library(survival)
library(survminer)
cox <- coxph(Surv(survival_time,status)~CPT1A,data = CPT1A_COX)
coxsummary <- summary(cox)
CPT1A_COX$group <- ifelse(CPT1A_COX$CPT1A>0,"High","Low")
sum(CPT1A_COX$group=="High")
sum(CPT1A_COX$group=="Low")
fit <- survfit(Surv(survival_time,status) ~group, data = CPT1A_COX)
ggsurvplot(fit,
           #data=rt,
           conf.int = F,
           pval = pValue,
           pval.size=6,
           legend.title=gene,
           legend.labs=c("High","Low"),
           xlab="Time(years)",
           break.time.by=2,
           xlim=c(0,10),
           palette = c("red","blue"),
           risk.table = F,
           risk.table.title="",
           risk.table.height=.25)
#Gene correlation analysis
rm(list=ls())
library(limma)
expFile <- "mRNAmatrix.txt"
surFile <- "time.txt"
gene <- "CPT1A"
gene1 <- "FDXR"
geneexp <- read.table(expFile,header = T,sep = "\t",check.names = F)
geneexp <- as.matrix(geneexp)
rownames(geneexp) <- geneexp[,1]
exp <- geneexp[,2:ncol(geneexp)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
data <- t(data)
group <- sapply(strsplit(rownames(data),"\\-"),"[",4)
group <- sapply(strsplit(group,""),"[",1)
group <- gsub("2","1",group)
data <- data[group==0,,drop=F]
rownames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
data <- avereps(data)
surTime <- read.table(surFile,header = T,sep = "\t",check.names = F,row.names = 1,quote = "")
surTime <- surTime[which(surTime$breast_carcinoma_estrogen_receptor_status=="Positive1"),]
sameSample <- intersect(row.names(data),row.names(surTime))
data <- data[sameSample,,drop=F]
write.table(data, file="mRNAmatrix-ER-Positive.csv",quote = F,col.names = T,row.names = F,sep = ",")
data1 <- t(data)
data2 <- intersect(row.names(data1),row.names(FDXR))
data3 <- data1[data2,,drop=F]
data4 <- scale(t(data3))
data5 <- as.data.frame(rowMeans(data4))
data6 <- data1[gene,,drop=F]
data7 <- scale(t(data6))
data7 <- as.data.frame(data7)
data9 <- data1[gene1,,drop=F]
data10 <- scale(t(data9))
data10 <- as.data.frame(data10)
test <- cor.test(as.numeric(unlist(data5)),as.numeric(unlist(data10)),type="spearman")
data10$ID <- rownames(data10)
data5$ID <- rownames(data5)
data8 <- merge(data5,data10,by="ID")
colnames(data8) <- c("ID","GENE","FDXR")
rownames(data8) <- data8$ID
data8 <- data8[,-1]
data8$GENE <- as.numeric(unlist(data8$GENE))
data8$FDXR <- as.numeric(unlist(data8$FDXR))
#Correlation Plot
library(ggstatsplot)
ggstatsplot::ggscatterstats(data = as.data.frame(data8), 
                            x = FDXR,
                            y = GENE,
                            ylab = "OXPHOS-related genes",
                            type = "pearson",
                            centrality.para = "mean",                              
                            margins = "both",                                         
                            xfill = "#CC79A7", 
                            yfill = "#009E73", 
                            marginal = F,
                            title = "Association of FDXR with genes involved in oxidative phosphorylation",
                            ggtheme = ggplot2::theme_classic())

