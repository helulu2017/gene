##################################################
rm(list=ls())
library(GEOquery)
library(limma)
gene <- getGEO('GSE43696', destdir=".",getGPL = F)
GPL_GENE=read.csv("GPL6480-9577.csv")
colnames(GPL_GENE)
##################################################
genename= GPL_GENE[,c(1,6)]
exprSet <- as.data.frame(exprs(gene[[1]]))
exprSet$ID <- rownames(exprSet)
express <- merge(x = genename, y = exprSet, by = "ID")
express$ID =NULL
#####################################################
omitnull<-function(cc){
  id=c()
  for(i in 1:nrow(cc))
  {
    
    if((length(which(cc[i,]==""))) >0){
      id=c(id,F) }else{
        id=c(id,T)
      }
  }
  return(cc[id,])
}

exprSet_symbol=(omitnull(express))
dim(exprSet_symbol)
dim(express)
exprSet_symbol <- aggregate(x = exprSet_symbol[,-1],by = list(GENE=exprSet_symbol$GENE_SYMBOL), FUN = max)
dim(exprSet_symbol)  
#################################################
pdata_g = pData(gene[[1]]) 
#################################################
# # exprSet_g = exprs(gene[[1]])
# write.csv(exprSet_g, paste("gene","_exprSet.csv")) 
# write.csv(pdata_g, paste("gene","_metadata.csv")) 
exprSet_g =as.matrix(exprSet_symbol[,-1])
rownames(exprSet_g)=exprSet_symbol$GENE
boxplot(exprSet_g)
#sampleNames(eSet_g)
n.sample=ncol(exprSet_g)
tmp=ifelse(pdata_g$`disease state:ch1`=="Control","Control","asthma")
pdata_g$disease=tmp
##################################################

##################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
# site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# install.packages(c("WGCNA", "stringi", "reshape2"), repos=site)
library(WGCNA)
library(reshape2)
library(stringi)
datTraits=pdata_g
WGCNA_matrix = t(exprSet_g[order(apply(exprSet_g,1,mad), decreasing = T)[1:5000],])
datExpr0 <- WGCNA_matrix  
datExpr <- datExpr0 
#核对行名，可不运行
# sampleNames = rownames(datExpr)
# traitRows = match(sampleNames, datTraits$geo_accession)  
# rownames(datTraits) = datTraits[traitRows, 2]
##################################################
##################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
power=sft$powerEstimate
power
##################################################
#     
##################################################
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)
#write.csv(table(mergedColors),"mergedColors.csv")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
##################################
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datExpr_tree<-hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
##################################################
#   
##################################################
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# design=model.matrix(~0+ datTraits$disease)
# colnames(design)=levels(as.factor(datTraits$disease))
moduleColors <- labels2colors(net$colors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

str(datTraits)
datTraits2=datTraits[,c("disease","disease state:ch1","gender:ch1","age:ch1")]
str(datTraits2)
colnames(datTraits2)=c("Asthma Status","Asthma Severity","Gender","Age")
datTraits2$Gender=ifelse(datTraits2$Gender=="female","0","1")
datTraits2$`Asthma Status`=ifelse(datTraits2$`Asthma Status`=="Control","0","1")
unique(datTraits2$`Asthma Severity`)
datTraits2$`Asthma Severity`=ifelse(datTraits2$`Asthma Severity`=="Control","0",
                                    ifelse(datTraits2$`Asthma Severity`=="Moderate Asthma","1","2"))

moduleTraitCor = cor(MEs, datTraits2, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3,8,0.2,0.5));
length(names(datTraits2))
names(MEs)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               xLabelsAngle = 30,
               cex.lab.x = 0.8,
               cex.lab.y = 0.8)
#main = paste("Module-trait relationships"))
##################################################

##################################################
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
Luminal = as.data.frame(datTraits2[,1]);
names(Luminal) = "Status"
geneTraitSignificance = as.data.frame(cor(datExpr, Luminal, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Luminal), sep="");
names(GSPvalue) = paste("p.GS.", names(Luminal), sep="");
##
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Asthma status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
##################################################
#     
##################################################
geneTree = net$dendrograms[[1]]
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)
plotTOM = dissTOM^7
diag(plotTOM) = NA 
#
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
nSelect = 400#
set.seed(10);
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
####
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
sizeGrWindow(9,9)
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, 
        selectTree, 
        selectColors, 
        main = "Network heatmap plot, selected genes")
##################################################
#   
##################################################
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
TYPE = as.data.frame(datTraits2[,1])
names(TYPE) = "Status"
MET = orderMEs(cbind(MEs, TYPE))
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle= 90)

##################################################
#    
##################################################
datKME=signedKME(datExpr, MEs, outputColumnName="kME_MM.")							
#write.csv(datKME, "kME_MM_test.csv")
##################################################
#    
##################################################
TOM = TOMsimilarityFromExpr(datExpr, power = 4)
module = "purple"
probes = colnames(datExpr) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
modProbes#
#write.csv(modProbes, paste("red-new.csv"))#
##################################################
#   
##################################################
# #
tmp=GPL_GENE[,c(5,6)]
tmp1=data.frame(modProbes,stringsAsFactors = FALSE)
colnames(tmp1)=c("GENE_SYMBOL")
modProbes_purple=merge(x = tmp1, y = tmp, by = "GENE_SYMBOL")
modProbes_purple$GENE_SYMBOL=NULL
modProbes_purple=unique(modProbes_purple)
# #zuotu:
library(clusterProfiler)
ego_cc=enrichGO(OrgDb="org.Hs.eg.db", gene = modProbes_purple$GENE,
                pvalueCutoff = 0.05,readable=TRUE)#
#barplot(ego_cc, showCategory=15,title="EnrichmentGO")#
#dotplot(ego_cc,title="EnrichmentGO")#
write.csv(as.data.frame(ego_cc),"G-enrich_purple.csv",row.names =F)
#go图可视化#
GO_BP <- read.csv("G-enrich_red+purple.csv")
GO_BP$point_shape<-"0"
GO_BP$point_size<-"10"
library(ggplot2)
ggplot(data=GO_BP)+
  geom_bar(aes(x=reorder(Description,Count),y=Count, fill=-log10(pvalue)), stat='identity') +
  coord_flip()+
  scale_fill_gradient(expression(-log["10"]("P.value")),low="red", high = "blue") +
  xlab("") +
  ylab("Gene count")+
  scale_y_continuous(expand=c(0, 0))+
  theme(
    axis.text.x=element_text(color="black",size=rel(1.5)),                                                                                                                                                                                                                                                                                       
    axis.text.y=element_text(color="black", size=rel(1.6)),                                                                                                                                                                                                                                                                                      
    axis.title.x = element_text(color="black", size=rel(1.6)),                                                                                                                                                                                                                                                                                      
    legend.text=element_text(color="black",size=rel(1.0)),
    legend.title = element_text(color="black",size=rel(1.1))
  ) 
dev.off()

##################################################
#  
##################################################
# TOM = TOMsimilarityFromExpr(datExpr, power = 4)
# module = "red"
# probes = colnames(datExpr) 
# inModule = (moduleColors==module)
# modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# cyt = exportNetworkToCytoscape(
#   modTOM,
#   edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
#   nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
#   weighted = TRUE,
#   threshold = 0.02,
#   nodeNames = modProbes, 
#   nodeAttr = moduleColors[inModule]
# )
# #
nTop = 100
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]
##

cyt = exportNetworkToCytoscape(
  modTOM,
  #filter,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02
  #nodeNames = modProbes, 
  #nodeAttr = moduleColors[inModule]
)
