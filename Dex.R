rm(list = ls())
library(tidyverse)
library(data.table)
library(DESeq2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(ggvenn)
library(sva)
library(openxlsx)

setwd("~/Desktop/pav104 R/")
customPlot = list(
  theme(plot.title = element_text(size=14, face="bold", vjust=1, hjust = 0.5, colour="#052049"), 
        axis.title.y = element_text(size=14, face="bold", vjust=1, colour="#052049"),
        axis.title.x = element_text(size=14, face="bold", vjust=0, colour="#052049"),
        axis.text = element_text(size=14, face="plain", colour="black"), 
        axis.text.x = element_text(size=14, face="plain", angle=0, hjust=0.5, vjust=1, colour="black"),
        strip.text.x = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size=14, face="plain"), legend.key.size = unit(4, "mm"),
        legend.title = element_text(size=14, face="plain"),
        legend.spacing = unit(0, "mm"), plot.margin = unit(c(2,0,2,0), "mm"),
        panel.background = element_rect(colour = "#052049", fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA), strip.text = element_text(size=9, face="plain"),
  )
)

genes = fread("Homo_sapiens.GRCh38.99.genes.tsv", header = F) %>% 
  rename("V1" = "chrom", "V2" = "start", "V3" = "end", "V4" = "strand", 
         "V5" = "geneID", "V6" = "symbol", "V7" = "biotype")

countFiles = list.files("counts/")
counts = do.call(cbind, 
        lapply(1:length(countFiles), function(x){
          id = gsub('"', "", strsplit((deparse(countFiles[x])), "\\.")[[1]][1])
          fread(file.path("counts", countFiles[x])) %>% filter(!grepl("_", V1)) %>%
            column_to_rownames(var = "V1") %>% select(V2) %>% rename("V2" = id)
        }))

samples = readxl::read_excel("samples.xlsx") %>% column_to_rownames(var = "sampleID")
all(rownames(samples)==colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ condition)

## 1. remove genes with very low counts
dim(dds)
keep <- rowSums(counts(dds)>=10) >= 5
dds <- dds[keep,]
dim(dds)

dds$condition <- factor(dds$condition, levels = c("uninfected","infected", "pav104"))
dds <- estimateSizeFactors(dds)

design(dds) <- formula(~ donor + condition)
dds = DESeq(dds)
test1 = results(dds, contrast=c("condition","infected","uninfected")) %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  merge(genes, ., by = "geneID") %>% arrange(padj)
test2 = results(dds, contrast=c("condition","pav104","uninfected")) %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>%
  merge(genes, ., by = "geneID") %>% arrange(padj)
test3 = results(dds, contrast=c("condition","pav104","infected")) %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>%
  merge(genes, ., by = "geneID") %>% arrange(padj)
list_of_sheets <- list("infected_uninfected" = test1, 
                       "pav104_uninfected" = test2, 
                       "pav104_infected" = test3)
#write.xlsx(list_of_sheets, file = "results/dex.xlsx")

#write.xlsx(test1, file="results/dex.xlsx", sheetName="infected_uninfected", row.names=FALSE)
#write.xlsx(test2, file="results/dex.xlsx", sheetName="pav104_uninfected", append=TRUE, row.names=FALSE)
#write.xlsx(test3, file="results/dex.xlsx", sheetName="pav104_infected", append=TRUE, row.names=FALSE)

p = ggplot() + 
  geom_point(data = test1 %>% filter(biotype=='protein_coding' & pvalue<=0.05 & log2FoldChange>=0.5),
             aes(x = log2FoldChange, y = -log10(pvalue)), color = "red")+
  geom_point(data = test1 %>% filter(biotype=='protein_coding' & pvalue<=0.05 & log2FoldChange<= -0.5), 
             aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue")+
  geom_point(data = test1 %>% filter(biotype=='protein_coding' & ((pvalue>0.05) | (abs(log2FoldChange)<0.5))), 
             aes(x = log2FoldChange, y = -log10(pvalue)), color = "black")+
  theme_bw()+ggtitle("Infected vs uninfected")+customPlot+
  geom_text(data = test1 %>% filter(biotype=='protein_coding' & (-log10(pvalue)>7|abs(log2FoldChange)>5)),
            aes(x = log2FoldChange, y = -log10(pvalue), label = symbol))
#ggsave("results/infected_uninfected.pdf", plot = p, width=8.27, height=8.27, units="in")
p = ggplot() + 
  geom_point(data = test2 %>% filter(pvalue<=0.05 & log2FoldChange>=0.5), aes(x = log2FoldChange, y = -log10(pvalue)), color = "red")+
  geom_point(data = test2 %>% filter(pvalue<=0.05 & log2FoldChange<= -0.5), aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue")+
  geom_point(data = test2 %>% filter((pvalue>0.05) | (abs(log2FoldChange)<0.5)), aes(x = log2FoldChange, y = -log10(pvalue)), color = "black")+
  theme_bw()+ggtitle("PAV104 vs uninfected")+customPlot+
  geom_text(data = test2 %>% filter(biotype=='protein_coding' & -log10(pvalue)>80),
            aes(x = log2FoldChange, y = -log10(pvalue), label = symbol))
#ggsave("results/pav104_uninfected.pdf", plot = p, width=8.27, height=8.27, units="in")
p = ggplot() + 
  geom_point(data = test3 %>% filter(pvalue<=0.05 & log2FoldChange>=0.5), aes(x = log2FoldChange, y = -log10(pvalue)), color = "red")+
  geom_point(data = test3 %>% filter(pvalue<=0.05 & log2FoldChange<= -0.5), aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue")+
  geom_point(data = test3 %>% filter((pvalue>0.05) | (abs(log2FoldChange)<0.5)), aes(x = log2FoldChange, y = -log10(pvalue)), color = "black")+
  theme_bw()+ggtitle("PAV104 vs infected")+customPlot+
  geom_text(data = test3 %>% filter(biotype=='protein_coding' & -log10(pvalue)>80),
            aes(x = log2FoldChange, y = -log10(pvalue), label = symbol))
#ggsave("results/pav104_infected.pdf", plot = p, width=8.27, height=8.27, units="in")


library("genefilter")
vsd <- vst(dds, blind = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
mat = mat %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>% merge(genes[,c("symbol","geneID")]) %>% 
  select(-geneID) %>% column_to_rownames(var = "symbol") %>% as.matrix()
annotation = data.frame(Condition = factor(c(rep("uninfected", 5), rep("infected", 5), rep("pav104", 5)),
                                            levels = c("uninfected", "infected", "pav104"),
                                            labels = c("Control", "SARS-CoV-2 infection", 
                                                       "SARS-CoV-2 infection + PAV-104")), 
                   Donor = factor(rep(1:5, 3), labels = c("1", "2", "3", "4", "5")),
                   row.names = colnames(mat))
annotation_colors = list(
  Condition = c("Control" = "#c1272d", "SARS-CoV-2 infection" = "#0000a7",
                "SARS-CoV-2 infection + PAV-104" = "#008176"),
  Donor = c("1"="#ffd700", "2"="#ffb14e", "3" = "#fa8775", "4" = "#ea5f94", "5" = "#cd34b5"))

paletteLength <- 1000
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(mat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat)/paletteLength, max(mat), length.out=floor(paletteLength/2)))

#pdf("results/heatmap.pdf", height = 8, width = 10)
pheatmap(mat, annotation = annotation, annotation_colors = annotation_colors, show_colnames = F,
         fontsize_row = 7, 
         #color=colorRampPalette(rev(brewer.pal(n=5,name="RdBu")))(255)
         color=myColor, breaks=myBreaks
         )
#dev.off()


### Show batch effects using sva package
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(counts(dds, normalized=TRUE), mod, mod0, n.sv=2)
qplot(svseq$sv[,1], svseq$sv[,2], color=factor(dds$donor), shape=dds$condition,
      xlab="SV1", ylab="SV2")

ddssva = dds
ddssva$SV1 = svseq$sv[,1]
ddssva$SV2 = svseq$sv[,2]
design(ddssva) = ~ SV1 + SV2 + condition
ddssva = DESeq(ddssva)
ggvenn(list(
  new = results(ddssva, contrast=c("condition","infected","uninfected")) %>% 
    as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
    merge(genes, ., by = "geneID") %>% arrange(padj) %>% 
    filter(padj<0.05) %>% pull(symbol),
  
  old = test1 %>% filter(padj<0.05) %>% pull(symbol)
))
###

####draw volcano with interested genes marked
infile = "~/Downloads/Copy of dex.xlsx"
sheetnames = excel_sheets(infile)
form.names = c("SARS-CoV-2 vs Control", "SARS-CoV-2+PAV-104 vs Control", "SARS-CoV-2+PAV-104 vs SARS-CoV-2")
ps = list()
for (i in 1:length(sheetnames)){
  dat = read_excel("~/Downloads/Copy of dex.xlsx", sheet = sheetnames[i])
  ps[[i]] = ggplot() + 
    geom_point(data = dat %>% filter(biotype=='protein_coding' & padj<=0.05 & log2FoldChange>=0.5),
               aes(x = log2FoldChange, y = -log10(padj)), color = "red", alpha = 0.3)+
    geom_point(data = dat %>% filter(biotype=='protein_coding' & padj<=0.05 & log2FoldChange<= -0.5), 
               aes(x = log2FoldChange, y = -log10(padj)), color = "blue", alpha = 0.3)+
    geom_point(data = dat %>% filter(biotype=='protein_coding' & ((padj>0.05) | (abs(log2FoldChange)<0.5))), 
               aes(x = log2FoldChange, y = -log10(padj)), color = "grey", alpha = 0.3)+
    theme_bw()+ggtitle(form.names[i])+customPlot+
    geom_text_repel(data = dat %>% filter(`show in the plot`=='Y'),
                    aes(x = log2FoldChange, y = -log10(padj), label = symbol), 
                    size = 3, fontface = "bold", max.overlaps = Inf)
  
}

#pdf("results/volcanos.pdf", height = 8, width = 14)
margin = theme(plot.margin = unit(c(0.5,0.5, 0.5, 0.5), "cm"))
grid.arrange(grobs = lapply(ps, "+", margin), ncol = 3)
#dev.off()

#pdf("results/venn.pdf", height = 8, width = 8)
ggvenn(list(
  `SARS-CoV-2 \nvs Control` = list_of_sheets$infected_uninfected %>% filter(padj<=0.05) %>% pull(geneID),
  `SARS-CoV-2+PAV-104 \nvs Control` = list_of_sheets$pav104_uninfected %>% filter(padj<=0.05) %>% pull(geneID),
  `SARS-CoV-2+PAV-104 \nvs SARS-CoV-2` = list_of_sheets$pav104_infected %>% filter(padj<=0.05) %>% pull(geneID)
))
#dev.off()

