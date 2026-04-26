
#_______0. Source Functions__________________________________________________ 
source("C:/Users/nelly/Desktop/University of Glasgow/R/Further OMICS/functions.r")
output_dir = output_dir = "C:/Users/nelly/Desktop/University of Glasgow/R/Further OMICS/P53_activation_AML/MDE"

library(eulerr)
library(org.Hs.eg.db)
library(ggplot2)

#set group colors - for metagene plot colors
group_colors <- c("DMSO"  = "#B0BEC5",  
                  "BETI"  = "#42A5F5", 
                  "MDM2I" = "#EF5350",  
                  "Combo" = "#AB47BC")  


genes_of_interest = c("CDKN1A", "MDM2", "BBC3", "PMAIP1", "BAX", "NOXA")

#______1. Load Files______________________________________________________________
setwd("C:/Users/nelly/Desktop/University of Glasgow/R/Further OMICS/P53_activation_AML")
em = read.table("em.csv", header= TRUE, row.names=1,sep="\t")
annotations = read.table("anno.csv", header= TRUE, row.names=1,sep="\t")
ss = read.table("ss.csv", header= TRUE, row.names=1,sep="\t")
de_m_v_d = read.table("DE_MDM2I_vs_DMSO.csv", header= TRUE, row.names=1,sep="\t")
de_c_v_d = read.table("DE_Combo_vs_DMSO.csv", header= TRUE, row.names=1,sep="\t")
de_b_v_d = read.table("DE_BETI_vs_DMSO.csv", header= TRUE, row.names=1,sep="\t")


#_____2. Data Wrangling_________________________________________________________

#Add significance column to all de tables 
de_m_v_d$sig = as.factor(de_m_v_d$p.adj < 0.05 & abs(de_m_v_d$log2fold) > 1)
de_c_v_d$sig = as.factor(de_c_v_d$p.adj < 0.05 & abs(de_c_v_d$log2fold) > 1)
de_b_v_d$sig = as.factor(de_b_v_d$p.adj < 0.05 & abs(de_b_v_d$log2fold) > 1)

#Merge the data (two first then the third as merge takes only 2 at a time)
master_temp = merge(de_m_v_d,de_c_v_d, by=0, suffixes=c(".m_v_d",".c_v_d"))
master      = merge(master_temp, de_b_v_d, by.x="Row.names", by.y=0, suffixes=c("", ".b_v_d"))

#Insert the suffixes to the columns from the last merge without the suffix
colnames(master)[colnames(master) %in% c("log2fold", "p", "p.adj", "sig")] = 
  c("log2fold.b_v_d", "p.b_v_d", "p.adj.b_v_d", "sig.b_v_d")

#Merge expression matrix and annotations
master = merge(em, master, by.x=0, by.y=1)
master = merge(annotations, master, by.x=0, by.y=1)
row.names(master) = master$SYMBOL

# extract chromosome from locus column into a new column and retain number only
master$chr = gsub(":.*", "", master$LOCUS)
master$chr = gsub("chr", "", master$chr)

#Get em symbols and em_symbols scaled
em_symbols = master[,row.names(ss)]
em_symbols.s = na.omit(data.frame(t(scale(t(em_symbols)))))

# get significant genes for each comparison
sig_m_v_d = row.names(subset(master, sig.m_v_d == TRUE))
sig_c_v_d = row.names(subset(master, sig.c_v_d == TRUE))
sig_b_v_d = row.names(subset(master, sig.b_v_d == TRUE))

# set the levels for the ss
ss$SAMPLE_GROUP = factor(ss$SAMPLE_GROUP, levels=c("DMSO","BETI","MDM2I","Combo"))

#______3. Correlation heatmap______________________________________________________

# calculate spearman correlation between all samples
cor_matrix = cor(em_symbols, method="spearman")

# melt the correlation matrix into long format for ggplot
cor_matrix.m = melt(cor_matrix)

# plot correlation heatmap
ggp = ggplot(cor_matrix.m, aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(low="#1565C0", mid="#FFFFFF", high="#C62828", midpoint=0.9)+
  labs(title="Sample Correlation Heatmap")+
  savanna_theme+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title   = element_blank())
ggp
save_emf(ggp, "correlation_heatmap.emf")

#_____4. Sig-Any Heat-map_______________________________________________________

# get genes significant in any group 
sig_any = row.names(subset(master, sig.m_v_d == TRUE | sig.c_v_d == TRUE | sig.b_v_d == TRUE))

#make heat map
ggp = make_heatmap(em_symbols.s, sig_any)+
  theme(axis.text.y = element_blank())
ggp
save_emf(ggp, "SigAny_Heatmap_plot.emf")


#_____5. Venn___________________________________________________________________
# prep data
venn_data = list("M_v_D" = sig_m_v_d, "B_v_D" = sig_b_v_d,"C_v_D" = sig_c_v_d)

# plot Venn
ggp = plot(euler(venn_data, shape = "ellipse"),fills = c("#FFD600","lightgreen","#FF6F00"),
     edges = TRUE,labels= list(fontsize=11),quantities = TRUE)
ggp
save_emf(ggp, "venn_plot.emf")


#_____6. Fold vs Fold___________________________________________________________
# MDM2I vs BETI
ggp = ggplot(master, aes(x=log2fold.m_v_d, y=log2fold.b_v_d)) + geom_point()+savanna_theme+labs(title = "MDM2I Vs BETI")
cor.test(master$log2fold.m_v_d, master$log2fold.b_v_d)
save_emf(ggp, "M_V_B_fold.emf")

# MDM2I vs Combo
ggp = ggplot(master, aes(x=log2fold.m_v_d, y=log2fold.c_v_d)) + geom_point()+savanna_theme+labs(title = "MDM2I Vs COMBO") 
cor.test(master$log2fold.m_v_d, master$log2fold.c_v_d)
save_emf(ggp, "M_v_C_fold.emf")

# BETI vs Combo
ggp = ggplot(master, aes(x=log2fold.b_v_d, y=log2fold.c_v_d)) + geom_point()+savanna_theme+labs(title = "BETI Vs COMBO") 
cor.test(master$log2fold.b_v_d, master$log2fold.c_v_d)
save_emf(ggp, "B_v_C.emf")

#_____7. Signatures by K means clustering____________________________________________________________
# scale the data
#em.s = na.omit(data.frame(t(scale(t(em)))))

# do the k-means
km = kmeans(em_symbols.s, 3, iter.max = 10, nstart = 1)
clusters = data.frame(km$cluster)

# get the genes in the clusters
cluster_1_genes = row.names(subset(clusters, clusters == 1))
cluster_2_genes = row.names(subset(clusters, clusters == 2))
cluster_3_genes = row.names(subset(clusters, clusters == 3))


#_______8. Heatmaps of the signatures___________________________________________
ggp = make_heatmap(em_symbols,cluster_1_genes[1:70])+labs(title = "Signature 1")
ggp
save_emf(ggp, "signature1_heatmap.emf")

ggp = make_heatmap(em_symbols,cluster_2_genes[1:70])+labs(title = "Signature 2")
ggp
save_emf(ggp, "signature2_heatmap.emf")

ggp = make_heatmap(em_symbols,cluster_3_genes[1:70])+labs(title = "Signature 3")
ggp
save_emf(ggp, "signature3_heatmap.emf")


#______9. Metagene Analysis_____________________________________________________
#Cluster1
cluster1_em = em_symbols.s[cluster_1_genes,]
ggp=make_metagene(em_symbols.s,cluster_1_genes,cluster1_em,ss$SAMPLE_GROUP)+savanna_theme+labs(title = "Cluster 1")
ggp
save_emf(ggp, "Metagene1.emf")

#Cluster2
cluster2_em = em_symbols.s[cluster_2_genes,]
ggp=make_metagene(em_symbols.s,cluster_2_genes,cluster2_em,ss$SAMPLE_GROUP)+savanna_theme+labs(title = "Cluster 2")
ggp
save_emf(ggp, "Metagene2.emf")

#Cluster3 
cluster3_em = em_symbols.s[cluster_3_genes,]
ggp=make_metagene(em_symbols.s,cluster_3_genes,cluster3_em,ss$SAMPLE_GROUP)+savanna_theme+labs(title = "Cluster 3")
ggp
save_emf(ggp, "Metagene3.emf")

#______10. Pathway analysis of the signatures____________________________________
#cluster1

pathway1_results = do_pathway(org.Hs.eg.db,cluster_1_genes,em_symbols.s,ss$SAMPLE_GROUP,"SYMBOL")

# Extract bar plot
ggp = pathway1_results$plots$ggp.barplot +
  labs(title = "Cluster 1", x = "Gene Count")+ savanna_theme
save_emf(ggp, "pathway1_barplot.emf")

#Extract dot plot
ggp = pathway1_results$plots$ggp.dotplot +
  labs(title = "Cluster 1")+ savanna_theme
save_emf(ggp, "pathway1_dotplot.emf")
#______________________________________________________
#Cluster2
pathway2_results = do_pathway(org.Hs.eg.db,cluster_2_genes,em_symbols.s,ss$SAMPLE_GROUP,"SYMBOL")

# Extract bar plot
ggp = pathway2_results$plots$ggp.barplot +
  labs(title = "Cluster 2", x = "Gene Count")+ savanna_theme
save_emf(ggp, "pathway2_barplot.emf")

#Extract dot plot
ggp = pathway2_results$plots$ggp.dotplot +
  labs(title = "Cluster 2")+ savanna_theme
save_emf(ggp, "pathway2_dotplot.emf")

#_______________________________________________________________________________
#Cluster3
pathway3_results= do_pathway(org.Hs.eg.db,cluster_3_genes,em_symbols.s,ss$SAMPLE_GROUP,"SYMBOL")

# Extract bar plot
ggp = pathway3_results$plots$ggp.barplot +
  labs(title = "Cluster 3", x = "Gene Count")+ savanna_theme
save_emf(ggp, "pathway3_barplot.emf")

#Extract dot plot
ggp = pathway3_results$plots$ggp.dotplot +
  labs(title = "Cluster 3")+ savanna_theme
save_emf(ggp, "pathway3_dotplot.emf")








