###############################################################################
# Differential Expression Analysis Pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Description : A generic, reusable pipeline for RNA-seq data analysis.
#               Covers: data loading, parsing, QC plots, volcano/MA/PCA,
#               heat maps, ORA pathway enrichment.
#
# Outputs : .emf plots saved to output_dir
#
# Author  : Nelly Nyaga
# Date    : 19/03/2026
###############################################################################

#___0. CONFIGURATIONS - to be changed before running the script_________________
setwd("C:/Users/nelly/Desktop/University of Glasgow/R/Further OMICS/P53_activation_AML")

output_dir = "C:/Users/nelly/Desktop/University of Glasgow/R/Further OMICS/P53_activation_AML/Mdm2I_v_dmso"

source("C:/Users/nelly/Desktop/University of Glasgow/R/Further OMICS/functions.r")


# Group labels - must match values in sample_group column of sample sheet
# Add or remove groups as needed
group_1 = "DMSO"
group_2 = "MDM2I"

#set group colors
group_colors <- c("DMSO"  = "#B0BEC5",  
                  "MDM2I" = "#EF5350")  


# colors for volcano and MA plots
direction_colors = c("unchanged"     = "#E0E0E0",
                     "upregulated"   = "#FF6F00",
                     "downregulated" = "#00838F")

#Set number of top genes needed for labeling on volcano/MA or making box plots 
top_genes = 5

#Genes of Interest
genes_of_interest = c("CDKN1A", "MDM2", "BBC3", "PMAIP1", "BAX", "NOXA")


#_____1. Load libraries_________________________________________________________
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(clusterProfiler)  
library(org.Hs.eg.db)
library(STRINGdb)
library(devEMF)



#____2. Load Data_______________________________________________________________
em = read.table("em.csv", header= TRUE, row.names=1,sep="\t")
de = read.table("DE_MDM2I_vs_DMSO.csv", header= TRUE, row.names=1,sep="\t")
annotations = read.table("anno.csv", header= TRUE, row.names=1,sep="\t")
ss = read.table("ss.csv", header= TRUE, row.names=1,sep="\t")

#subset to keep groups
ss = ss[-c(4:6, 10:12),,drop=FALSE ]
em = em[, -c(4:6, 10:12)]


#____3. Parse Data______________________________________________________________
#Build master: expression + de + annotations
master_temp = merge(em, annotations, by.x=0, by.y=0) 
master  = merge(master_temp,de,by.x=1, by.y=0)

#Rename columns 
row.names(master) = master[,"Row.names"]
row.names(master) = master[,"SYMBOL"]
names(master)[1] = "ENSEMBL"

#Create em_symbols
n_samples = ncol(em)  # gets number of expression columns
em_symbols = master[, rownames(ss)]

# Remove rows with NA values from the master table
master = na.omit(master)

# Sort the master table by p.adj value
master = master[order(master[, "p.adj"], decreasing = FALSE), ]

#Add new derived columns to master 
master$means = rowMeans(master[, rownames(ss)])
master$mlog10p = -log10(master$p.adj)
master$sig = as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1.0)#flags significance as true or false 

# Create a scaled expression matrix from em_symbols
em_scaled = data.frame(t(scale(t(em_symbols))))
em_scaled =na.omit(em_scaled)

# Create subsets of significant genes
master_sig      = subset(master, p.adj < 0.05 & abs(log2fold) >1)
master_sig_up   = subset(master,p.adj <0.05 & log2fold > 1 )
master_sig_down = subset(master, p.adj <0.05 & log2fold < -1)
master_non_sig  = subset(master, sig == FALSE)

#Add a new column called direction to master tagging every gene based on its subset
master_non_sig$direction   = "unchanged" 
master_sig_up$direction    = "upregulated"         
master_sig_down$direction  = "downregulated"  
master=rbind(master_non_sig,master_sig_up,master_sig_down)

# extract chromosome from locus column into a new column and retain number only
master$chr = gsub(":.*", "", master$LOCUS)
master$chr = gsub("chr", "", master$chr)

#Count number of up and down genes to add to plots downstream
n_up   = sum(master$direction == "upregulated")
n_down = sum(master$direction == "downregulated")

#Get the list of top N genes that shall be displayed on plots. 
top_up =  master_sig_up[seq_len(top_genes), ]
top_down = master_sig_down[seq_len(top_genes), ]

#Create a vector of significant gene names by row names
sig_genes = row.names(master_sig)

# Create an expression table of significant genes only
em_symbols_sig = em_symbols[sig_genes,]

# Create a scaled expression table of significant genes only
em_scaled_sig = em_scaled[sig_genes,]

#_______4. Plots________________________________________________________________
#_______4.1 Volcano Plot________________________________________________________
ggp = ggplot(master, aes(x=log2fold, y=mlog10p, color=direction))+
  geom_point()+
  labs(title=" MDM2i Volcano Plot", x="log2(Fold Change)", y="-log10(p.adj)")+
  savanna_theme+
  xlim(c(-10,10))+
  ylim(c(0,30))+
  scale_colour_manual(values=direction_colors, name="Legend")+
  annotate("text", x=10,  y=30, label=paste(n_up),   color="black", fontface="bold", size=3, hjust=2)+
  annotate("text", x=-7, y=30, label=paste(n_down), color="black", fontface="bold", size=3, hjust=2)+
  geom_label_repel(data=top_up,   aes(label=SYMBOL), color="black", fill=direction_colors["upregulated"],   label.padding=0.2, label.size=0.25, force=2, max.overlaps=20, show.legend=FALSE)+
  geom_label_repel(data=top_down, aes(label=SYMBOL), color="black", fill=direction_colors["downregulated"], label.padding=0.2, label.size=0.25, force=2, max.overlaps=20, show.legend=FALSE)
ggp
save_emf(ggp, "volcano_plot.emf")


#_______4.2 MA plot_____________________________________________________________
ggp = ggplot(master, aes(x=log10(means), y=log2fold, colour=direction))+
  geom_point()+
  labs(title="MA Plot", x="Mean Expression(log10)", y="Log2 fold change")+
  savanna_theme+
  xlim(c(1,5))+
  ylim(c(-10,+8))+
  scale_colour_manual(values=direction_colors, name="Legend")
ggp
save_emf(ggp, "ma_plot.emf")



#____4.3  PCA plot______________________________________________________________ 
ggp = make_pc1_pc2(colour_groups = ss$SAMPLE_GROUP, e_data = em_symbols)
ggp = ggp + savanna_theme + scale_colour_manual(values = group_colors)
ggp
save_emf(ggp, "pca_plot.emf")


#_____4.4 Expression density __________________________________________________
# derive expanded colors automatically from sample sheet
group_colors_expanded = group_colors[ss$SAMPLE_GROUP]
names(group_colors_expanded) = rownames(ss)

em.m = melt(em)
em.m$sample_group = gsub("_[0-9]+$", "", em.m$variable)

ggp = ggplot(em.m, aes(x=log10(value+0.01)))+
  geom_density(aes(fill=variable), linewidth=0.5, alpha=0.75, show.legend=FALSE)+
  labs(x="Expression(log10)", y="Density")+
  facet_wrap(~variable, ncol=3, scales="free")+
  scale_x_continuous(breaks=seq(-2, 6, by=2))+
  scale_fill_manual(values=group_colors_expanded)+
  savanna_theme
ggp
save_emf(ggp, "expression_density_plot.emf")



#_______4.5 Top 10 Upregulated genes Box Plot_________________________________________
candidate_genes_up = rownames(master_sig_up)
candidate_genes_up10 = rownames(master_sig_up)[1:10]
ggp = multi_boxplot(em_symbols, candidate_genes_up10, ss$SAMPLE_GROUP)
ggp = ggp + 
  scale_fill_manual(values = group_colors)+
  labs(title = "Top 10 Upregulated Genes")+
  savanna_theme
ggp
save_emf(ggp, "top10_up_boxplot.emf")


#______4.6 Top 10 Downregulated genes Boxplot_________________________________________________
candidate_genes_down = rownames(master_sig_down)
candidate_genes_down10 = rownames(master_sig_down)[1:10]
ggp = multi_boxplot(em_symbols, candidate_genes_down, ss$SAMPLE_GROUP)
ggp = ggp + 
  scale_fill_manual(values = group_colors)+
  labs(title = "Top 10 Downregulated Genes")+
  savanna_theme
ggp
save_emf(ggp, "top10_down_boxplot.emf")


#______4.7 Heatmap of All significant genes_____________________________________
ggp = make_heatmap(em_symbols, sig_genes)
ggp = ggp +
  savanna_theme
ggp
save_emf(ggp, "heatmap_sig.emf")

#_______4.8 Heatmap of significantly up-regulated genes__________________________
ggp = make_heatmap(em_symbols, candidate_genes_up)
ggp = ggp +
  labs(title = "Upregulated Genes")+
  savanna_theme
ggp
save_emf(ggp, "heatmap_upsig.emf")

#_______4.9 Heatmap of significantly down-regulated genes_______________________
ggp = make_heatmap(em_symbols, candidate_genes_down)
ggp = ggp +
  labs(title = "Downregulated Genes")+
  savanna_theme
ggp
save_emf(ggp, "heatmap_downsig.emf")

#_______4.10 Single Gene Plots__________________________________________________
for (gene in genes_of_interest)
{
  ggp = single_boxplot(em_symbols, gene, ss$SAMPLE_GROUP)
  ggp = ggp + 
    scale_fill_manual(values = group_colors)+
    savanna_theme
  ggp
  save_emf(ggp, paste(gene, "_boxplot.emf", sep=""))
}

#_________4.11 Box Plot to compare target genes_____ 
ggp = multi_boxplot(em_symbols, genes_of_interest, ss$SAMPLE_GROUP)
ggp = ggp +
  scale_fill_manual(values = group_colors)+
  labs(title = "Key P53 Pathway Genes")+
  savanna_theme 
ggp
save_emf(ggp, "p53_pathway_genes.emf")


#______5.0 Rug _________________________________________________________________
ss$SAMPLE_GROUP = factor(ss$SAMPLE_GROUP, levels = names(group_colors))
rug_data = melt(as.matrix(as.numeric(ss$SAMPLE_GROUP)))

rug = ggplot(rug_data, aes(x=Var1, y=Var2, fill=value))+
  geom_tile(show.legend = FALSE)+
  scale_fill_gradientn(colours = group_colors)+
  theme_void()
rug
save_emf(rug, "rug.emf")


#______6.0 Pathway Analysis_____________________________________________________

#______6.1 Pathway Analysis on All significant genes______________________
pathway_all = do_pathway(org.Hs.eg.db,sig_genes,em_symbols,ss$SAMPLE_GROUP,"SYMBOL")

# Extract bar plot
ggp = pathway_all$plots$ggp.barplot +
  labs(title = "MDM2i Top 10 Enriched Biological Processes", x = "Gene Count")+ savanna_theme
ggp
save_emf(ggp, "pathway_barplot.emf")

#Extract dot plot
ggp = pathway_all$plots$ggp.dotplot +
  labs(title = "MDM2i Top 10 Enriched Biological Processes")+ savanna_theme
ggp
save_emf(ggp, "pathway_dotplot.emf")

#Extract genes for apoptotic pathway 
ggp= pathway_all$plots$ontology2$ggp.boxplot+
  scale_fill_manual(values = group_colors)+
  labs(title = "Apoptotic signalling by P53")
ggp
save_emf(ggp, "apoptosis_genes.emf")

#______6.2 Pathway Analysis on significantly Up-regulated genes______________________
pathway_up = do_pathway(org.Hs.eg.db,candidate_genes_up,em_symbols,ss$SAMPLE_GROUP,"SYMBOL")

# Extract bar plot
ggp = pathway_up$plots$ggp.barplot +
  labs(title = "MDM2i Top 10 UPregulated Biological Processes", x = "Gene Count")+ savanna_theme
ggp
save_emf(ggp, "pathway_Up_barplot.emf")

# Extract dot plot
ggp = pathway_up$plots$ggp.dotplot +
  labs(title = "MDM2i Top 10 Upregulated Biological Processes", x = "Gene Count")+savanna_theme
ggp
save_emf(ggp, "pathway_up_dotplot.emf")


#______6.3 Pathway Analysis on significantly Down-regulated genes____________________
pathway_down = do_pathway(org.Hs.eg.db,candidate_genes_down,em_symbols,ss$SAMPLE_GROUP,"SYMBOL")

# Extract bar plot
ggp = pathway_down$plots$ggp.barplot +
  labs(title = "MDM2i Top 10 Downregulated Biological Processes")+ savanna_theme
ggp
save_emf(ggp, "pathway_down_barplot.emf")

# Extract dot plot
ggp = pathway_down$plots$ggp.dotplot +
  labs(title = "MDM2i Top 10 Downregulated Biological Processes")+savanna_theme
save_emf(ggp, "pathway_down_dotplot.emf")



#_____6.4 chromosomal distribution of significant genes for each comparison
ggp = make_chr_plot(master[sig_genes,], "Chromosomal Distribution - MDM2I vs DMSO")
ggp
save_emf(ggp, "chr_distribution_mdm2i.emf")
























