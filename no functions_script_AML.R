##############################################
#Acute Myeloid Leukemia RNA-Seq Data Analysis
##############################################

###############
# 1. Load Data 
###############
em = read.table("C:/Users/nelly/Desktop/University of Glasgow/R/AML/EM.csv", header= TRUE, row.names=1,sep="\t")
de = read.table("C:/Users/nelly/Desktop/University of Glasgow/R/AML/DE_control_vs_combo.csv", header= TRUE, row.names=1,sep="\t")
annotations = read.table("C:/Users/nelly/Desktop/University of Glasgow/R/AML/annotations.csv", header= TRUE, row.names=1,sep="\t")
ss = read.table("C:/Users/nelly/Desktop/University of Glasgow/R/AML/sample_sheet.csv", header= TRUE, row.names=1,sep="\t")

#combine the data tables to form master 
master_temp = merge(em, annotations, by.x=0, by.y=0) 
master  = merge(master_temp,de,by.x=1, by.y=0)

#delete column with Ensemble IDs and set new row names 
master = master [-1]
row.names(master)=master[,"symbol"]

#em table with gene symbols as row names not Ensemble ID
em_symbols = master[,1:6]

######################
# 2. Advanced Parsing 
######################

# Remove rows with NA values from the master table
master = na.omit(master)

# Sort the master table by p.adj value
sorted = order(master[,"p.adj"],decreasing = FALSE)
master = master [sorted,]

# Add a column to master containing mean expression per gene
master$means =rowMeans(master[, 1:6])

# Add a column to master containing −log10(p)
master$mlog10p = -log10(master$p.adj)

# Add a column to master flagging significant genes as T/F
master$sig = as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1.0) 

# Create a scaled expression matrix from em_symbols
em_scaled = data.frame(t(scale(t(em_symbols))))
em_scaled =na.omit(em_scaled)

# Subset the master table to retain significant genes only
master_sig = subset(master, p.adj < 0.05 & abs(log2fold) >1)

# Create a vector of significant gene names by row names
sig_genes = master_sig$gene_name
sig_genes = row.names(master_sig)

# Create an expression table of significant genes only
em_symbols_sig = em_symbols[sig_genes,]

# Create a scaled expression table of significant genes only
em_scaled_sig = em_scaled[sig_genes,]

# Save all key tables to disk
write.table(em_symbols,file = "C:/Users/nelly/Desktop/University of Glasgow/R/AML/em_symbols.csv",sep="\t")
write.table(em_scaled,file = "C:/Users/nelly/Desktop/University of Glasgow/R/AML/em_scaled.csv", sep="\t")
write.table(master_sig,file = "C:/Users/nelly/Desktop/University of Glasgow/R/AML/master_sig.csv",sep="\t")
write.table(master,file = "C:/Users/nelly/Desktop/University of Glasgow/R/AML/master.csv",sep="\t")
write.table(em_symbols_sig,file = "C:/Users/nelly/Desktop/University of Glasgow/R/AML/em_symbols_sig.csv",sep="\t")
write.table(em_scaled_sig,file = "C:/Users/nelly/Desktop/University of Glasgow/R/AML/em_scaled_sig.csv",sep="\t")



#########################
# 3. Plots
#########################
#Load Libraries
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(clusterProfiler)  
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(STRINGdb)
#install.packages("devEMF")
library(devEMF)

#creates theme
savanna_theme = theme(
  text = element_text(family = "Roboto"),
  plot.title = element_text(size = 12, face = "bold", color = "black", hjust = 0.5, margin = margin(b=20)),
  axis.title = element_text(size = 12, face = "bold", color = "#3E2723"),
  axis.text = element_text(size = 10, color = "#3E2723"),
  strip.background = element_rect(fill = "white", color = "black"),
  strip.text = element_text(face = "bold", margin = margin(5, 5, 5, 5)),
  panel.background = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "#3E2723", linewidth = 0.6),
  legend.position = "right",
  panel.spacing = unit(1.5, "lines")
)

savanna_colors = c("control" = "#4E342E","combo" = "#F57F17")


#Significantly Upregulated Genes
master_sig_up = subset(master,p.adj <0.05 & log2fold > 1 )
#master_sig_up = master_sig_up[order(master_sig_up$log2fold, decreasing = TRUE), ]
master_sig_up_5 = master_sig_up[1:5, ]

#Significantly Downregulated genes
master_sig_down = subset(master, p.adj<0.05 & log2fold < -1)
master_sig_down_5 = master_sig_down[1:5, ]
  
#Add new columns to master
master_non_sig =subset(master,sig == FALSE)
master_non_sig$direction = "a" 
master_sig_up$direction = "b" 
master_sig_down$direction ="c" 
master=rbind(master_non_sig,master_sig_up,master_sig_down)

#count number of upregulated and downregulated genes to add to volcano plot.
n_up= sum(master$direction == "b")
n_down =sum(master$direction == "c")


#Volcano Plot
ggp1 = ggplot(master, aes(x=log2fold, y=mlog10p,color=direction))+
  geom_point()+
  labs(title="Volcano Plot: Treatment Combo vs Control", x="log2(Fold Change)",y="-log10(p.adj)")+
  savanna_theme+
  xlim(c(-10,10))+
  ylim(c(0,70))+
  scale_colour_manual(values=c("#4E342E", "#F57F17", "#FFB300"),labels=c("Unchanged","Upregulated","Downregulated"),name = "Legend")+
  annotate("text", x = 7, y = 65, label = paste(n_up), color = "black", fontface="bold",size = 3, hjust = 2)+
  annotate("text",x = -7, y = 65,label = paste(n_down),color = "black",fontface="bold",size = 3, hjust = 2)+
  geom_label_repel(data=master_sig_up_5,aes(label=symbol),color="black",fill="#F57F17",label.padding=0.2,label.size=0.25,force =2,max.overlaps = 20, show.legend =FALSE)+ 
  geom_label_repel(data=master_sig_down_5,aes(label=symbol),color="black",fill= "#FFB300",label.padding=0.2,label.size=0.25,force =2,max.overlaps = 20,show.legend = FALSE)
ggp1

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/volcano_plot.emf", height = 7, width = 7)
print(ggp1) 
dev.off() 

#MA plot
ggp2 =ggplot(master,aes(x=log10(means),y=log2fold,colour = direction))+
  geom_point()+ 
  labs(title="MA Plot: Treatment Combo Vs. Control", x="Mean Expression(log10)",y="Log2 fold change")+
  savanna_theme+
  xlim(c(0,5))+
  ylim(c(-10,+10))+
  scale_colour_manual(values=c("#4E342E", "#F57F17", "#FFB300"),labels=c("Unchnaged","Upregulated","Downregulated"),name = "Legend")+
  geom_text_repel(data=master_sig_up_5,aes(label=symbol),color="black", show.legend =FALSE)+ 
  geom_text_repel(data=master_sig_down_5,aes(label=symbol),color="black",show.legend = FALSE)
ggp2

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ma_plot.emf", height = 7, width = 7) 
print(ggp2) 
dev.off() 

  
# PCA plot (Global sample structure)

em_symbols = master[ ,rownames(ss)]

#convert em_symbols into a matrix and scale
em_symbols = as.matrix(sapply(em_symbols,as.numeric))
pca = prcomp(t(em_symbols), scale. = TRUE)

pca_coordinates =data.frame(pca$x)
pca_coordinates=cbind(pca_coordinates, ss) 
vars = apply(pca$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100

x_axis_label = paste("PC1"," (",prop_x, "%)",sep="")
y_axis_label = paste('PC2'," (",prop_y, "%)",sep="")

pca_coordinates$sample_group = factor(pca_coordinates$sample_group,levels = c("control", "combo"))

ggp3=ggplot(pca_coordinates,aes(x=PC1,y=PC2,color= sample_group))+
  geom_point(size=5)+ 
  scale_color_manual(values = savanna_colors)+ #color dots
  savanna_theme+
  labs(title = "PCA", x = x_axis_label, y= y_axis_label)
ggp3

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/pca_plot.emf", height = 7, width = 7) 
print(ggp3) 
dev.off() 

# Expression density plots (Sample quality)
em.m = melt((em)) #melt the em table. 
em.m$sample_group = gsub("_[1-3]", "", em.m$variable)

savanna_colors_expanded = c("CONTROL_1" = "#4E342E", "CONTROL_2" = "#4E342E", "CONTROL_3" = "#4E342E",
  "RITA_CPI_COMBO_1" = "#F57F17", "RITA_CPI_COMBO_2" = "#F57F17", "RITA_CPI_COMBO_3" = "#F57F17")

ggp4 = ggplot(em.m,aes(x=log10(value+0.01)))+
  geom_density(aes(fill = variable), linewidth = 0.5, alpha=0.75, show.legend = FALSE)+
  labs(x="Expression(log10)",y="Density")+
  facet_wrap(~variable,ncol = 3,scales ="free")+
  scale_x_continuous(breaks = seq(-2, 6, by = 2))+ #for all x-axis to start at -2
  scale_fill_manual(values = savanna_colors_expanded)+
  savanna_theme
ggp4

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/expression_density_plot.emf", height = 7, width = 7) 
print(ggp4) 
dev.off() 


# Single-gene plots
em_symbols = master[ ,rownames(ss)]
gene_data = em_symbols["FLT3",]
gene_data=data.frame(t(gene_data)) #transform and make it into a dataframe
gene_data$sample_group=ss$sample_group #add a column to the data called sample group from the ss table 
names(gene_data) = c("expression","sample_group")#renames the columns in the table 
gene_data$sample_group=factor(gene_data$sample_group,levels = c("control","combo"))#sets the order of the groups 


##FLT3 Violin-Box Plot
ggp5 = ggplot(data=gene_data,aes(x=sample_group,y=expression,fill=sample_group))+
  geom_violin(width=0.7,alpha=0.5,trim=TRUE,show.legend = FALSE)+
  geom_boxplot(color="black", width=0.1,alpha=0.5,show.legend = FALSE)+
  scale_fill_manual(values = savanna_colors)+
  labs(title="FLT3 Expression in Control vs Treatment Groups", x= "Sample_group", y="Expression")+
  savanna_theme
ggp5

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/violin_box_plot.emf", height = 7, width = 7) 
print(ggp5) 
dev.off() 


#FLT3 Box-Jitter combo
ggp6 = ggplot(data=gene_data,aes(x=sample_group,y=expression,fill=sample_group))+
  geom_boxplot(width=0.5,alpha=0.5,show.legend = FALSE)+
  geom_jitter(aes(color=sample_group),size=4,show.legend = FALSE)+
  scale_fill_manual(values = savanna_colors)+
  scale_color_manual(values = savanna_colors)+
  labs(title="FLT3 Expression in Control vs Treatment Groups", x= "Sample_group", y="Expression")+
  savanna_theme
ggp6

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/FLT3_box_jitter_plot.emf", height = 7, width = 7) 
print(ggp6) 
dev.off() 

#DNMT3A
em_symbols = master[ ,rownames(ss)]
gene_data = em_symbols["DNMT3A",]
gene_data=data.frame(t(gene_data)) #transform and make it into a dataframe
gene_data$sample_group=ss$sample_group #add a column to the data called sample group from the ss table 
names(gene_data) = c("expression","sample_group")#renames the columns in the table 
gene_data$sample_group=factor(gene_data$sample_group,levels = c("control","combo"))#sets the order of the groups right cz if you plot without its is not right

#plot
set.seed(42) #ensures the jitter points remain the same
ggp7 = ggplot(data=gene_data,aes(x=sample_group,y=expression,fill=sample_group))+
  geom_boxplot(width=0.5,alpha=0.5,show.legend = FALSE)+
  geom_jitter(aes(color=sample_group),size=4,show.legend = FALSE)+
  scale_fill_manual(values = savanna_colors)+
  scale_color_manual(values = savanna_colors)+
  labs(title="DNMT3A Expression in Control vs Treatment Groups",x= "Sample_group", y="Expression")+
  savanna_theme
ggp7

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/dnmt3a_box_jitter_plot.emf", height = 7, width = 7) 
print(ggp7) 
dev.off() 

#NPM1
em_symbols = master[ ,rownames(ss)]
gene_data = em_symbols["NPM1",]
gene_data=data.frame(t(gene_data)) #transform and make it into a dataframe
gene_data$sample_group=ss$sample_group #add a column to the data called sample group from the ss table 
names(gene_data) = c("expression","sample_group")#renames the columns in the table 
gene_data$sample_group=factor(gene_data$sample_group,levels = c("control","combo"))#sets the order of the groups right cz if you plot without its is not right

#plot
ggp8 = ggplot(data=gene_data,aes(x=sample_group,y=expression,fill=sample_group))+
  geom_boxplot(width=0.5,alpha=0.5,show.legend = FALSE)+
  geom_jitter(aes(color=sample_group),size=4,show.legend = FALSE)+
  scale_fill_manual(values= savanna_colors)+
  scale_color_manual(values = savanna_colors)+
  labs(title="NPM1 Expression in Control vs Treatment Groups",x= "Sample_group", y="Expression")+
  savanna_theme
ggp8

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/NPM1_box_jitter_plot.emf", height = 7, width = 7) 
print(ggp8) 
dev.off() 

################################
# Multigene Plots 

#top10 Up-regulated Genes
top10_up = master_sig_up[1:10,]
candidate_genes = rownames(top10_up)

gene_data =em_symbols[candidate_genes,]
gene_data=data.frame(t(gene_data)) #transform and make it into a dataframe
gene_data$sample_group=ss$sample_group #add a column to the data called sample group from the ss table (Check table, has 10 genes, is transposed and has a column with the group names)
gene_data$sample_group=factor(gene_data$sample_group,levels = c("control","combo"))
#melting
gene_data.m = melt(gene_data,id.vars="sample_group")

#plotting
ggp9 = ggplot(data=gene_data.m,aes(x=variable,y=value,fill=sample_group))+
  geom_boxplot(width=1,alpha=0.5,show.legend = TRUE)+
  scale_fill_manual(values = savanna_colors)+
  labs(title="Top 10 Upregulated Genes",x ="Gene",y="Gene Expression Level")+
  savanna_theme+
  theme(axis.text.x =element_text(angle=45,hjust=1))
ggp9

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/top10_up_box_plot.emf", height = 7, width = 7) 
print(ggp9) 
dev.off() 

#Top 10 downregulated genes
top10_down = master_sig_down[1:10,]
candidate_genes = rownames(top10_down)

#transforming and adding column
gene_data =em_symbols[candidate_genes,]
gene_data=data.frame(t(gene_data)) #transform and make it into a dataframe
gene_data$sample_group=ss$sample_group #add a column to the data called sample group from the ss table (Check table, has 10 genes, is transposed and has a column with the group names)
gene_data$sample_group=factor(gene_data$sample_group,levels = c("control","combo"))

#melting
gene_data.m = melt(gene_data,id.vars="sample_group")

#plotting
ggp10 = ggplot(data=gene_data.m,aes(x=variable,y=value,fill=sample_group))+
  geom_boxplot(width=1,alpha=0.5,show.legend = TRUE)+
  scale_fill_manual(values = savanna_colors)+
  labs(title="Top 10 Downregulated Genes",x ="Gene",y="Gene Expression Level")+
  savanna_theme+
  theme(axis.text.x =element_text(angle=45,hjust=1))
ggp10

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/top10_down_box_plot.emf", height = 7, width = 7) 
print(ggp10) 
dev.off() 

#=====================================================
#Faceted Boxplot for the 10 most significant genes
top10 = master_sig[1:10,]
candidate_genes = rownames(top10)

#transforming and adding column
gene_data =em_symbols[candidate_genes,]
gene_data=data.frame(t(gene_data)) 
gene_data$sample_group=ss$sample_group 
gene_data$sample_group=factor(gene_data$sample_group,levels = c("control","combo"))

#Melting
gene_data.m = melt(gene_data,id.vars="sample_group")

#plotting
ggp11 = ggplot(data=gene_data.m,aes(x=sample_group,y=value, fill=sample_group))+
  geom_boxplot(width=0.5,alpha=0.5,show.legend = FALSE)+
  facet_wrap(~variable,ncol=5)+
  scale_fill_manual(values = savanna_colors)+
  labs(title="Top 10 Significant Genes",x ="Sample Group",y="Gene Expression Level")+
  savanna_theme+
  theme(axis.text.x =element_text(angle=45,hjust=1))
ggp11

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/faceted_box_plot.emf", height = 7, width = 7) 
print(ggp11) 
dev.off() 


# Heatmap for all significant genes
#matrix
hm.matrix=as.matrix(em_scaled_sig)

#distance between the genes 
y.dist=Dist(hm.matrix,method="spearman")

#cluster the genes
y.cluster = hclust(y.dist, method="average")

#pull out the dendogram
y.dd =as.dendrogram(y.cluster)

#untangle the dendogram 
y.dd.reorder = reorder(y.dd,0,FUN="average")

#get the untagled gene order from the dendogran
y.order = order.dendrogram(y.dd.reorder)

#reorder the original matrix in the new order
hm.matrix_clustered =hm.matrix[y.order,]

#make color palette
colors =c("#00F5FF","#121212","#FF4D00")
palette =colorRampPalette(colors)(100)

#melt and plot
hm.matrix_clustered=melt(hm.matrix_clustered)

ggp12 = ggplot(hm.matrix_clustered,aes(x=Var2, y=Var1, fill=value))+
  geom_tile()+
  scale_fill_gradientn(colors=palette)+
  savanna_theme+
  theme(plot.margin=unit(c(0,1,1,1), "cm"), 
        axis.line=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 2),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="right",panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())
ggp12

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/heat_map.emf", emfPlus = FALSE, height = 4, width = 6) 
print(ggp12) 
dev.off() 

# Rug
ss$sample_group = factor(ss$sample_group, levels = c("control", "combo"))
rug_colours = savanna_colors
rug_data = as.matrix(as.numeric(factor(ss$sample_group)))
rug_data = melt(rug_data)

ggp13 = ggplot(rug_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(show.legend = FALSE)+
  scale_fill_gradientn(colours = rug_colours)+
  theme_void()
ggp13

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/rug.emf", height = 7, width = 7) 
print(ggp13) 
dev.off() 

####################################################################################
# Pathway enrichment plot ORA
#####################################################################################
#Convert the IDs
sig_genes_entrez =bitr(row.names(master_sig),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Over Representation Analysis - Biological Process
ora_results = enrichGO(gene=sig_genes_entrez$ENTREZID,OrgDb = org.Hs.eg.db,readable = T,ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10)

ggp14 = barplot(ora_results,showCategory=10)+
  scale_fill_gradient(low = "#4E342E", high = "#F57F17")+
  labs(title = "Top 10 Enriched Biological Processes", x = "Gene Count")+
  geom_text(aes(label = Count), hjust = -0.2, size = 3, color = "black")+
  savanna_theme
ggp14
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_bp_barplot.emf", height = 7, width = 7) 
print(ggp14) 
dev.off() 

ggp15 = dotplot(ora_results,showCategory=10)+
  scale_fill_gradient(low = "#4E342E", high = "#F57F17")+
  labs(title = "Top 10 Enriched Biological Processes")+
  savanna_theme
ggp15
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_bp_dotplot.emf", height = 7, width = 7) 
print(ggp15) 
dev.off() 

ggp16 = goplot(ora_results,showCategory=5)
ggp16
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_bp_goplot.emf", height = 7, width = 7) 
print(ggp16) 
dev.off() 


ggp = cnetplot(ora_results, showCategory = 5)+
  labs(title = "Top 5 Enriched Biological Processes Network")+
  savanna_theme
ggp
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/cnet.emf", height = 7, width = 7) 
print(ggp) 
dev.off()

#ORA - Molecular Functions 
ora_results = enrichGO(gene=sig_genes_entrez$ENTREZID,OrgDb = org.Hs.eg.db,readable = T,ont="MF",pvalueCutoff = 0.05,qvalueCutoff = 0.10)

#Plot
ggp17 = barplot(ora_results,showCategory=10)+
  labs(title = "Top 10 Enriched Molecular Functions", x = "Number of Genes")+
  scale_fill_gradient(low = "#FFB300", high = "#4E342E")+
  savanna_theme
ggp17
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_mf_barplot.emf", height = 7, width = 7) 
print(ggp17) 
dev.off() 

#plot
ggp18 = dotplot(ora_results,showCategory=10)+
  labs(title = "Top 10 Enriched Molecular Functions")+
  scale_fill_gradient(low = "#FFB300", high = "#4E342E")+
  savanna_theme
ggp18
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_mf_dotplot.emf", height = 7, width = 7) 
print(ggp18) 
dev.off() 

#plot
ggp19 = goplot(ora_results,showCategory=10)
ggp19
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_mf_goplot.emf", height = 7, width = 7) 
print(ggp19) 
dev.off() 

ggp20 = cnetplot(ora_results, showCategory = 5)+
  labs(title = "Top 5 Enriched Molecular Functions Network")+
  savanna_theme
ggp20
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/cnet_mf.emf", height = 7, width = 7) 
print(ggp20) 
dev.off() 

#ORA for Significantly Up-regulated and Down regulated Biological Processes

#Upregulated
upsig_genes_entrez =bitr(row.names(master_sig_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ora_results_up = enrichGO(gene=upsig_genes_entrez$ENTREZID,OrgDb = org.Hs.eg.db,readable = T,ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10)

ggp21 = barplot(ora_results_up,showCategory=10)+
  labs(title = "Top 10 Upregulated Biological Processes", x = "Number of Genes")+
  geom_text(aes(label = Count), hjust = -0.2, size = 3, color = "black")+
  scale_fill_gradient(low = "#FFB300", high = "#2E7D32")+
  savanna_theme
ggp21
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_bp_upreg_barplot.emf", height = 7, width = 7) 
print(ggp21) 
dev.off() 
  
#Downregulated
downsig_genes_entrez =bitr(row.names(master_sig_down),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ora_results_down = enrichGO(gene=downsig_genes_entrez$ENTREZID,OrgDb = org.Hs.eg.db,readable = T,ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10)

ggp22 = barplot(ora_results_down,showCategory=10)+
  labs(title = "Top 10 Down-regulated Biological Processes", x = "Number of Genes")+
  scale_fill_gradient(low = "#FFB300", high = "#2E7D32")+
  geom_text(aes(label = Count), hjust = -0.2, size = 3, color = "black")+
  savanna_theme
ggp22
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/ora_bp_downreg_barplot.emf", height = 7, width = 7) 
print(ggp22) 
dev.off() 

#####################################################################################
#GSEA
######################################################################################
#master into a vector of genes sorted by log2fold change 
gsea_input =master$log2fold

#add names to the vector 
names(gsea_input) = row.names(master)

#omitNA values 
gsea_input=na.omit(gsea_input)

#sort list in decreasing order the clusterProfiler
gsea_input = sort(gsea_input,decreasing = TRUE)
gsea_input

#Running the GSEA function 
gse_results = gseGO(geneList=gsea_input,
                    ont ="BP",
                    keyType = "SYMBOL",
                    nPerm = 10000,
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    OrgDb = org.Hs.eg.db,
                    pAdjustMethod = "none")


#Plot
ggp23 =ridgeplot(gse_results,showCategory = 15)+
  labs(title="GSEA:Distribution of Fold Change",x ="Log2 Fold Change")+
  savanna_theme
ggp23
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/GSEA.emf", height = 7, width = 7) 
print(ggp23) 
dev.off() 


#######################################################################################
#GENE Lists
######################################################################################
#List in Downregulated Mono nuclear cell proliferation
gene_sets = ora_results_down$geneID
description = ora_results_down$Description
p.adj = ora_results_down$p.adjust

ora_results_table = data.frame(cbind(gene_sets,p.adj), row.names =description)

#extract the list of most enriched genes from the ontology 
enriched_gene_set = as.character(ora_results_table[3,"gene_sets"])
#split into vector 
candidate_genes = unlist(strsplit(enriched_gene_set,"/"))
candidate_genes

#create the table
gene_data = em_scaled[candidate_genes[1:10], ]
gene_data = data.frame(t(gene_data))
gene_data$groups = ss$sample_group

#set the order of the groups 
gene_data$groups=factor(gene_data$groups,levels = c("control","combo"))

# melt
gene_data.m = melt(gene_data, id.vars = "groups")


#Pathway-specific gene expression Box Plot - Cell-cell Adhesion

ggp24 = ggplot(gene_data.m,aes(x=groups,y=value, fill=groups))+
  geom_boxplot(width=0.5,alpha=0.5,show.legend = FALSE)+
  facet_wrap(~variable,ncol=5)+
  scale_fill_manual(values = savanna_colors)+
  labs(title = "Top 10 Genes in Mononuclear Cell Proliferation", y = "Expression")+
  savanna_theme+
  theme(axis.text.x =element_text(angle=45,hjust=1))
ggp24

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/mononuclear_genes.emf", height = 7, width = 7) 
print(ggp24) 
dev.off() 

#Heat Map of the Mononuclear cell proliferation genes
hm.matrix = as.matrix(em_scaled[candidate_genes[1:20],])

#y clustering
y.dist = Dist(hm.matrix, method="spearman")
y.cluster = hclust(y.dist, method="average")
y.dd = as.dendrogram(y.cluster)
y.dd.reorder = reorder(y.dd,0,FUN="average")
y.order = order.dendrogram(y.dd.reorder)

# re-order
hm.matrix_clustered = hm.matrix[y.order,]

# melt
hm.matrix_clustered = melt(hm.matrix_clustered)
# colour palette
colours = c("#3E2723", "#F5F5DC", "#D35400")
palette = colorRampPalette(colours)(100)

# plot
ggp25 = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value))+
  geom_tile()+
  scale_fill_gradientn(colours = palette)+ 
  labs(subtitle = paste("Top Genes in:", description[3]))+ 
  savanna_theme+
  theme(plot.margin=unit(c(0,1,1,1), "cm"), 
      axis.line=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 2),axis.title.x=element_blank(),axis.text.y=element_text(size=5),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="right",panel.background=element_blank(),panel.border=element_blank(),
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())
ggp25
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/mononuclear_heat_map.emf", height = 7, width = 7) 
print(ggp25) 
dev.off() 

#Faceted Box Plot of Signal Transduction by P53 genes
gene_sets = ora_results_up$geneID
description = ora_results_up$Description
p.adj = ora_results_up$p.adjust

ora_results_table = data.frame(cbind(gene_sets,p.adj), row.names =description)

#extract the list of most enriched genes from the ontology 
enriched_gene_set = as.character(ora_results_table[8,"gene_sets"])
#split into vector 
candidate_genes = unlist(strsplit(enriched_gene_set,"/"))
candidate_genes

#create the table
gene_data = em_scaled[candidate_genes[1:10], ]
gene_data = data.frame(t(gene_data))
gene_data$groups = ss$sample_group

#set the order of the groups 
gene_data$groups=factor(gene_data$groups,levels = c("control","combo"))

# melt
gene_data.m = melt(gene_data, id.vars = "groups")


#Plot Faceted Box Plot 

ggp27 = ggplot(gene_data.m,aes(x=groups,y=value, fill=groups))+
  geom_boxplot(width=0.5,alpha=0.5,show.legend = FALSE)+
  facet_wrap(~variable,ncol=5)+
  scale_fill_manual(values = savanna_colors)+
  labs(title = "Top 10 Genes in Signal Transduction by P53 ", y = "Expression")+
  savanna_theme+
  theme(axis.text.x =element_text(angle=45,hjust=1))
ggp27

emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/P53_upregulated_genes.emf", height = 7, width = 7) 
print(ggp27) 
dev.off() 

#P53 Signal Transduction Heat Map
hm.matrix = as.matrix(em_scaled[candidate_genes[1:10],])

#y clustering
y.dist = Dist(hm.matrix, method="spearman")
y.cluster = hclust(y.dist, method="average")
y.dd = as.dendrogram(y.cluster)
y.dd.reorder = reorder(y.dd,0,FUN="average")
y.order = order.dendrogram(y.dd.reorder)

# re-order
hm.matrix_clustered = hm.matrix[y.order,]

# melt
hm.matrix_clustered = melt(hm.matrix_clustered)
# colour palette
colours = c("darkgreen", "white", "#D35400")
palette = colorRampPalette(colours)(100)

# plot
ggp28 = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile() +
  scale_fill_gradientn(colours = palette) + 
  labs(subtitle = paste("Top Genes in:", description[8])) + 
  savanna_theme+
  theme(legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'),axis.text.x =element_text(angle=45,hjust=1), axis.text.y = element_text(size = 8),axis.title.x = element_blank(),axis.title.y = element_blank(), 
        axis.ticks=element_blank())
ggp28
emf("C:/Users/nelly/Desktop/University of Glasgow/R/AML/P53_genes_heat_map.emf", height = 7, width = 7) 
print(ggp28) 
dev.off() 

################################################################################
# STRING interaction network
################################################################################
candidate_genes_table = data.frame(candidate_genes)

names(candidate_genes_table) = "gene"

# load the database, 9606 is human
string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=700, network_type="full", input_directory="")

# map genes to the database
string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE )

#plot
string_db$plot_network(string_mapped)


save.image(file = "AML_workspace.RData")
save.image("C:/Users/nelly/Desktop/University of Glasgow/R/AML/environment.rdata")


