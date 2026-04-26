###############################################################################
#DE and MDE Analysis Pipeline - Functions
# ─────────────────────────────────────────────────────────────────────────────
# Description : A collection of reusable functions for RNA-seq data analysis.
#               Covers: plot saving, PCA, heatmaps, boxplots, pathway 
#               enrichment, metagene and signature analysis.
#
# Functions included:
#   savanna_theme  - custom theme with earthy colors
#   save_emf        – Saves a ggplot object as an EMF file
#   make_pc1_pc2    – Generates a PCA plot from an expression matrix
#   make_heatmap    – Generates a clustered heatmap from an expression matrix
#   multi_boxplot   – Generates a boxplot for multiple genes
#   single_boxplot  – Generates a boxplot for a single gene
#   do_pathway      – Runs ORA pathway enrichment analysis
#   make_metagene   – Generates a metagene boxplot from a scaled expression matrix
#   plot_signature  – Runs a full signature analysis
#   make_chr_plot   - Generates chromosome signature
#
#
# Author  : Nelly Nyaga
# Date    : 19/03/2026
###############################################################################

#_________Custom theme___________________________________________________________

savanna_theme = theme(
  text                = element_text(family = "Roboto"),
  plot.title          = element_text(size = 12, face = "bold", color = "black",
                                     hjust = 0.5, margin = margin(b=20)),
  axis.title          = element_text(size = 12, face = "bold", color = "#3E2723"),
  axis.text           = element_text(size = 10, color = "#3E2723"),
  strip.background    = element_rect(fill = "white", color = "black"),
  strip.text          = element_text(face = "bold", margin = margin(5, 5, 5, 5)),
  panel.background    = element_blank(),
  panel.grid.major.y  = element_blank(),
  panel.grid.major.x  = element_blank(),
  panel.grid.minor    = element_blank(),
  axis.line           = element_line(color = "#3E2723", linewidth = 0.6),
  legend.position     = "right",
  panel.spacing       = unit(1.5, "lines")
)

#__________________PCA1_PCA2 Function___________________________________________ 
make_pc1_pc2 = function(colour_groups, e_data)
{
  #load libraries 
  library("ggplot2")
  # scale data
  e_data_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  # run PCA
  xx = prcomp(t(e_data_scaled))
  pca_coordinates = data.frame(xx$x)
  
  # get % variation
  vars = apply(xx$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
  y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")
  # plot 
  ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = colour_groups)) +
    geom_point(size = 3) +
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) +
    theme_bw()
  return(ggp)
}
#______________Save plot function_______________________________________________
save_emf = function(plot_obj, file_name, h = 7, w = 7) 
  {
  library(devEMF)
  path = file.path(output_dir, file_name)
  emf(path, height = h, width = w)
  print(plot_obj)
  dev.off()
  message("Saved: ", path)
}
#____________Heatmap function___________________________________________________ 
make_heatmap = function(e_data, candidate_genes)
  
{
  # load libraries
  library(amap)
  library(reshape2)
  library(ggplot2)
  
  # parses
  em_scaled_candidates = na.omit(data.frame(t(scale(t(e_data[candidate_genes,])))))
  hm.matrix = as.matrix(em_scaled_candidates)
  # does the y clustering
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  
  # melt and plot
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  # create heatmap
  ggp <- ggplot(hm.matrix_clustered, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low="#001F7A", mid="#000000", high="#FFD600", midpoint=0) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 4),
      axis.title = element_blank()
    )
  
  # return plot
  return(ggp)
}
#______Single Gene Box Plot_____________________________________________________
single_boxplot = function(e_data, gene, g_data)
{
  library(reshape2)
  library(ggplot2)
  
  # create the table
  gene_data = data.frame(expression = as.numeric(e_data[gene,]))
  gene_data$groups = g_data
  
  # plot
  ggp = ggplot(gene_data, aes(x=groups, y=expression, fill=groups))+
    geom_boxplot()+
    theme_classic()+
    labs(title = gene, x = "Group", y = "Expression")
  
  return(ggp)
}
#__________________Multi-gene box plot function__________________________________ 
multi_boxplot = function(e_data,candidate_genes,g_data)
{
  library(reshape2)
  library(ggplot2)
  # create the table
  gene_data = na.omit(data.frame(t(scale(t(e_data[candidate_genes,])))))
  gene_data = data.frame(t(gene_data))
  gene_data$groups = g_data
  # melt
  
  gene_data.m = melt(gene_data, id.vars = "groups")
  
  # makes the plot
  ggp = ggplot(gene_data.m, aes(x=variable, y=value, fill = groups)) +
    geom_boxplot() +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #return Plot
  return(ggp)

}

#___________Pathway analysis function___________________________________________
do_pathway = function(organism_db, genes, e_data, g_data, id_type)
{
  # pathway libraries
  library(clusterProfiler)
  library(STRINGdb)

  print("libraries done")
 
  # converts from ensembl Symbols to Entrez
  sig_genes_entrez = bitr(genes, 
                          fromType = id_type, 
                          toType = c("ENTREZID"),
                          OrgDb = organism_db)
  print("gene IDs converted")
  
  # run enrichment
  pathway_data = enrichGO(
    gene = sig_genes_entrez$ENTREZID,
    OrgDb = organism_db,
    ont = "BP",
    pvalueCutoff = 0.05,  
    qvalueCutoff = 0.1)
  
  print("enrichment done")
  
  # extract the enrichment table from go enrich
  ora_results = data.frame(
    cbind(pathway_data$geneID,
          pathway_data$Description, 
          pathway_data$p.adjust ))   
  print("enrichment table exracted")
  
  ## Get the basic plots 
  ggp.barplot = barplot(pathway_data, showCategory=10)+ 
    geom_text(aes(label = Count),hjust = -0.2, size = 3)+theme_classic()
  
  ggp.dotplot = dotplot(pathway_data, showCategory=10)+theme_classic() 
  
  ## Create the list to store the results
  pathway_results = list("tables" = list(),"plots"= list()) 
  print("list to store resulsts created")
  
  # store tables and plots
  pathway_results$tables$pathway_data = pathway_data
  pathway_results$tables$ora_results = ora_results
  pathway_results$plots$ggp.barplot = ggp.barplot 
  pathway_results$plots$ggp.dotplot = ggp.dotplot 
  
   ## Take the top 10 ontologies, and make the boxplot & heatmap,
  for (row_index in 1:10)
  {
    # skip if row doesn't exist
    if (row_index > nrow(ora_results)) next
    
    # get the genes
    enriched_gene_set = as.character(ora_results[row_index,1])
    candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
    
    # convert ENTREZ → original ID type
    converted = bitr(candidate_genes,
                     fromType = "ENTREZID",
                     toType = id_type,
                     OrgDb = organism_db)
    
    candidate_genes = unique(converted[[id_type]])
    
    # ensure genes exist in expression matrix
    candidate_genes = intersect(candidate_genes, rownames(e_data))
    
    # skip if fewer than 2 genes
    if (length(candidate_genes) < 2) next
   
    # make the plots
    ggp.heatmap = make_heatmap(e_data, candidate_genes)
    ggp.boxplot = multi_boxplot(e_data, candidate_genes, g_data)
    ## store in a list
    ontology_result = list()
    ontology_result$candidate_genes = candidate_genes
    ontology_result$ggp.heatmap = ggp.heatmap
    ontology_result$ggp.boxplot = ggp.boxplot
    # put list into results
    pathway_results$plots[[paste0("ontology",row_index)]] = ontology_result
  }
  return(pathway_results)
}
#______________Metagene Function________________________________________________
make_metagene = function(scaled_e_data,signature_genes,signature_edata,g_data)
{  
  signature_em = scaled_e_data[signature_genes,]
  signature_metagene = data.frame(colMeans(signature_edata))
  names(signature_metagene) = "meta_expression"
  signature_metagene$group = g_data
  ggp = ggplot(signature_metagene, aes(x=group, y=meta_expression, fill=group)) +geom_violin(trim=FALSE, alpha=0.7)+geom_jitter(width=0.05, size=2)+
    scale_fill_manual(values=group_colors)
  
  return(ggp)
}
#__________Plot Signature function______________________________________________
plot_signature = function() 
{ 
ggp = make_heatmap(e_data, candidate_genes) 
pathway_results = do_pathway(organism_db, genes, e_data, g_data, id_type) 
ggp2 = make_metagene(scaled_e_data,signature_genes,signature_edata,g_data) 
results = list(list("tables" = list(),"plots"= list())) 
return(results) 
} 
#________Chromosome plot function_____________________
make_chr_plot = function(sig_data, title)
{
  # count how many significant genes are on each chromosome
  chr_counts = as.data.frame(table(sig_data$chr))
  colnames(chr_counts) = c("chromosome", "count")
  
  # order chromosomes correctly 1,2,3...22,X,Y to avoid automatic alphabetic ordering
  chr_counts$chromosome = factor(chr_counts$chromosome,
                                 levels=c(as.character(1:22), "X", "Y"))
  # plot
  ggplot(chr_counts, aes(x=chromosome, y=count, fill=count))+
    geom_bar(stat="identity")+
    scale_fill_gradient(low="blue", high="red")+
    labs(title=title, x="Chromosome", y="Number of Significant Genes")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1))
}






