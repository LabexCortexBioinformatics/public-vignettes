# library(tidyverse)
# library(Seurat) #Seurat V4 Required for this script because of avg_log2FC
# library(enrichplot)
# library(clusterProfiler)
# library(DOSE)
# require(org.Mm.eg.db)
# organism = org.Mm.eg.db
# setwd("D:/LabEx_Teams/Peyron/Narco_Project/clusterProfiler/")
# 
# 
# ###################################################################################################################################
# # CLUSTER PROFILER ANALYSIS #######################################################################################################
# 
# ################################################################################
# #                                   ORA                                        #
# ################################################################################
# 
# #Obtain findmarkers table 
# Idents(clb) = "neuron_clusters"
# clb_markers_ORA = FindMarkers(clb, ident.1 = "LPS", ident.2 = "base", group.by = "condition",
#                               subset.ident = "clb_neurons", logfc.threshold = 0.7, assay = "RNA")
# 
# #If the table is present in the environment
# dataORA = clb_markers_ORA
# dataORA <- dataORA %>% rownames_to_column(var = "X") #To create a column gene and remove rownames
# 
# #If load from csv file
# dataORA = read.csv("clb_markers_ORA.csv", header = T) 
# 
# ##1st column is ID (gene name named X)
# ##3rd column is FC 
# 
# dataORA_qNSCs <- dataORA[dataORA$avg_log2FC > 0 & dataORA$p_val < 0.05,]
# dataORA_TAPs <- dataORA[dataORA$avg_log2FC < 0 & dataORA$p_val < 0.05,]
# 
# 
# ############################## ORA qNSCs ##########################################
# namesqNSCs <- as.character(dataORA_qNSCs$X)
# egoqNSCs = bitr(namesqNSCs, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
# head(egoqNSCs)
# geneListORA_qNSCs <- egoqNSCs$ENTREZID
# geneListORA_qNSCs
# 
# ############################## ORA TAPs ########################################
# namesTAPs <- as.character(dataORA_TAPs$X)
# egoTAPs = bitr(namesTAPs, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
# head(egoTAPs)
# geneListORA_TAPs <- egoTAPs$ENTREZID
# geneListORA_TAPs
# 
# ################################################################################
# #                               BP analysis                                  #
# ################################################################################
# egoqNSCs <- enrichGO(gene          = geneListORA_qNSCs,
#                   universe      = names(geneListORA_qNSCs),
#                   OrgDb         = org.Mm.eg.db,
#                   ont           = "BP",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.01,
#                   qvalueCutoff  = 0.05,
#                   readable      = TRUE)
# head(summary(egoqNSCs))
# write.csv(egoqNSCs, row.names = TRUE, file = "ORA_qNSCs.csv")
# 
# egoqNSCs2 <- simplify(egoqNSCs, measure="Wang", semData=NULL) # The simplify method allows to remove redundant terms based on GOSemSim (select one representative term from redundant ones based on higher p.adjust, which ahve similarity higher than "cutoff)
# write.csv(egoqNSCs2, row.names = T, file = "ORA_qNSCs_simplify.csv")
# 
# egoTAPs <- enrichGO(gene          = geneListORA_TAPs,
#                     universe      = names(geneListORA_TAPs),
#                     OrgDb         = org.Mm.eg.db,
#                     ont           = "BP",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05,
#                     readable      = TRUE)
# head(summary(egoTAPs))
# write.csv(egoTAPs, row.names = T, file = "ORA_TAPs.csv")
# 
# egoTAPs2 <- simplify(egoTAPs, measure="Wang", semData=NULL) # The simplify method allows to remove redundant terms based on GOSemSim (select one representative term from redundant ones based on higher p.adjust, which ahve similarity higher than "cutoff)
# write.csv(egoTAPs2, row.names = T, file = "ORA_TAPs_simplify.csv")
# 
# ################################################################################
# #                               DotPlot ORA/BP analysis                        #
# ################################################################################
# setwd("D:/LabEx_Teams/Raineteau/Pseudotime/TradeSeq/csv/")
# jpeg(filename = "LPS_ORA_qNSCs.jpg", width=1000, height=500)
# p<-dotplot(egoqNSCs, x= "GeneRatio", orderBy = "GeneRatio", showCategory=25)
# 
# print(p)
# dev.off()
# 
# jpeg(filename = "LPS_ORA_qNSCs_simplify.jpg", width=1000, height=500)
# p<-dotplot(egoqNSCs2, x= "GeneRatio", orderBy = "GeneRatio", showCategory=25)
# print(p)
# dev.off()
# 
# jpeg(filename = "LPS_ORA_TAPs.jpg", width=1000, height=500)
# p<-dotplot(egoTAPs, x= "GeneRatio", orderBy = "GeneRatio", showCategory=25)
# print(p)
# dev.off()
# 
# jpeg(filename = "LPS_ORA_TAPs_simplify.jpg", width=1000, height=500)
# p<-dotplot(egoTAPs2, x= "GeneRatio", orderBy = "GeneRatio", showCategory=25)
# print(p)
# dev.off()

# # ###############################################################################
# #                               Heatmap ORA/BP analysis                        #
# # ###############################################################################
# 
# ## Heatmap ORA/BP qNSCs #########################################################
# df1 <- read.csv('ORA_qNSCs.csv', row.names = 3) # number is the column you want to use as names for heatmap. Default is 3 as it is the full name of the GO term
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("vesicle", "cycle", "synapse")
# position = unique(grep(paste(interest,collapse="|"), row.names(df1)))
# df1 = df1[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("muscle")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df1), invert = TRUE))
# df1 = df1[position, ]
# 
# genes = lapply(df1$geneID, FUN = function(x) strsplit(x, "/")[[1]])
# uniques = unique(unlist(as.list(genes), recursive = FALSE)) # find every genes in the analysis
# names(genes) = row.names(df1) # name each list of gene with GO associated
# 
# data <- data.frame(matrix(NA,
#                           nrow = length(uniques),
#                           ncol = nrow(df1)), row.names = uniques)
# colnames(data) <- row.names(df1) # create empty dataframe
# 
# for (i in row.names(df1)){
#   data[i] = as.integer(as.logical(uniques %in% genes[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered1 <- data[order(apply(data, 1, FUN=sum), decreasing = TRUE), ] # reorder genes by sum, comment if clustering should be done
# 
# pheatmap(data_ordered1, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20,
#          main = "Genes associated with GO terms", angle_col = 45,
#          border_color = "white", legend = FALSE, filename = "ORA_BP_qNSCs_Heatmap.jpg", silent = TRUE)
# 
# ## Heatmap ORA/BP TAPs #########################################################
# df2 <- read.csv('ORA_TAPs.csv', row.names = 3) # number is the column you want to use as names for heatmap. Default is 3 as it is the full name of the GO term
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("vesicle", "cycle", "synapse")
# position = unique(grep(paste(interest,collapse="|"), row.names(df2)))
# df2 = df2[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("muscle")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df2), invert = TRUE))
# df2 = df2[position, ]
# 
# genes = lapply(df2$geneID, FUN = function(x) strsplit(x, "/")[[1]])
# uniques = unique(unlist(as.list(genes), recursive = FALSE)) # find every genes in the analysis
# names(genes) = row.names(df2) # name each list of gene with GO associated
# 
# data <- data.frame(matrix(NA,
#                           nrow = length(uniques),
#                           ncol = nrow(df2)), row.names = uniques)
# colnames(data) <- row.names(df2) # create empty dataframe
# 
# for (i in row.names(df2)){
#   data[i] = as.integer(as.logical(uniques %in% genes[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered2 <- data[order(apply(data, 1, FUN=sum), decreasing = TRUE), ] # reorder genes by sum, comment if clustering should be done
# 
# pheatmap(data_ordered2, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20,
#          main = "Genes associated with GO terms", angle_col = 45,
#          border_color = "white", legend = FALSE, filename = "ORA_BP_TAPs_Heatmap.jpg", silent = TRUE)

# # ###############################################################################
# #                               ORA/KEGG analysis                              #
# # ###############################################################################
# search_kegg_organism('mmu', by='kegg_code')
# kkqNSCs <- enrichKEGG(gene         = geneListORA_qNSCs,
#                    organism     = "mmu",
#                    pvalueCutoff = 0.01,
#                    pAdjustMethod = "BH",
#                    minGSSize = 5,
#                    maxGSSize = 500,
#                    use_internal_data = FALSE)
# head(kkqNSCs)
# write.csv(kkqNSCs, row.names = T, file = "KEGGAnalysis_LPS_ORA_qNSCs.csv")
# 
# search_kegg_organism('mmu', by='kegg_code')
# kkTAPs <- enrichKEGG(gene         = geneListORA_TAPs,
#                      organism     = "mmu",
#                      pvalueCutoff = 0.01,
#                      pAdjustMethod = "BH",
#                      minGSSize = 5,
#                      maxGSSize = 500,
#                      use_internal_data = FALSE)
# head(kkTAPs)
# write.csv(kkTAPs, row.names = T, file = "KEGGAnalysis_LPS_ORA_TAPs.csv")

# ################################################################################
# #                              DotPlot ORA/KEGG analysis                       #
# ################################################################################
# jpeg(filename = "KEGGAnalysis_LPS_ORA_qNSCs.jpg", width=600, height=500)
# p<-dotplot(kkqNSCs, x= "GeneRatio", orderBy = "GeneRatio", showCategory=25)
# print(p)
# dev.off()
# 
# jpeg(filename = "KEGGAnalysis_LPS_ORA_TAPs.jpg", width=600, height=500)
# p<-dotplot(kkTAPs, x= "GeneRatio", orderBy = "GeneRatio", showCategory=25)
# print(p)
# dev.off()

# ################################################################################
# #                              Heatmap ORA/KEGG analysis                       #
# ################################################################################
# library(pheatmap)
# 
# ## Heatmap ORA/KEGG qNSCs #########################################################
# df3 <- read.csv('KEGGAnalysis_LPS_ORA_qNSCs.csv', row.names = 3) # number is the column you want to use as names for heatmap. Default is 3 as it is the full name of the GO term
# genes = lapply(df3$geneID, FUN = function(x) strsplit(x, "/")[[1]])
# uniques = unique(unlist(as.list(genes), recursive = FALSE)) # find every genes in the analysis
# 
# #Convert uniques (entrezID) to gene symbol
# uniques = bitr(uniques, fromType <- "ENTREZID", toType <- "SYMBOL", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
# head(uniques)
# uniques <- uniques$SYMBOL
# uniques
# 
# names(genes) = row.names(df3) # name each list of gene with GO associated
# 
# for (i in rownames(df3)) {
#   genes[[i]] = bitr(genes[[i]], fromType <- "ENTREZID", toType <- "SYMBOL", OrgDb <- "org.Mm.eg.db")
#   genes[[i]] <- genes[[i]]$SYMBOL
# }
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("vesicle", "recycling")
# position = unique(grep(paste(interest,collapse="|"), row.names(df3)))
# df3 = df3[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("vesicle", "recycling")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df3)))
# df3 = df3[-position, ]
# 
# data <- data.frame(matrix(NA,
#                           nrow = length(uniques),
#                           ncol = nrow(df3)), row.names = uniques)
# colnames(data) <- row.names(df3) # create empty dataframe
# 
# for (i in row.names(df3)){
#   data[i] = as.integer(as.logical(uniques %in% genes[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered3 <- data[order(apply(data, 1, FUN=sum), decreasing = TRUE), ] # reorder genes by sum, comment if clustering should be done
# 
# pheatmap(data_ordered3, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20, 
#          main = "Genes associated with GO terms", angle_col = 45, 
#          border_color = "white", legend = FALSE, filename = "ORA_KEGG_qNSCs_Heatmap.jpg", silent = TRUE)
# 
# ## Heatmap ORA/KEGG TAPs #######################################################
# df4 <- read.csv('KEGGAnalysis_LPS_ORA_TAPs.csv', row.names = 3) # number is the column you want to use as names for heatmap. Default is 3 as it is the full name of the GO term
# genes = lapply(df4$geneID, FUN = function(x) strsplit(x, "/")[[1]])
# uniques = unique(unlist(as.list(genes), recursive = FALSE)) # find every genes in the analysis
# 
# #Convert uniques (entrezID) to gene symbol
# uniques = bitr(uniques, fromType <- "ENTREZID", toType <- "SYMBOL", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
# head(uniques)
# uniques <- uniques$SYMBOL
# uniques
# 
# names(genes) = row.names(df4) # name each list of gene with GO associated
# 
# for (i in rownames(df4)) {
#   genes[[i]] = bitr(genes[[i]], fromType <- "ENTREZID", toType <- "SYMBOL", OrgDb <- "org.Mm.eg.db")
#   genes[[i]] <- genes[[i]]$SYMBOL
# }
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("vesicle", "recycling")
# position = unique(grep(paste(interest,collapse="|"), row.names(df4)))
# df4 = df4[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("vesicle", "recycling")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df4)))
# df4 = df4[-position, ]
# 
# data <- data.frame(matrix(NA,
#                           nrow = length(uniques),
#                           ncol = nrow(df4)), row.names = uniques)
# colnames(data) <- row.names(df4) # create empty dataframe
# 
# for (i in row.names(df4)){
#   data[i] = as.integer(as.logical(uniques %in% genes[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered4 <- data[order(apply(data, 1, FUN=sum), decreasing = TRUE), ] # reorder genes by sum, comment if clustering should be done
# 
# pheatmap(data_ordered4, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20, 
#          main = "Genes associated with GO terms", angle_col = 45, 
#          border_color = "white", legend = FALSE, filename = "ORA_KEGG_TAPs_Heatmap.jpg", silent = TRUE)


# ################################################################################
# #                               GSEA                                           #
# ################################################################################
# 
# #in case you work with a table present in your environment
# df = clb_markers_GSEA
# df <- df %>% rownames_to_column(var = "X")
# original_gene_list <- df$avg_log2FC # we want the log2 fold change 
# names(original_gene_list) <- df$X # name the vector
# gene_list<-na.omit(original_gene_list) # omit any NA values
# gene_list = sort(gene_list, decreasing = TRUE)
# 
# #in case you load a table
# df = read.csv("clb_markers_GSEA.csv", header=TRUE, row.names = NULL) #This will add automatically the X column
# original_gene_list <- df$avg_log2FC # we want the log2 fold change 
# names(original_gene_list) <- df$X # name the vector
# gene_list<-na.omit(original_gene_list) # omit any NA values 
# gene_list = sort(gene_list, decreasing = TRUE)


# ################################################################################
# #                               BP analysis                                    #
# ################################################################################
# 
# gse <- gseGO(geneList=gene_list, 
#              ont ="BP", 
#              keyType = "SYMBOL", 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = org.Mm.eg.db, 
#              pAdjustMethod = "BH")
# 
# head(summary(gse))
# write.csv(gse, row.names = TRUE, file = "clb_GSEA_GO.csv")
# 
# gse2 <- simplify(gse, measure="Wang", semData=NULL) 
# # The simplify method allows to remove redundant terms based on GOSemSim (select one representative term from redundant ones based on higher p.adjust, which ahve similarity higher than "cutoff)
# write.csv(gse2, row.names = T, file = "clb_GSEA_GO_simplified.csv")

# ################################################################################
# #                               Dotplot BP analysis                            #
# ################################################################################
# 
# jpeg(filename = "clb_GSEA_GO.jpg", width=800, height=500)
# p<-dotplot(gse, showCategory=10, split=".sign",orderBy = "GeneRatio", x="GeneRatio") + facet_grid(.~.sign)
# print(p)
# dev.off()
# 
# jpeg(filename = "clb_GSEA_GO_simplified.jpg", width=800, height=500)
# p<-dotplot(gse2, showCategory=20, split=".sign",orderBy = "GeneRatio", x="GeneRatio") + facet_grid(.~.sign)
# print(p)
# dev.off()

# ################################################################################
# #                               Heatmap GSEA/BP analysis                       #
# ################################################################################
# df5 <- read.csv('clb_GSEA_GO.csv', row.names = 3)
# 
# ## Heatmap GSEA/BP qNSCs ###############################################################
# df5.qNSCs <- df5 %>% filter(NES>0)
# genes2.qNSCs = lapply(df5.qNSCs$core_enrichment, FUN = function(x) strsplit(x, "/")[[1]])
# 
# uniques2.qNSCs = unique(unlist(as.list(genes2.qNSCs), recursive = FALSE)) # find every genes in the analysis
# names(genes2.qNSCs) = row.names(df5.qNSCs) # name each list of gene with GO associated
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("synapse", "synaptic")
# position = unique(grep(paste(interest,collapse="|"), row.names(df5.qNSCs)))
# df5.qNSCs = df5.qNSCs[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("assembly", "cycle", "synapse")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df5.qNSCs)))
# df5.qNSCs = df5.qNSCs[-position, ]
# 
# data2.qNSCs <- data.frame(matrix(NA,
#                            nrow = length(uniques2.qNSCs),
#                            ncol = nrow(df5.qNSCs)), row.names = uniques2.qNSCs)
# colnames(data2.qNSCs) <- row.names(df5.qNSCs) # create empty dataframe
# 
# for (i in row.names(df5.qNSCs)){
#   data2.qNSCs[i] = as.integer(as.logical(uniques2.qNSCs %in% genes2.qNSCs[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered5.up <- data2.qNSCs[order(apply(data2.qNSCs, 1, FUN=sum), decreasing = TRUE), ]
# 
# pheatmap(data_ordered5.up, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20, 
#          main = "Genes associated with GO terms", angle_col = 45, 
#          border_color = "white", legend = FALSE, filename = "GSEA_BP_qNSCs_Heatmap.jpg", silent = TRUE)

# ## Heatmap GSEA/BP TAPs ###############################################################
# df6.TAPs <- df6 %>% filter(NES<0)
# genes2.TAPs = lapply(df6.TAPs$core_enrichment, FUN = function(x) strsplit(x, "/")[[1]])
# 
# uniques2.TAPs = unique(unlist(as.list(genes2.TAPs), recursive = FALSE)) # find every genes in the analysis
# names(genes2.TAPs) = row.names(df6.TAPs) # name each list of gene with GO associated
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("synapse", "synaptic")
# position = unique(grep(paste(interest,collapse="|"), row.names(df6.TAPs)))
# df6.TAPs = df6.TAPs[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("assembly", "cycle", "synapse")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df6.TAPs)))
# df6.TAPs = df6.TAPs[-position, ]
# 
# data2.TAPs <- data.frame(matrix(NA,
#                               nrow = length(uniques2.TAPs),
#                               ncol = nrow(df6.TAPs)), row.names = uniques2.TAPs)
# colnames(data2.TAPs) <- row.names(df6.TAPs) # create empty dataframe
# 
# for (i in row.names(df6.TAPs)){
#   data2.TAPs[i] = as.integer(as.logical(uniques2.TAPs %in% genes2.TAPs[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered6.TAPs <- data2.TAPs[order(apply(data2.TAPs, 1, FUN=sum), decreasing = TRUE), ]
# 
# pheatmap(data_ordered6.TAPs, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20, 
#          main = "Genes associated with GO terms", angle_col = 45, 
#          border_color = "white", legend = FALSE, filename = "GSEA_TAPs_GO_Heatmap.jpg", silent = TRUE)

# ################################################################################
# #                               GSEA/KEGG analysis                             #
# ################################################################################
# 
# # Convert gene IDs for gseKEGG function We will lose some genes here because not all IDs will be converted
# ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),] # remove duplicate IDS
# df3 = df[df$X %in% dedup_ids$SYMBOL,] # Create a new dataframe df3 which has only the genes which were successfully mapped using the bitr function above
# df3$Y = dedup_ids$ENTREZID # Create a new column in df3 with the corresponding ENTREZ IDs
# kegg_gene_list <- df3$avg_log2FC # Create a vector of the gene universe
# names(kegg_gene_list) <- df3$Y # Name vector with ENTREZ ids
# kegg_gene_list<-na.omit(kegg_gene_list) # omit any NA values 
# kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE) # sort the list in decreasing order 
# 
# kk <- gseKEGG(geneList     = kegg_gene_list,
#               organism     = "mmu",
#               minGSSize    = 5,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               pAdjustMethod = "BH")
# write.csv(kk, row.names = TRUE, file = "KEGGAnalysis_clb_GSEA.csv")
# 
# ################################################################################
# #                               DotPlot GSEA/KEGG analysis                     #
# ################################################################################
# jpeg(filename = "KEGGAnalysis_clb_GSEA.jpg", width=800, height=500)
# p <- dotplot(kk, showCategory = 20, title = "" , split=".sign", orderBy = "GeneRatio", x="GeneRatio") + facet_grid(.~.sign)
# print(p)
# dev.off()


# ################################################################################
# #                               Heatmap GSEA/KEGG analysis                     #
# ################################################################################
# df7 <- read.csv('KEGGAnalysis_clb_GSEA.csv', row.names = 3)
# 
# # Heatmap GSEA/KEGG qNSCs ##############################################################
# df7.qNSCs <- df7 %>% filter(NES>0)
# genes2.qNSCs = lapply(df7.qNSCs$core_enrichment, FUN = function(x) strsplit(x, "/")[[1]])
# 
# uniques2.qNSCs = unique(unlist(as.list(genes2.qNSCs), recursive = FALSE)) # find every genes in the analysis
# names(genes2.qNSCs) = row.names(df7.qNSCs) # name each list of gene with GO associated
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("synapse", "synaptic")
# position = unique(grep(paste(interest,collapse="|"), row.names(df7.qNSCs)))
# df7.qNSCs = df7.qNSCs[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("assembly", "cycle", "synapse")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df7.qNSCs)))
# df7.qNSCs = df7.qNSCs[-position, ]
# 
# data2.qNSCs <- data.frame(matrix(NA,
#                               nrow = length(uniques2.qNSCs),
#                               ncol = nrow(df7.qNSCs)), row.names = uniques2.qNSCs)
# colnames(data2.qNSCs) <- row.names(df7.qNSCs) # create empty dataframe
# 
# for (i in row.names(df7.qNSCs)){
#   data2.qNSCs[i] = as.integer(as.logical(uniques2.qNSCs %in% genes2.qNSCs[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered7.up <- data2.qNSCs[order(apply(data2.qNSCs, 1, FUN=sum), decreasing = TRUE), ]
# 
# pheatmap(data_ordered7.up, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20, 
#          main = "Genes associated with GO terms", angle_col = 45, 
#          border_color = "white", legend = FALSE, filename = "GSEA_KEGG_qNSCs_Heatmap.jpg", silent = TRUE)
# 
# # Heatmap GSEA/KEGG TAPs #############################################################
# df8.TAPs <- df8 %>% filter(NES<0)
# genes2.TAPs = lapply(df8.TAPs$core_enrichment, FUN = function(x) strsplit(x, "/")[[1]])
# 
# uniques2.TAPs = unique(unlist(as.list(genes2.TAPs), recursive = FALSE)) # find every genes in the analysis
# names(genes2.TAPs) = row.names(df8.TAPs) # name each list of gene with GO associated
# 
# #Option1 - Not mandatory
# #Keep specific terms
# interest = c("synapse", "synaptic")
# position = unique(grep(paste(interest,collapse="|"), row.names(df8.TAPs)))
# df8.TAPs = df8.TAPs[position, ]
# 
# #Option2 - Not mandatory
# #Remove specific terms
# toremove = c("assembly", "cycle", "synapse")
# position = unique(grep(paste(toremove,collapse="|"), row.names(df8.TAPs)))
# df8.TAPs = df8.TAPs[-position, ]
# 
# data2.TAPs <- data.frame(matrix(NA,
#                                 nrow = length(uniques2.TAPs),
#                                 ncol = nrow(df8.TAPs)), row.names = uniques2.TAPs)
# colnames(data2.TAPs) <- row.names(df8.TAPs) # create empty dataframe
# 
# for (i in row.names(df8.TAPs)){
#   data2.TAPs[i] = as.integer(as.logical(uniques2.TAPs %in% genes2.TAPs[i][[1]])) # fill dataframe with 0 and 1 whether the gene is associated with go term or not
# }
# 
# data_ordered8.TAPs <- data2.TAPs[order(apply(data2.TAPs, 1, FUN=sum), decreasing = TRUE), ]
# 
# pheatmap(data_ordered8.TAPs, color = c("#f2f2f2", "#ff5640"), # first color is for 0 and second for 1, here a very light gray and soft red
#          cluster_rows = FALSE, treeheight_row = 20, treeheight_col = 20, 
#          main = "Genes associated with GO terms", angle_col = 45, 
#          border_color = "white", legend = FALSE, filename = "GSEA_KEGG_TAPs_Heatmap.jpg", silent = TRUE)




