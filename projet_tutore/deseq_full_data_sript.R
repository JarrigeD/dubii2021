#-----------------------------------------------------------------
# Objectif : Traiter les données de comptages entières par gène  
#            avec DESeq2 pour les normaliser avec la méthode rlog.
#
# Auteur : Domitille Jarrige
#
# Date : 2021-06-11
#-----------------------------------------------------------------


## Enregistrement du dossier de travail

setwd("/shared/ifbstor1/projects/dubii2021/djarrige/projet_scientifique/scripts")


## Chargement librairies :


required_lib <- c("FactoMineR",
                  "factoextra",
                  "knitr",
                  "tidyverse",
                  "data.table",
                  "SummarizedExperiment",
                  "tibble",
                  "DESeq2",
                  "dplyr",
                  "ClassDiscovery",
                  "pheatmap")

required_bioc <- c("DiscoRhythm",
                   "gprofiler2")

### Libraries CRAN
for (lib in required_lib) {
    if (!require(lib, character.only = TRUE)) {
        install.packages(lib)
    }
    require(lib, character.only = TRUE)
}


### Libraries Biocmanager

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

for (lib in required_bioc){
    if (!require(lib, character.only = TRUE)) {
        BiocManager::install(lib, character.only = TRUE)
    }
    require(lib, character.only = TRUE)
}

    

kable(as.data.frame(c(required_lib, required_bioc)),
      col.names = "libraries",
      caption = "Loaded required libraries")


# Essai DESeq2

## Chargement données

meta_temp <- read.csv("../data/sample_DiscoRythm_format.tsv", sep="\t")
counts_rund_df <- read.csv("../results/mapping/full_data/total_counts_int.csv", sep="\t", row.names = 1)

## Ajout de couleurs personnalisées pour les différents temps

meta_temp$color <- c(rep("#000f3b", 3), rep("#00144e", 3), rep("#001962", 3), 
                    rep("#001162", 3), rep("#001789", 3), rep("#001eb1", 3),
                    rep("#143bff", 3), rep("#86bcf9", 3), rep("#ffc576", 3), 
                    rep("#ffd400", 3), rep("#ffe900", 3), rep("#fff04e", 3),
                    rep("#fff589", 3), rep("#cec031", 3), rep("#ffbf00", 3),
                    rep("#001162", 3))

## Chargement données d'annotation

annotations <- read.csv("../data/phytozome/Creinhardtii/v5.6/annotation/gene_annotation.tsv",
                        sep="\t")
annotations_organelles <- read.csv("../data/phytozome/Creinhardtii/v5.6/annotation/organelles_annotation.tsv",
                                   sep="\t")

# /!\ certains gènes sont dupliqués dans l'annotation phytozome !!

annotations$gene_id  <- make.unique(annotations$gene_id, sep=".")


## Transformation des données de comptage en une matrice

matrix_counts_DESeq2 <- as.matrix(counts_rund_df[,2:49]) 
rownames(matrix_counts_DESeq2) <- counts_rund_df[,1]



#### Normalisation et transformation des données par DESeq2 ####
# /!\ pas pour analyse différencielle !


## Total counts par échantillon
colSums(matrix_counts_DESeq2)


# Pour la formule de design j'ai essayé avec ~ condition et ~ Time

ddsFull <- DESeqDataSetFromMatrix(countData = matrix_counts_DESeq2, 
                                  colData = meta_temp, 
                                  design = ~ time)

## Estimation des size factors

ddsFull <- estimateSizeFactors(ddsFull)

## retrait des gènes jamais exprimés

ddsFull <- ddsFull[ rowSums(counts(ddsFull)) > 0, ]


## rlog normalisation

rlog_data <- rlog(ddsFull, blind=TRUE)

## vst normalisation

vst_data <- varianceStabilizingTransformation(ddsFull, 
                                              blind = TRUE, 
                                              fitType = "parametric")



png('../doc/images/boxplot_log2_full.png', 
    width = 1500, height = 1500, res=150)
boxplot(log2(matrix_counts_DESeq2 + 0.1),
        col = meta_temp$color,
        horizontal = TRUE, 
        las = 1,
        cex.axis = 0.45,
        cex = 0.5,
        xlab = "log2(value)")
dev.off()



png('../doc/images/boxplot_rlog_full.png', 
    width = 700, height = 1000, res=150)
boxplot(assay(rlog_data),
        col = meta_temp$color,
        horizontal = TRUE, 
        las = 1,
        cex.axis = 0.45,
        cex = 0.5,
        xlab = "rlog(value)",
        main="Boxplots of rlog transformed 
gene counts")
dev.off()


png('../doc/images/boxplot_vst_full.png', 
    width = 1500, height = 1500, res=150)
boxplot(assay(vst_data),
        col = meta_temp$color,
        horizontal = TRUE, 
        las = 1,
        cex.axis = 0.45,
        cex = 0.5,
        xlab = "vst(value)")
dev.off()



## r rlog : visualisation des différentes normalisations

png(filename = "../doc/images/normalisation_comparison_full_data.png",
    width = 1300, height = 700, res=150)
par(mfrow=c(1,4))
plot(counts(ddsFull, normalized=TRUE)[,11:12],
pch=16, cex=0.3, xlim=c(0,20e3), ylim=c(0,20e3), main="normalized counts")

plot(log2(counts(ddsFull, normalized=TRUE)[,11:12] + 1),
pch=16, cex=0.3, main="log2 normalized counts")

plot(assay(rlog_data)[,11:12],
pch=16, cex=0.3, main="rlog normalized counts")

plot(assay(vst_data)[,11:12],
     pch=16, cex=0.3, main="vst normalized counts")
dev.off()


## PCA des rlog

## extract rlog matrix
counts_rlog <- assay(rlog_data)
write.csv(counts_rlog, file = "../results/mapping/full_data/counts_rlog.csv", sep="\t")

## Select most variable genes
#gvar <- apply(counts_rlog, 1, var)
#mostvargenes <- order(gvar, decreasing=TRUE)[1:1000]

## Run PCA

res_pca <- PCA(t(counts_rlog), graph=FALSE)

## Graphes de l'ACP

png('../doc/images/eigen_full.png', 
    width = 1000, height = 700, res=150)
fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()


png('../doc/images/indiv_full.png', 
    width = 1000, height = 700, res=200)
fviz_pca_ind(res_pca, label="none",
             habillage=as.factor(meta_temp$time),
             mean.point=FALSE,
             pointshape = 19) 
dev.off()


## Clustering hierarchique

dists <- dist(t(assay(rlog_data)))
par(mfrow=c(1,1))
tree_rlog <- hclust(dists)
plot(tree_rlog)

dists_vst <- dist(t(assay(vst_data)))
tree_vst <- hclust(dists_vst)
#plot(tree_vst)

my_colors = c("#ffc576", "#445cd7", "#fff04e", "#001789", "#faf5c6")

png(filename = "../doc/images/clustering_rlog_euclidian.png", 
    width = 1000, height = 700, res=150)
par(bg = "darkgrey", mfrow=c(1, 1))
plotColoredClusters(tree_rlog, labs = meta_temp$time,
                    ylab = NA, xlab = NA, cex = , las = 1,
                    cols = meta_temp$color, col = "white",
                    main = "Samples Euclidian distance hierarchical 
clustering, rlog, complete linkage")
rect.hclust(tree_rlog, k=5, border=my_colors)
dev.off()

plotColoredClusters(tree_rlog, labs = meta_temp$time,
                    ylab = NA, xlab = NA, cex = 1 , las = 1,
                    cols = meta_temp$color, col = "white",
                    main = "Samples Euclidian distance hierarchical 
clustering, VST, complete linkage")

rect.hclust(tree_rlog, k=5, border=my_colors)

par(bg="white")


# Enregistrement des résultats de clustering rlog

clusters_names <- c("early_night",
                    "night",
                    "dawn",
                    "day",
                    "dusk")
clusters_rlog <- cutree(tree_rlog, k=5)
meta_temp$cluster <- clusters_names[clusters_rlog]


#### Analyse différentielle avec DESeq2 selon les clusters ####

ddsClusters <- DESeqDataSetFromMatrix(countData = matrix_counts_DESeq2, 
                                      colData = meta_temp, 
                                      design = ~ cluster) # essai avec time ou cluster ?

## retrait des gènes jamais exprimés

ddsClusters <- ddsClusters[ rowSums(counts(ddsClusters)) > 60, ]
write.csv(counts(ddsClusters), file="../results/mapping/full_data/full_data_counts_for_SARTools_nettoye.csv", sep="\t")

ddsClusters <- DESeq(ddsClusters)

resDESeq <- results(ddsClusters)
print(resDESeq)
plotMA(resDESeq, ylim = c(-1, 1))

par(bg="white", mfrow=c(2, 2))


for (clu in clusters_names){
    for (ster in clusters_names){
        if (clu != ster){
            text <- paste0(clu, "/", ster)
            resDESeq <- results(ddsClusters, contrast = c("cluster", clu, ster),
                                independentFiltering = TRUE, alpha=0.01)
            print(text)
            print( sum( resDESeq$padj < 0.01, na.rm=TRUE ) )
            plotMA(resDESeq, ylim = c(-9, 9), main = text)
        }
    }
}


### Essai sur les temps contigus seulement 

pair_wise_cluster = list(c("early_night", "night"),
                         c("night", "dawn"),
                         c("dawn", "day"),
                         c("day", "dusk"),
                         c("dusk", "early_night"))

# initialisation de significant_genes
significant_genes <- c()

for (ele in pair_wise_cluster){
    print(ele[1])
    resDESeq <- results(ddsClusters, contrast = c("cluster", ele[1], ele[2]),
                        independentFiltering = TRUE, alpha=0.01)
    temp_list <- row.names(resDESeq[which(resDESeq$padj < 0.01),])
    print((length(temp_list)/17823 * 100))
    significant_genes <- union(temp_list, significant_genes)

}


#file_name <- paste0("../doc/images/DESeq_", clu, ".png")
#print(file_name)
#png(filename = file_name, width = 700, height = 700, res=150)

message("Gènes différentiellement exprimés :")
print((length(significant_genes)/17823 * 100))

sign_df <- data.frame(row.names = significant_genes)
sign_df$diff_expr <- "yes"
sign_df$gene_id <- row.names(sign_df)






#### Analyse rhythmique ####


# Essai DiscoRhythm

# Je convertis la matrice des rlog en dataframe...
gene_id <- row.names(counts_rlog)
rlog_df <- as.data.frame(counts_rlog)
rlog_df <- add_column(rlog_df, gene_id, .before=1)
#write_csv(rlog_df, file = "../rlog_full.csv")

# Puis en SE (Submarised Experiment)
input_data <- discoDFtoSE(rlog_df)

discoDesignSummary(input_data)

# Ici les données ne sont pas équidistantes on ne peut utiliser
# que cosinor ou Lomb-Scargle

# Détection de rythmes avec la méthode cosinor
rythms_genes_CS <- discoODAs(input_data, period = 24, method = "CS",
                             circular_t = TRUE, ncores = 4)

# Détection de rythmes avec la méthode Lomb-Scargle
rythms_genes_LS <- discoODAs(input_data, period = 24, method = "LS", 
                                circular_t = TRUE, ncores = 4)

png(filename = "../doc/images/qvalues_cs_ls.png",
    width = 1300, height = 700, res=150)
par(mfrow=c(1,2))
hist(data.frame(rythms_genes_CS)$CS.qvalue, breaks = 100, 
     main = "qvalue cosinor",
     xlab = "qvalue")
hist(data.frame(rythms_genes_LS)$LS.qvalue, breaks = 100, 
     main = "qvalue Lomb-Scargle",
     xlab = "qvalue")
dev.off()

hist(data.frame(rythms_genes_CS)$CS.acrophase, breaks = 100, 
     main = "acrophases cosinor",
     xlab = "hour")
hist(data.frame(rythms_genes_LS)$LS.acrophase, breaks = 100, 
     main = "acrophases Lomb-Scargle",
     xlab = "hour")



# Extracting just the dataframes from the discorythm output
CS_df <- data.frame(rythms_genes_CS)
LS_df <- data.frame(rythms_genes_LS) 


# Retreiving only most significant genes 
# qvalue inferior to 0.05
CS_df <- CS_df[CS_df$CS.qvalue < 0.05,]
LS_df <- LS_df[LS_df$LS.qvalue < 0.05,]

#png(filename = "../doc/images/qvalues_cs_ls.png",
#    width = 1300, height = 700, res=150)
par(mfrow=c(1,2))
hist(CS_df$CS.qvalue, breaks = 100, 
     main = "qvalues cosinor 
after selection",
     xlab = "qvalue")
hist(LS_df$LS.qvalue, breaks = 100, 
     main = "qvalues Lomb-Scargle 
after selection",
     xlab = "qvalue")
#dev.off()

hist(CS_df$CS.acrophase, breaks = 100, 
     main = "acrophases cosinor 
after selection",
     xlab = "hour")
hist(LS_df$LS.acrophase, breaks = 100, 
     main = "acrophases Lomb-Scargle 
after selection",
     xlab = "hour")


## Extraction des annotations sur les OPR

row_CS <- row.names(CS_df)
row_LS <- row.names(LS_df)

CS_df$gene_id <- row_CS
LS_df$gene_id <- row_LS

CS_df <- left_join(CS_df, annotations)
LS_df <- left_join(LS_df, annotations)

id_organelle = annotations_organelles$gene_id
row.names(annotations_organelles) <- id_organelle

row_CS -> row.names(CS_df)
row_LS -> row.names(LS_df)

for (id in id_organelle) {
    CS_df[id, "gene_symbol"] <- annotations_organelles[id,"gene_symbol"]
    LS_df[id, "gene_symbol"] <- annotations_organelles[id,"gene_symbol"]
    
}

# Je continue avec les données cosinor

CS_df$gene_id <- row.names(CS_df)

CS_df$encoded <- "Nucleus"

CS_df[which(startsWith(CS_df$gene_id, "CreMt.")), "encoded"] <- "Mitochondrion"
CS_df[which(startsWith(CS_df$gene_id, "CreCp.")), "encoded"] <- "Chloroplast"
CS_df[which(startsWith(CS_df$gene_id, "CreMt.")), "subcellular_location"] <- "Mitochondrion"
CS_df[which(startsWith(CS_df$gene_id, "CreCp.")), "subcellular_location"] <- "Chloroplast"
CS_df[which(startsWith(CS_df$gene_id, "CreMt.")), "simplified_subcellular_location"] <- "Mitochondrion"
CS_df[which(startsWith(CS_df$gene_id, "CreCp.")), "simplified_subcellular_location"] <- "Chloroplast"

CS_df[which(startsWith(CS_df$gene_description, "OctotricoPeptide") == TRUE), ] -> CS_subset_opr
LS_df[which(startsWith(LS_df$gene_description, "OctotricoPeptide") == TRUE), ] -> LS_subset_opr
CS_df[which(startsWith(CS_df$gene_description, "PentatricoPeptide") == TRUE), ] -> CS_subset_ppr
LS_df[which(startsWith(LS_df$gene_description, "PentatricoPeptide") == TRUE), ] -> LS_subset_ppr

final_data <- left_join(CS_df, sign_df, by = "gene_id")
final_data <- final_data[which(final_data$diff_expr == "yes"),]
final_data <- final_data[ which(!is.na(final_data$CS.acrophase)) , ]



#### Analyse générale des gènes rhythmiques ####

par(mfrow=c(1,1))
#plot(tree_gene, label=FALSE)


ggplot(data = final_data, mapping = aes(x=CS.acrophase, fill=encoded))+
    #geom_histogram(binwidth = 1)
    geom_density() +
    facet_wrap(~encoded) +
    scale_fill_manual(values = c("#99e55c", "#eb957c", "#0fc6e1"))



c(paste0("Chloroplast", " (", nrow(final_data[final_data$simplified_subcellular_location == "Chloroplast",]), ")"),
  paste0("Chromosome", " (", nrow(final_data[final_data$simplified_subcellular_location == "Chromosome",]), ")"),
  paste0("Cilium", " (", nrow(final_data[final_data$simplified_subcellular_location == "Cilium",]), ")"),
  paste0("Cytoplasm", " (", nrow(final_data[final_data$simplified_subcellular_location == "Cytoplasm",]), ")"),
  paste0("Cytoskeleton", " (", nrow(final_data[final_data$simplified_subcellular_location == "Cytoskeleton",]), ")"),
  paste0("Endoplasmic reticulum", " (", nrow(final_data[final_data$simplified_subcellular_location == "Endoplasmic reticulum",]), ")"),
  paste0("Golgi apparatus", " (", nrow(final_data[final_data$simplified_subcellular_location == "Golgi apparatus",]), ")"),   
  paste0("Membrane", " (", nrow(final_data[final_data$simplified_subcellular_location == "Membrane",]), ")"),
  paste0("Mitochondrion", " (", nrow(final_data[final_data$simplified_subcellular_location == "Mitochondrion",]), ")"),
  paste0("Nucleus", " (", nrow(final_data[final_data$simplified_subcellular_location == "Nucleus",]), ")"),
  paste0("Other", " (", nrow(final_data[final_data$simplified_subcellular_location == "Other",]), ")"),
  paste0("unknown", " (", nrow(final_data[final_data$simplified_subcellular_location == "unknown",]), ")")) -> locations
   

png(filename = "../doc/images/acrophase_locations.png", width=1400, height=900, res=150)
ggplot(data = final_data, mapping = aes(x=CS.acrophase, fill=simplified_subcellular_location))+
    annotate("rect", fill = "darkgray", xmin = 0, xmax = 11, ymin =0, ymax = Inf, alpha =0.3) +
    annotate("rect", fill = "gold", xmin = 11, xmax = 23, ymin =0, ymax = Inf, alpha =0.3) +
    annotate("rect", fill = "darkgray", xmin = 23, xmax = Inf, ymin =0, ymax = Inf, alpha =0.3)+
    geom_density() +
    facet_wrap(~simplified_subcellular_location)+
    labs( x = "Acrophase (heures)", y = "Densité (% gènes)",
            title ="Répartition des acrophases d'après la localisation cellulaire",
            subtitle = "Méthode cosinor")+
    scale_fill_discrete(name = "Localisation cellulaire", labels = locations)
dev.off()

png(filename="../doc/images/acrophases_general.png", 
    width=1000, height=700, res=150)
ggplot(data = final_data, mapping = aes(x=CS.acrophase)) +
    annotate("rect", fill = "darkgray", xmin = 0, xmax = 11, ymin =0, ymax = Inf, alpha =0.3) +
    annotate("rect", fill = "gold", xmin = 11, xmax = 23, ymin =0, ymax = Inf, alpha =0.3) +
    annotate("rect", fill = "darkgray", xmin = 23, xmax = Inf, ymin =0, ymax = Inf, alpha =0.3) +
    geom_density() +
    labs( x = "Acrophase (heures)", y = "Densité (% gènes)",
          title ="Répartition des acrophases",
          subtitle = "Gènes rhythmiques (Méthode cosinor) et différentiellement exprimés")
    #scale_x_continuous(limits = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24))
dev.off()
    
    

ggplot(data = final_data, mapping = aes(x=CS.acrophase, y=simplified_subcellular_location))+
    #geom_histogram(binwidth = 1)
    geom_count()
    #facet_wrap(~simplified_subcellular_location)



#### Etude des OTAFs ####

# Récupération des OTAFs rhytmiques :
opr <- final_data[ which(startsWith(final_data$gene_description, "OctotricoPeptide Repeat")) , ]                                                           
ppr <- final_data[ which(startsWith(final_data$gene_description, "PentatricoPeptide Repeat")) , ]
tpr <- final_data[ which(startsWith(final_data$gene_description, "TetratricoPeptide Repeat")) , ]

otaf <- union_all(opr, ppr)
otaf <- union_all(otaf, tpr)

# Ajout des données d'expression :
columns <- colnames(otaf)
otaf <- left_join(otaf, rlog_df)
colnames(otaf) <- make.unique(c(columns, rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)))



lab <- c(paste0("Chloroplaste (", nrow(otaf[which(otaf$subcellular_location == "Chloroplast"),]), ")"),
         paste0("Mitochondrie (", nrow(otaf[which(otaf$simplified_subcellular_location == "Mitochondrion"),]), ")"),
         paste0("Inconnu (", nrow(otaf[which(otaf$simplified_subcellular_location == "unknown"),]), ")"))
ggplot(data = otaf, mapping = aes(x=CS.acrophase, fill=simplified_subcellular_location))+
    geom_density() +
    #geom_histogram() +
    facet_wrap(~simplified_subcellular_location)+
    labs( x = "Acrophase (h)", y = "Densité (% gènes)",
          title ="Répartition des acrophases d'OTAF",
          subtitle = "d'après la localisation cellulaire, Méthode cosinor")+
    scale_fill_manual(name = "Localisation cellulaire", 
                      values = c("#99e55c", "#eb957c", "#b1aeae"),
                      labels = lab)

##### Inférences de régulations #####

## Extraction des informations par OTAF : ##

# OTAF prédites importées dans le chloroplaste
multi_otaf_chloro <- data.frame(time=0, gene_expression=0.0, gene="")

for (prot in otaf[which(otaf$subcellular_location == "Chloroplast"),]$gene_symbol){
    otaf[otaf$gene_symbol == prot, 19:66] %>%
        pivot_longer(cols = c(1:48), names_to = "time", values_to = "gene_expression") -> tmp_df
    tmp_df$time <- rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)
    tmp_df$gene <- prot
    multi_otaf_chloro <- rbind(multi_otaf_chloro, tmp_df)
}

multi_otaf_chloro <- multi_otaf_chloro[-1,]



# OTAF prédites importées dans la mitochondrie
multi_otaf_mito <- data.frame(time=0, gene_expression=0.0, gene="")

for (prot in otaf[which(otaf$subcellular_location == "Mitochondrion"),]$gene_symbol){
    otaf[otaf$gene_symbol == prot, 19:66] %>%
        pivot_longer(cols = c(1:48), names_to = "time", values_to = "gene_expression") -> tmp_df
    tmp_df$time <- rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)
    tmp_df$gene <- prot
    multi_otaf_mito <- rbind(multi_otaf_mito, tmp_df)
}

multi_otaf_mito <- multi_otaf_mito[-1,]



## Extraction des informations des organelles : ##

chloro <- final_data[ which(startsWith(final_data$gene_id, "CreCp")) , ]
mito <- final_data[ which(startsWith(final_data$gene_id, "CreMt")) , ]


# Ajout des données d'expression :
columns <- colnames(chloro)
chloro <- left_join(chloro, rlog_df)
colnames(chloro) <- make.unique(c(columns, rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)))
chloro$gene_symbol <- make.unique(chloro$gene_symbol)

columns <- colnames(mito)
mito <- left_join(mito, rlog_df)
colnames(mito) <- make.unique(c(columns, rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)))


# Chloroplaste
multi_chloro <- data.frame(time=0, gene_expression=0.0, gene="")

for (prot in chloro$gene_symbol){
    chloro[chloro$gene_symbol == prot, 19:66] %>%
        pivot_longer(cols = c(1:48), names_to = "time", values_to = "gene_expression") -> tmp_df
    tmp_df$time <- rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)
    tmp_df$gene <- prot
    multi_chloro <- rbind(multi_chloro, tmp_df)
}

multi_chloro <- multi_chloro[-1,]


# Mitochondrie
multi_mito <- data.frame(time=0, gene_expression=0.0, gene="")

for (prot in mito$gene_symbol){
    mito[mito$gene_symbol == prot, 19:66] %>%
        pivot_longer(cols = c(1:48), names_to = "time", values_to = "gene_expression") -> tmp_df
    tmp_df$time <- rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)
    tmp_df$gene <- prot
    multi_mito <- rbind(multi_mito, tmp_df)
}

multi_mito <- multi_mito[-1,]



## Plots ! ##

periode=24

png(filename = "../doc/images/otaf_chloro.png", width = 500, height = 500)
ggplot(data = multi_otaf_chloro, mapping = aes(x=time, y=gene_expression, group=gene, color=gene)) + 
    #facet_wrap(~gene) +
    geom_smooth(method = "lm", se = FALSE, level = 0.95,
                formula = y ~ sin(x / periode * 2 * pi) + cos(x / periode * 2 * pi),
                fullrange = TRUE) +
    theme(legend.position = "none") +
    geom_point() +
    labs(y = "rlog(counts)", x = "Temps (h)",
         title ="Modèles d'expression OTAF chloroplastiques",
         subtitle = "Principalement prédiction de localisation, Méthode cosinor")
dev.off()


png(filename = "../doc/images/ARN_chloro.png", width = 500, height = 500)  
ggplot(data = multi_chloro, mapping = aes(x=time, y=gene_expression, group=gene, color=gene)) + 
    #facet_wrap(~gene) +
    geom_smooth(method = "lm", se = FALSE, level = 0.95,
                formula = y ~ sin(x / periode * 2 * pi) + cos(x / periode * 2 * pi),
                fullrange = TRUE) +
    theme(legend.position = "none") +
    geom_point() +
    labs(y = "rlog(counts)", x = "Temps (h)",
         title ="Modèles d'expression gènes chloroplastiques",
         subtitle = "Méthode cosinor")
dev.off()


png(filename = "../doc/images/otaf_mito.png", width = 500, height = 500)
ggplot(data = multi_otaf_mito, mapping = aes(x=time, y=gene_expression, group=gene, color=gene)) + 
    #facet_wrap(~gene) +
    geom_smooth(method = "lm", se = FALSE, level = 0.95,
                formula = y ~ sin(x / periode * 2 * pi) + cos(x / periode * 2 * pi),
                fullrange = TRUE) +
    theme(legend.position = "none") +
    geom_point() +
    labs(y = "rlog(counts)", x = "Temps (h)",
         title ="Modèles d'expression OTAF mitochondriaux",
         subtitle = "Principalement prédiction de localisation, Méthode cosinor")
dev.off()


png(filename = "../doc/images/RNA_mito.png", width = 500, height = 500)
ggplot(data = multi_mito, mapping = aes(x=time, y=gene_expression, group=gene, color=gene)) + 
    #facet_wrap(~gene) +
    geom_smooth(method = "lm", se = FALSE, level = 0.95,
                formula = y ~ sin(x / periode * 2 * pi) + cos(x / periode * 2 * pi),
                fullrange = TRUE) +
    theme(legend.position = "none") +
    geom_point() +
    labs(y = "rlog(counts)", x = "Temps (h)",
         title ="Modèles d'expression gènes mitochondriaux",
         subtitle = "Méthode cosinor")
dev.off()

########################################
## Correlation of expression patterns ##
########################################


acro_day_chloro <- data.frame(time=0, gene_expression=0.0, gene="")

for (prot in chloro[which(chloro$CS.acrophase <14 & chloro$CS.acrophase >12 ),]$gene_symbol) {
    if (!startsWith(prot, "I-")) {
        chloro[chloro$gene_symbol == prot, 19:66] %>%
            pivot_longer(cols = c(1:48), names_to = "time", values_to = "gene_expression") -> tmp_df
        tmp_df$time <- rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)
        tmp_df$gene <- prot
        acro_day_chloro <- rbind(acro_day_chloro, tmp_df)
    }
}

acro_day_chloro <- acro_day_chloro[-1,]


otaf_chloro_acro_day <- data.frame(time=0, gene_expression=0.0, gene="")

for (prot in otaf[which(otaf$subcellular_location == "Chloroplast" & otaf$CS.acrophase <10.6 & otaf$CS.acrophase >8),]$gene_symbol){
    otaf[otaf$gene_symbol == prot, 19:66] %>%
        pivot_longer(cols = c(1:48), names_to = "time", values_to = "gene_expression") -> tmp_df
    tmp_df$time <- rep(c(0, 2, 4, 6, 8, 10, 10.5, 11, 11.5, 12, 14, 16, 18, 20, 22, 24), each=3)
    tmp_df$gene <- prot
    otaf_chloro_acro_day <- rbind(otaf_chloro_acro_day, tmp_df)
}

otaf_chloro_acro_day <- otaf_chloro_acro_day[-1,]




png(filename = "../doc/images/otaf_regu.png", width = 700, height = 700)
ggplot(data = otaf_chloro_acro_day, mapping = aes(x=time, y=gene_expression, group=gene, color=gene)) + 
    annotate("rect", fill = "darkgray", xmin = 0, xmax = 11, ymin =6, ymax = 13, alpha =0.3) +
    annotate("rect", fill = "gold", xmin = 11, xmax = 23, ymin =6, ymax = 13, alpha =0.3) +
    annotate("rect", fill = "darkgray", xmin = 23, xmax = Inf, ymin =6, ymax = 13, alpha =0.3)+
    facet_wrap(~gene) +
    ylim(6, 13) +
    geom_smooth(method = "lm", se = FALSE, level = 0.95,
                formula = y ~ sin(x / periode * 2 * pi) + cos(x / periode * 2 * pi),
                fullrange = TRUE) +
    theme(legend.position = "none") +
    geom_point() +
    labs(y = "rlog(counts)", x = "Temps (h)",
         title ="OTAF chloroplastiques, fin de nuit",
         subtitle = "Méthode cosinor") +
    theme(plot.title = element_text(size = 20))+
    theme(strip.text.x = element_text(size = 20))
dev.off()

png(filename = "../doc/images/otaf_target.png", width = 700, height = 700)
ggplot(data = acro_day_chloro, mapping = aes(x=time, y=gene_expression, group=gene, color=gene)) + 
    ylim(10, 19.5) +
    annotate("rect", fill = "darkgray", xmin = 0, xmax = 11, ymin =10, ymax = 19.5, alpha =0.3) +
    annotate("rect", fill = "gold", xmin = 11, xmax = 23, ymin =10, ymax = 19.5, alpha =0.3) +
    annotate("rect", fill = "darkgray", xmin = 23, xmax = Inf, ymin =10, ymax = 19.5, alpha =0.3)+
    facet_wrap(~gene) +
    geom_smooth(method = "lm", se = FALSE, level = 0.95,
                formula = y ~ sin(x / periode * 2 * pi) + cos(x / periode * 2 * pi),
                fullrange = TRUE) +
    theme(legend.position = "none") +
    geom_point() +
    labs(y = "rlog(counts)", x = "Temps (h)",
         title ="Gènes chloroplastiques, aube",
         subtitle = "Méthode cosinor")+
    theme(plot.title = element_text(size = 20))+
    theme(strip.text.x = element_text(size = 20))
dev.off()




#####################
#### Protéomique ####
#####################


# Je charge les données d'abondance protéique
proteo_df <- read.csv("../data/proteomics/proteo.tsv", sep="\t", row.names = 1)

prot_stat <- data.frame(row.names=row.names(proteo_df))
prot_stat$zero <- apply(proteo_df == 0, 1, sum, na.rm = TRUE)
undetected_prot <- prot_stat$zero >= 1

proteo_df <- proteo_df[!undetected_prot,]
#proteo_df <- na.omit(proteo_df)

proteo_meta <- read.csv("../data/protein_meta.csv")

# Transfo en matrice
matrix_proteo <- as.matrix(proteo_df)
#matrix_proteo <- round(matrix_proteo)
#matrix_proteo <- gsub(",", "\t", matrix_proteo)    
#matrix_proteo <- as.integer(matrix_proteo)
#matrix_proteo <- na.omit(matrix_proteo)

# Transformation en log2

proteo_log2_df <- log2(proteo_df)

png(filename = "../doc/images/boxplot_prot_log2.png", 
    width = 700, height = 1000, res=150)
boxplot(proteo_log2_df, col = proteo_meta$color,
        horizontal = TRUE, 
        las = 1,
        cex.axis = 0.45,
        cex = 0.5,
        xlab = "log2(value)",
        main="Boxplots of log2 transformed 
protein abundance")
dev.off()


# ACP
pca_prot <- PCA(t(proteo_log2_df), scale.unit = TRUE, graph = FALSE)
fviz_eig(pca_prot)
png(filename = "../doc/images/pca_prot.png", width = 1000, height = 700, res=200)
fviz_pca_ind(pca_prot, habillage = as.factor(proteo_meta$time),
             mean.point = FALSE)
dev.off()


# Clustering hierarchique

dist_prot_pearson <- as.dist(1 - cor(proteo_log2_df, use = "everything", method = "pearson"))
tree_prot <- hclust(dist_prot_pearson, method = "average")
clusters_prot <- cutree(tree_prot, k=4) 

clusters_names <- c("dawn",
                    "night",
                    "dusk",
                    "day")
proteo_meta$cluster <- clusters_names[clusters_prot]

png(filename = "../doc/images/clustering_prot_pearson.png", 
    width = 1000, height = 700, res=150)
par(bg = "darkgrey", mfrow=c(1, 1))
plotColoredClusters(tree_prot, labs = proteo_meta$time,
                    ylab = NA, xlab = NA, cex = , las = 1,
                    cols = proteo_meta$color, col = "white",
                    main = "Samples Pearson distance hierarchical 
clustering, proteins")
rect.hclust(tree_prot, k=4, 
            border = c("#faf5c6", "#001789", "#ffc576", "#fff04e"))
dev.off()
par(bg = "white", mfrow=c(1, 1))

### Analyse differentielle DESeq2 ###
# NE MARCHE PAS À CAUSE DE LA CONVERTION EN INTEGER !!! >:-c

#ddsProt <- DESeqDataSetFromMatrix(countData = matrix_proteo, 
#                                  colData = proteo_meta, 
#                                  design = ~ cluster) # essai avec time ou cluster ?

#ddsProt <- DESeq(ddsProt)


### Essai sur les clusters contigus seulement 

#pair_wise_cluster = list(c("night", "dawn"),
#                         c("dawn", "day"),
#                         c("day", "dusk"),
#                         c("dusk", "night"))

# initialisation de significant_prot
#significant_prot <- c()

#for (ele in pair_wise_cluster){
#    print(ele[1])
#    resDESeq <- results(ddsProt, contrast = c("cluster", ele[1], ele[2]),
#                        independentFiltering = TRUE, alpha=0.01)
#    plotMA(resDESeq, ylim = c(-10, 10))
#    temp_list <- row.names(resDESeq[which(resDESeq$padj < 0.01),])
#    print((length(temp_list)/17823 * 100))
#    significant_prot <- union(temp_list, significant_prot)
#}

#message("Protéines différentiellement exprimés :")
#print((length(significant_prot)/1930 * 100))

#sign_prot_df <- data.frame(row.names = significant_prot)
#sign_prot_df$diff_expr <- "yes"
#sign_prot_df$gene_id <- row.names(sign_prot_df)


### Essai DiscoRhythm ###

proteo_log2_df <- add_column(.data = proteo_log2_df, .before = 1,
                       gene_id = row.names(proteo_df))
# Puis en SE (Submarised Experiment)
proteo_input_data <- discoDFtoSE(proteo_log2_df)

discoDesignSummary(proteo_input_data)

# Ici les données ne sont pas équidistantes on ne peut utiliser
# que cosinor ou Lomb-Scargle

# Détection de rythmes avec la méthode cosinor
rythms_prot_CS <- discoODAs(proteo_input_data, period = 24, method = "CS",
                             circular_t = FALSE, ncores = 4)

# Détection de rythmes avec la méthode Lomb-Scargle
rythms_prot_LS <- discoODAs(proteo_input_data, period = 24, method = "LS", 
                             circular_t = FALSE, ncores = 4)

# Enregistrement en dataframe

prot_CS_df <- data.frame(rythms_prot_CS)
prot_LS_df <- data.frame(rythms_prot_LS)

par(mfrow=c(1,2))
hist(prot_CS_df$CS.qvalue, breaks = 100)
hist(prot_LS_df$LS.qvalue, breaks = 100) # LS est nulle ici !

prot_CS_df <- prot_CS_df[ which(prot_CS_df$CS.qvalue < 0.05),]
prot_LS_df <- prot_LS_df[ which(prot_CS_df$LS.qvalue < 0.05),] # LS est nulle ici !


par(mfrow=c(1,1))
hist(prot_CS_df$CS.qvalue, breaks = 100, xlim=c(0,1))
#hist(prot_LS_df$LS.qvalue, breaks = 100) # LS est nulle ici !

# Analyse rhythmique peu concluante ! Seulement 82 protéines rhytmiques détectées.

prot_id <- row.names(prot_CS_df)

prot_CS_df <- add_column(.data = prot_CS_df, gene_id = prot_id, .before = 1)

prot_CS_df <- left_join(prot_CS_df, annotations)

row.names(prot_CS_df) <- prot_id

for (id in (intersect(prot_id, id_organelle))) {
    prot_CS_df[id, "gene_symbol"] <- annotations_organelles[id,"gene_symbol"]
}

ggplot(data=prot_CS_df, mapping = aes(x=CS.acrophase))+
    geom_histogram()


#### Clustering hierarchique ####

# Je ne garde que les gènes rhythmiques et DE
to_keep <- intersect(row.names(counts_rlog), final_data$gene_id)
counts_rlog_filt <- counts_rlog[to_keep,]

#distance_gene <- dist(counts_rlog_filt)
#tree_gene <- hclust(distance_gene, method = "complete")

dist_pearson <- as.dist(1 - cor(t(counts_rlog_filt), use = "everything", method = "pearson"))
tree_gene_pearson <- hclust(dist_pearson, method = "complete")

plot(tree_gene_pearson, lab = FALSE, hang=-1)
clust <- cutree(tree_gene_pearson, k = 4)
rect.hclust(tree_gene_pearson, k=4, border=rainbow(4))
clust_names <- c("clust1", "clust2", "clust3", "clust4")
final_data$cluster <- clust_names[clust]

annot_clust <- data.frame(as.factor(final_data$cluster))
colnames(annot_clust) <- "cluster"
row.names(annot_clust) <- final_data$gene_id

annot_sample <- data.frame(as.factor(meta_temp$cluster))
colnames(annot_sample) <- "periode"
row.names(annot_sample) <- colnames(counts_rlog_filt)

ann_colors = list(periode = c(dawn="#faf5c6", day="#fff04e", dusk="#ffc576", early_night="#445cd7", night="#001789"),
                  cluster = c(clust1="cyan", clust2="darkgreen", clust3="red", clust4="gold"))

pheatmap(counts_rlog_filt,
         filename = "../doc/images/biclustering.png",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = dist_pearson,
         clustering_distance_cols = dists,
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = annot_clust,
         annotation_names_row = FALSE,
         annotation_col = annot_sample,
         annotation_names_col = FALSE,
         annotation_colors = ann_colors,
         cutree_rows = 4,
         cutree_cols = 5,
         scale = "row",
         use_raster = TRUE,
         #angle_col = "45",
         main = "Biclustering gènes et échantillons")


#### Analyse enrichissement fonctionnel ####

list_cluster = list(clust1 = final_data[which(final_data$cluster == "clust1"), "gene_ontology_id"],
                    clust2 = final_data[which(final_data$cluster == "clust2"), "gene_ontology_id"],
                    clust3 = final_data[which(final_data$cluster == "clust3"), "gene_ontology_id"],
                    clust4 = final_data[which(final_data$cluster == "clust4"), "gene_ontology_id"])

gost(query = final_data$gene_ontology_id,
     organism = "creinhardtii",
     significant = TRUE,
     multi_query = FALSE,) -> gostres

gost(query = list_cluster[2],
     organism = "creinhardtii",
     significant = TRUE,
     multi_query = FALSE) -> gostres
gostplot(gostres, capped = FALSE, interactive = TRUE)

sessionInfo()





