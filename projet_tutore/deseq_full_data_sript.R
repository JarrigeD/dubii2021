#-----------------------------------------------------------------
# Objectif : Traiter les données de comptages entières par gène  
#            avec DESeq2 pour les normaliser avec la méthode rlog.
#
# Auteur : Domitille Jarrige
#
# Date : 2021-06-02
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
                  "ClassDiscovery")

required_bioc <- c("DiscoRhythm",
                   "gProfiler2")

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
    width = 1500, height = 1500, res=150)
boxplot(assay(rlog_data),
        col = meta_temp$color,
        horizontal = TRUE, 
        las = 1,
        cex.axis = 0.45,
        cex = 0.5,
        xlab = "rlog(value)")
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
    width = 1000, height = 700, res=150)
fviz_pca_ind(res_pca, label="none",
             habillage=as.factor(meta_temp$time),
             mean.point=FALSE,
             pointshape = 19)
dev.off()


## Clustering hierarchique

dists <- dist(t(assay(rlog_data)))
par(mfrow=c(1,1))
tree_rlog <- hclust(dists)
#plot(tree_rlog)

dists_vst <- dist(t(assay(vst_data)))
tree_vst <- hclust(dists_vst)
#plot(tree_vst)

my_colors = c("blue", "gold", "red", "darkgreen", "purple")

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
                                      design = ~ cluster)

## retrait des gènes jamais exprimés

ddsClusters <- ddsClusters[ rowSums(counts(ddsClusters)) > 0, ]

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
                                independentFiltering = TRUE, alpha=0.05)
            print(text)
            print( sum( resDESeq$padj < 0.05, na.rm=TRUE ) )
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
                        independentFiltering = TRUE, alpha=0.05)
    temp_list <- row.names(resDESeq[which(resDESeq$padj < 0.05),])
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


par(mfrow=c(1,2))
hist(data.frame(rythms_genes_CS)$CS.qvalue, breaks = 100, 
     main = "qvalue cosinor",
     xlab = "qvalue")
hist(data.frame(rythms_genes_LS)$LS.qvalue, breaks = 100, 
     main = "qvalue Lomb-Scargle",
     xlab = "qvalue")

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

png(filename = "../doc/images/qvalues_cs_ls.png",
    width = 1300, height = 700, res=150)
par(mfrow=c(1,2))
hist(CS_df$CS.qvalue, breaks = 100, 
     main = "qvalues cosinor 
after selection",
     xlab = "qvalue")
hist(LS_df$LS.qvalue, breaks = 100, 
     main = "qvalues Lomb-Scargle 
after selection",
     xlab = "qvalue")
dev.off()

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


#### Essai analyse des gènes + acrophase

#distance_gene <- dist(matrix_counts_DESeq2)
#tree_gene <- hclust(distance_gene, method = "complete")

par(mfrow=c(1,1))
#plot(tree_gene, label=FALSE)




ggplot(data = final_data, mapping = aes(x=CS.acrophase, fill=encoded))+
    #geom_histogram(binwidth = 1)
    geom_density() +
    facet_wrap(~encoded) +
    scale_fill_manual(values = c("#99e55c", "#eb957c", "#0fc6e1"))

ggplot(data = final_data, mapping = aes(x=CS.acrophase, fill=simplified_subcellular_location))+
    #geom_histogram(binwidth = 1)
    geom_density() +
    facet_wrap(~simplified_subcellular_location)


#### Analyse enrichissement fonctionnel :

final_data[]
list_cluster = list()

gost(query = ,
     organism = "creinhardtii",
     significant = TRUE,
     multi_query = TRUE,)

sessionInfo()





