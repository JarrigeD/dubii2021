######################################################################################################
# Script purpose : Takes a time series counts table file and produce a list of rhythmic and DE genes #
#                                                                                                    #
#   Author: Domitille Jarrige                                                                        #
#                                                                                                    #
#   Date:   2021-06-15                                                                               #
######################################################################################################

## Librairies :

required_lib <- c("tidyverse",
                  "data.table",
                  "tibble",
                  "dplyr")

required_bioc <- c("DiscoRhythm",
                   "SummarizedExperiment",
                   "DESeq2")

### Libraries CRAN
for (lib in required_lib) {
    if (!require(lib, character.only = TRUE)) {
        install.packages(lib, quiet = TRUE)
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

message("#################################################")
message("#### DESeq2 normalisation and transformation ####")
message("#################################################")

## data loading
#meta_temp <- read.csv("../data/sample_DiscoRythm_format.tsv", sep="\t")
meta_temp <- read.csv( snakemake@input[["meta"]] , sep="\t")
#counts_rund_df <- read.csv("../results/mapping/full_data/total_counts_int.csv", sep="\t", row.names = 1)
counts_rund_df <- read.csv( snakemake@input[["counts"]] , sep="\t", row.names = 1)
n <- ncol(counts_rund_df)

## Conversion of counts data into matrix
matrix_counts_DESeq2 <- round(as.matrix(counts_rund_df[,1:n]))
rownames(matrix_counts_DESeq2) <- row.names(counts_rund_df)

# Building DESeq dataset from the matrix, design is time

ddsFull <- DESeqDataSetFromMatrix(countData = matrix_counts_DESeq2, 
                                  colData = meta_temp, 
                                  design = ~ time_factor)

## Size factors estimation
ddsFull <- estimateSizeFactors(ddsFull)

## elimination of never expressed genes
ddsFull <- ddsFull[ rowSums(counts(ddsFull)) > 0, ]

## rlog normalisation
rlog_data <- rlog(ddsFull, blind=TRUE)

## extract rlog matrix
counts_rlog <- assay(rlog_data)

message("###############################################################")
message("#### DESeq2 differential expression analysis by time point ####")
message("###############################################################")

ddsTimes <- DESeqDataSetFromMatrix(countData = matrix_counts_DESeq2, 
                                   colData = meta_temp, 
                                   design = ~ time_factor)

## Elimination of lowly expressed genes

ddsTimes <- ddsTimes[ rowSums(counts(ddsTimes)) > 60, ]

## Differential expression analysis
try( {
     ddsTimes <- DESeq(ddsTimes)
     time_points <- unique(meta_temp$time_factor)
     time_rep <- c(rep(time_points, c(1, rep(2, length(time_points)-1))), time_points[1])
     pairs <- matrix(as.factor(time_rep), ncol = 2, byrow = TRUE)
     
     # Significant genes initialisation
     significant_genes <- c()
     
     # Extraction of DESeq2 results
     for (i in (1 : nrow(pairs))) {
         ele <- pairs[i,]
         text <- paste0(ele[1], "/", ele[2], " DEG:")
         resDESeq <- results(ddsTimes, contrast = c("time_factor", ele[1], ele[2]),
                             independentFiltering = TRUE, alpha=0.01)
         message(text)
         message( sum( resDESeq$padj < 0.01, na.rm=TRUE ) )
         temp_list <- row.names(resDESeq[which(resDESeq$padj < 0.01),])
         message(paste0(c(length(temp_list)/nrow(counts_rund_df) * 100), " %"))
         significant_genes <- union(temp_list, significant_genes)
         }
     
     print(paste0("Total DEG in at least one time transition: ", (length(significant_genes)/nrow(counts_rund_df) * 100), " %"))
     
     sign_df <- data.frame(row.names = significant_genes)
     sign_df$diff_expr <- "yes"
     sign_df$gene_id <- row.names(sign_df)
     }
   )


message("##########################################")
message("#### Rhythm analysis with DiscoRhythm ####")
message("##########################################")


## DiscoRhythm

# convertion of rlog matrix to dataframe...
gene_id <- row.names(counts_rlog)
rlog_df <- as.data.frame(counts_rlog)
rlog_df <- add_column(rlog_df, gene_id, .before=1)

# convertion in SE (Submarised Experiment)
input_data <- discoDFtoSE(rlog_df)

discoDesignSummary(input_data)

# Cosinor regression method

# Détection de rythmes avec la méthode cosinor
rythms_genes_CS <- discoODAs(input_data, period = 24, method = "CS",
                             circular_t = TRUE, ncores = 4)

file_path <- paste0("./", snakemake@params[["dir"]], "qvalues_cs.png")
png(filename = file_path,
    width = 700, height = 700, res=150)
hist(data.frame(rythms_genes_CS)$CS.qvalue, breaks = 100, 
     main = "qvalue cosinor",
     xlab = "qvalue")
dev.off()


# Extracting just the dataframes from the discorythm output
CS_df <- data.frame(rythms_genes_CS)

# Retreiving only most significant genes 
# qvalue inferior to 0.05
CS_df <- CS_df[CS_df$CS.qvalue < 0.05,]

CS_df$gene_id <- row.names(CS_df)

if (exists("sign_df", mode = "any")) {
    final_data <- left_join(CS_df, sign_df, by = "gene_id")
    final_data <- final_data[which(final_data$diff_expr == "yes"),]
    final_data <- final_data[ which(!is.na(final_data$CS.acrophase)) , ]
    row.names(final_data) <- final_data$gene_id
    
}else{
    message("WARNING: differential expression analysis could not be performed on the dataset.")
    final_data <- CS_df
    if (nrow(final_data)>0) {
        final_data$diff_expr <- NA
        row.names(final_data) <- final_data$gene_id}
    else{message("WARNING: no rhytmic genes could be found in the dataset")}
    }

file_path <- paste0("./", snakemake@params[["dir"]], "rhythmic_genes_cosinor.csv")
write_csv(final_data, file = file_path)


sessionInfo()





