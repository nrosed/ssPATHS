
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# pathway activation score functions


get_hypoxia_genes <- function(){

    HIF_GENES <- c("ENSG00000148926", "ENSG00000109107", "ENSG00000176171",
                  "ENSG00000104765", "ENSG00000074410", "ENSG00000107159",
                  "ENSG00000130635", "ENSG00000047457", "ENSG00000168209",
                  "ENSG00000129521", "ENSG00000111674", "ENSG00000104812",
                  "ENSG00000159399", "ENSG00000100292", "ENSG00000134333",
                  "ENSG00000113083", "ENSG00000123384", "ENSG00000185499",
                  "ENSG00000119950", "ENSG00000104419", "ENSG00000185633",
                  "ENSG00000124785", "ENSG00000152256", "ENSG00000114268",
                  "ENSG00000204531", "ENSG00000119938", "ENSG00000139832",
                  "ENSG00000141526", "ENSG00000117394", "ENSG00000103257",
                  "ENSG00000113739", "ENSG00000265972", "ENSG00000112715",
                  "ENSG00000186918", "ENSG00000117289")
    return(HIF_GENES)
}

get_gene_weights <- function(expression_se){

    sample_info_names <- colnames(colData(expression_se))
    if(sum(sample_info_names %in% c("Y", "sample_id")) != 2){
        stop("Need column names Y and sample_id")
    }
    if(length(unique(colData(expression_se)$Y)) != 2){
        stop("Y must be binary")
    }
    if(!(1 %in% colData(expression_se)$Y) | !(0 %in% colData(expression_se)$Y)){
        stop("Y must have both a 1 and 0 entry")
    }

    gene_ids <- rownames(expression_se)

    dca_matr <- t(assay(expression_se))
    chunks <- rep(1, ncol(expression_se))
    chunks[colData(expression_se)$Y!=0] <- 2
    neglinks <- matrix(c(0, 1, 1, 0), 2, 2)

    # normalize
    dca_data <- as.matrix(t(scale(log10(t(dca_matr+1)))))
    dca_data[is.nan(dca_data)] <- 0
    row.names(dca_data) <- colData(expression_se)$sample_id

    # get weights
    dca_res <- dml::dca(data=dca_data, chunks=chunks, neglinks=neglinks)
    proj_vector <- t(as.matrix(dca_res$DCA))

    # get projection
    dca_proj <- dca_data %*% proj_vector
    dca_proj <- data.frame(pathway_score=dca_proj, sample_id=row.names(dca_proj))
    dca_proj <- merge(colData(expression_se), dca_proj)
    dca_proj <- dca_proj[order(dca_proj$pathway_score),]

    # check direction
    lower_score <- mean(dca_proj$pathway_score[dca_proj$Y==0])
    upper_score <- mean(dca_proj$pathway_score[dca_proj$Y==1])
    if(lower_score > upper_score){

        proj_vector <- proj_vector * -1
        dca_proj <- dca_data %*% proj_vector
        dca_proj <- data.frame(pathway_score=dca_proj, sample_id=row.names(dca_proj))
        dca_proj <- merge(colData(expression_se), dca_proj)
        dca_proj <- dca_proj[order(dca_proj$pathway_score),]

    }

    proj_vector_df <- data.frame(gene_weight=proj_vector, gene_id=gene_ids)

    return(list(proj_vector_df, dca_proj))

}


get_classification_accuracy <- function(sample_scores, positive_val){

    if(sum(colnames(sample_scores) %in% c("sample_id", "pathway_score")) != 2){
        stop("sample_scores need column names pathway_score and sample_id")
    }

    pred_dca <- ROCR::prediction(sample_scores$pathway_score, sample_scores$Y == positive_val)
    perf_dca_roc <- ROCR::performance(pred_dca, "tpr", "fpr")
    auc_dca <- ROCR::performance(pred_dca, "auc")
    auc_roc <- unlist(auc_dca@y.values)

    perf_dca_pr <- ROCR::performance(pred_dca, "prec", "rec")

    x <- perf_dca_pr@x.values[[1]] # Recall values
    y <- perf_dca_pr@y.values[[1]]

    auc_pr <- try(MESS::auc(x,y, type = 'spline'), TRUE)
    if(inherits(auc_pr, "try-error")){
        auc_pr <- NA
        warning("AUC for PR curve could not be calculated")
    }


    return(list(auc_pr=auc_pr, auc_roc=auc_roc, perf_pr=perf_dca_pr, perf_roc=perf_dca_roc))

}

get_new_samp_score <- function(gene_weights, expression_se, run_normalization=TRUE){

    if(sum(colnames(colData(expression_se)) %in% c("sample_id")) != 1){
        stop("Need column name sample_id")
    }

    if(sum(colnames(gene_weights) %in% c("gene_weight", "gene_id")) != 2){
        stop("Need column names gene_weight and gene_id")
    }
    has_Y <- "Y" %in% colnames(colData(expression_se))

    gene_ids <- rownames(expression_se)

    if(length(gene_ids) != length(gene_weights$gene_id)){
        warning("Genes missing in gene_weights or expression_se")
    }
    # normalize
    dca_matr <- t(assay(expression_se))
    if(run_normalization){
        dca_data <- as.matrix(t(scale(log10(t(dca_matr+1)))))
    }else{
        dca_data <- as.matrix(t(scale((t(dca_matr+1)))))
    }
    dca_data[is.nan(dca_data)] <- 0
    row.names(dca_data) <- colData(expression_se)$sample_id

    # format the projection vector
    proj_vector <- gene_weights$gene_weight
    names(proj_vector) <- gene_weights$gene_id
    proj_vector <- proj_vector[gene_ids]

    if(sum(names(proj_vector) != colnames(dca_matr)) > 0){
        stop("Error Matching Gene ids")
    }

    # get projection
    dca_proj <- dca_data %*% proj_vector
    dca_proj <- data.frame(pathway_score=dca_proj, sample_id=row.names(dca_proj))
    dca_proj <- merge(colData(expression_se), dca_proj, by="sample_id")
    dca_proj <- dca_proj[order(dca_proj$pathway_score),]

    if(has_Y){
        dca_proj <- dca_proj[,c("sample_id", "Y", "pathway_score")]
    }else{
        dca_proj <- dca_proj[,c("sample_id", "pathway_score")]
    }

    return(dca_proj)

}
