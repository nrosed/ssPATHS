test_that("Gene weight on reference dataset works", {
    library("dml")
    library("data.table")
    library("ROCR")
    library("MESS")
    library("ggplot2")
    library("ssPATHS")
    library("plyr")

    data(tcga_expr_df)


    hypoxia_gene_ids = get_hypoxia_genes()
    hypoxia_gene_ids = intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
    hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]

    colnames(hypoxia_df)[1:2] = c("sample_id", "Y")
    hypoxia_df$Y = 0
    hypoxia_df$Y[tcga_expr_df$is_normal==TRUE] = 0
    hypoxia_df$Y[tcga_expr_df$is_normal==FALSE] = 1

    # now we can get the gene weightings
    res = get_gene_weights(hypoxia_df)
    gene_weights_test = res[[1]]
    sample_scores = res[[2]]

    data("gene_weights_reference")
    expect_equal(gene_weights_test, gene_weights_reference)
})

test_that("test classification on reference dataset works", {
    library("dml")
    library("data.table")
    library("ROCR")
    library("MESS")
    library("ggplot2")
    library("ssPATHS")
    library("plyr")

    data(tcga_expr_df)


    hypoxia_gene_ids = get_hypoxia_genes()
    hypoxia_gene_ids = intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
    hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]

    colnames(hypoxia_df)[1:2] = c("sample_id", "Y")
    hypoxia_df$Y = 0
    hypoxia_df$Y[tcga_expr_df$is_normal==TRUE] = 0
    hypoxia_df$Y[tcga_expr_df$is_normal==FALSE] = 1

    # now we can get the gene weightings
    res = get_gene_weights(hypoxia_df)
    sample_scores = res[[2]]

    training_res = get_classification_accuracy(sample_scores, positive_val=1)
    print(training_res[[2]])

    expect_equal(training_res[[2]], 0.91112)
})


test_that("test classification on new dataset works", {
    library("dml")
    library("data.table")
    library("ROCR")
    library("MESS")
    library("ggplot2")
    library("ssPATHS")
    library("plyr")

    data(tcga_expr_df)


    hypoxia_gene_ids = get_hypoxia_genes()
    hypoxia_gene_ids = intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
    hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]

    colnames(hypoxia_df)[1:2] = c("sample_id", "Y")
    hypoxia_df$Y = 0
    hypoxia_df$Y[tcga_expr_df$is_normal==TRUE] = 0
    hypoxia_df$Y[tcga_expr_df$is_normal==FALSE] = 1

    # now we can get the gene weightings
    res = get_gene_weights(hypoxia_df)
    gene_weights = res[[1]]
    sample_scores = res[[2]]

    data(new_samp_df)
    new_score_df_calculated = get_new_samp_score(gene_weights, new_samp_df)

    data(expected_score_output)
    expect_equal(new_score_df_calculated, expected_score_output)

})

