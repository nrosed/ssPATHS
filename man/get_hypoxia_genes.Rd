\name{get_hypoxia_genes}
\alias{get_hypoxia_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Ensembl ids of hypoxia related genes.
}
\description{
Returns a vector of Ensembl ids of hypoxia related genes.
}
\usage{
get_hypoxia_genes()
}

\value{
Vector of ensembl ids.
}

\author{
Natalie R. Davidson
}


\examples{


# read in the reference expression data for hypoxia score generation
data(tcga_expr_df)

# let's get the expression of hypoxia associated genes
hypoxia_gene_ids <- get_hypoxia_genes()
hypoxia_gene_ids <- intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
hypoxia_df <- tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]


}
