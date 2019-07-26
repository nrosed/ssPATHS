\name{get_gene_weights}
\alias{get_gene_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get Gene Weights from Reference Data
}
\description{
  This method performs linear discriminant analysis on a reference dataset
  using a pre-defined set of genes related to a pathway of interest.
}
\usage{
get_gene_weights(expression_matr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_matr}{
   This is an expression matrix of the reference samples. Rows are
   samples and columns are Y (whether this is a positive or negative sample),
   sample_id (unique sample id), gene expression. All columns not titled "Y"
   or "sample_id" are assumed to be gene names in the pathway of interest.
}
}

\value{
  A list containing the gene weights and estimated scores of the reference
  samples.
  \item{proj_vector_df }{A dataframe containing the gene weights and gene ids}
  \item{dca_proj }{A dataframe containing the sample scores and sample ids.}
%% ...
}
\references{
 Steven C.H. Hoi, W. Liu, M.R. Lyu and W.Y. Ma (2006). Learning Distance
 Metrics with Contextual Constraints for Image Retrieval. Proceedings IEEE
 Conference on Computer Vision and Pattern Recognition (CVPR2006).
}
\author{
  Natalie R. Davidson
}
\note{
%%  ~~further notes~~
}


\examples{

library("dml")
library("data.table")
library("ROCR")
library("MESS")

# read in the reference expression data for hypoxia score generation
print(getwd())
data(tcga_expr_df)

# take a look at the format, we need a Y column and a sample_id column.
tcga_expr_df[1:10,1:10]

# let's get the info and columns we need to learn a hypoxia score
hypoxia_gene_ids = get_hypoxia_genes()
hypoxia_gene_ids = intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]
colnames(hypoxia_df)[1:2] = c("sample_id", "Y")

# format the Y column such that normal samples are 0 and tumor samples are 1
hypoxia_df$Y = 0
hypoxia_df$Y[tcga_expr_df$is_normal==TRUE] = 0
hypoxia_df$Y[tcga_expr_df$is_normal==FALSE] = 1

# now we can get the gene weightings
res = get_gene_weights(hypoxia_df)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
