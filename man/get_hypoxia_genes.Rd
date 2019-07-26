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
%- maybe also 'usage' for other objects documented here.

\value{
Vector of ensembl ids.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Natalie R. Davidson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


library("dml")
library("data.table")
library("ROCR")
library("MESS")

# read in the reference expression data for hypoxia score generation
data(tcga_expr_df)

# let's get the expression of hypoxia associated genes
hypoxia_gene_ids = get_hypoxia_genes()
hypoxia_gene_ids = intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
