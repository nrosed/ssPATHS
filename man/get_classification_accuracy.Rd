\name{get_classification_accuracy}
\alias{get_classification_accuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Classification Accuracy
}
\description{
Get the AUC-ROC, AUC-PR, and ROC/PR curves for plotting.
}
\usage{
get_classification_accuracy(sample_scores, positive_val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_scores}{
   This is a data.frame containing the sample id, score, and true label Y. This
   object is returned by the method get_gene_weights.
}
  \item{positive_val}{
   This is the value that will denote a true positive. It must be one of the two
   values in the Y column in sample_scores.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This returns a list of performance metrics
  \item{auc_pr }{Area under the PR-curve}
  \item{auc_roc }{Area under the ROC-curve}
  \item{perf_pr }{ROCR object for plotting the PR-curve}
  \item{perf_roc }{ROCR object for plotting the ROC-curve}

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
gene_weights = res[[1]]
sample_scores = res[[2]]

# now lets see how well we did in seperating the two classes defined by Y
training_res = get_classification_accuracy(sample_scores, positive_val=1)

plot(training_res[[3]], col="orange", ylim=c(0, 1))
legend(0.1,0.8,c(training_res$auc_pr,"\n"),border="white",cex=1.7,box.col = "white")

plot(training_res[[4]], col="blue", ylim=c(0, 1))
legend(0.1,0.8,c(training_res$auc_roc,"\n"),border="white",cex=1.7,box.col = "white")

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
