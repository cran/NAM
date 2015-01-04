\name{manhattan}
\alias{manhattan}
\title{
Plotting Genome-Wide Association Study
}
\description{
Generates a Manhattan plot for a \eqn{gwas} object.
}
\usage{
manhattan(gwas, colA = 2, colB = 4, GenDist = FALSE, OtherDist = NULL, ...)
}
\arguments{
  \item{gwas}{
An output of the function gwas.
}
  \item{colA}{
Color of even chromosomes.
}
  \item{colB}{
Color of odd chromosomes.
}
  \item{GenDist}{
Logical. If true, it displays the Manhattan plot with the estimated genetic map for NAM populations.
}
  \item{OtherDist}{
Numeric vector of length (\eqn{m}), equal to the number of columns in \eqn{gen}. If provided, it changes the X coordinates used for the Manhattan plot.
}
  \item{...}{
Other arguments for plotting function.
}
}
\author{
Alencar Xavier
}
\examples{
data(tpod)
test=gwas(y=y,gen=gen,fam=fam,chr=chr)
manhattan(gwas=test,colA=3,colB=1,type="h",lwd=2)
manhattan(gwas=test,GenDist=TRUE,pch=20)
}