\name{Manhattan Plot}
\alias{manhattan}
\title{
Plotting Genome-Wide Association Study
}
\description{
Generates a Manhattan plot for a \eqn{gwas} object.
}
\usage{
manhattan(gwas, colA = 2, colB = 4, alpha=0.05, ...)
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
  \item{alpha}{
If specified, it provides the negative log p-value based on a chi squared distribution with a significance threshold for the given alpha.
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
test=gwas(y=y,gen=gen[,1:240],fam=fam,chr=chr[1:12])
manhattan(gwas=test,colA=3,colB=1,type="h",lwd=2)
ii = 1+(2*test$PolyTest$lrt>4.9)
manhattan(gwas=test,pch=ii+3,lwd=ii)
}