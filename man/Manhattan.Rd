\name{Manhattan}
\alias{plot.NAM}
\title{
Manhattan plot for Association Studies 
}
\description{
Generates a graphical visualization for the output of the function gwas/gwas2. 
}
\usage{
\method{plot}{NAM}(x,...,alpha=0.05,colA=2,colB=4,find=NULL,FDR=NULL,gtz=FALSE)
}
\arguments{
   \item{x}{
Output of the gwas/gwas2 function.
}
  \item{...}{
Further arguments passed to or from other methods.
}
  \item{alpha}{
Numberic. Significance threshold to display in the Manhattan plot.
}
  \item{colA}{
Color of odd chromosomes in the Manhattan plot.
}
  \item{colB}{
Color of even chromosomes in the Manhattan plot.
}
  \item{find}{
Integer. If provided, you can click on the specified number of hits in the Manhattan plot to obtain the name of the markers.
}
  \item{FDR}{
Null or numeric between zero and one. If provided, it will display the Manhattan plot with Bonferroni threshold by chromosome, adjusted for the specified false discovery rate (FDR). Thus, zero provides the Bonferroni correction.
}
  \item{gtz}{
Logical. If TRUE, the argument FDR will just take into account markers with p-value Greater Than Zero (GTZ).
}
}
\author{
Alencar Xavier and William Beavis
}
\examples{
data(tpod)
test=gwas2(y=y,gen=gen,fam=fam,chr=chr)

par(mfrow=c(2,1))

# Example Manhattan 1
plot(x=test,colA=3,colB=1,alpha=0.05/ncol(gen),type="h",main="Example 1: Genome-Wide Bonferroni")

# Example Manhattan 2
Title = "Example 2: FDR 0.25 by Chromosome"
plot(x=test,alpha=0.05,FDR=0.25,gtz=TRUE,pch=20,main=Title)
}
\keyword{ NAM }
\keyword{ GWAS }
