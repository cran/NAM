\name{SNP QC}
\alias{snpQC}
\alias{cleanREP}
\title{
SNP Quality Control
}
\description{
A function for quality control. It may be used to count/remove neighbor repeated SNPs and markers with MAF lower than a given threshold. This function is also used for imputations.
}
\usage{
snpQC(gen,psy=1,MAF=0.05,remove=TRUE,impute=FALSE)
}
\arguments{
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozigous, heterozigous and reference homozigous. \code{NA} is allowed.
}
  \item{psy}{
Tolerance parameter for repeated markers. Default is 1, which removes only SNPs 100\% equal to its following neighbor.
}
  \item{MAF}{
Minor Allele Frequency. Default is 0.05. Useful to inform or remove markers below the MAF threshold. 
}
  \item{remove}{
Remove SNPs that are redundant or pursue low MAF: TRUE/FALSE.
}
  \item{impute}{
If TRUE, impute missing values using Random Forest implemented in the package missForest. Methods are described in Rutkoski et al (2013).
}
}
\value{
Returns the genomic matrix without missing, redundant or low MAF markers.
}
\references{
Rutkoski, J. E., Poland, J., Jannink, J. L., & Sorrells, M. E. (2013). Imputation of unordered markers and the impact on genomic selection accuracy. G3: Genes| Genomes| Genetics, 3(3), 427-439.
}
\author{
Alencar Xavier, Katy Rainey and William Muir
}
\examples{
data(tpod)
gen=snpQC(gen=gen,psy=1,MAF=0.05,remove=TRUE,impute=FALSE)
}