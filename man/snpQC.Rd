\name{SNP QC}
\alias{snpQC}
\title{
SNP Quality Control
}
\description{
A function for quality control. It may be used to count/remove neighbor repeated SNPs and markers with MAF lower than a given threshold. This function is also used for imputations.
}
\usage{
snpQC(gen,psy=1,MAF=0.05,misThr=0.8,remove=TRUE,impute=FALSE)
}
\arguments{
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozigous, heterozigous and reference homozigous. \code{NA} is allowed.
}
  \item{psy}{
Tolerance parameter for markers in Perfect SYymmetry (psy). This QC remove identical markers (aka. full LD) that carry the same information. Default is 1, which removes only SNPs 100\% equal to its following neighbor.
}
  \item{MAF}{
Minor Allele Frequency. Default is 0.05. Useful to inform or remove markers below the MAF threshold. Markers with standard deviation below the MAF threshold will be also removed.
}
  \item{misThr}{
Missing value threshold. Default is 0.8, removing markers with more than 80 percent missing values. 
}
  \item{remove}{
Logical. Remove SNPs due to PSY or MAF.
}
  \item{impute}{
If TRUE, impute missing values using Random Forest adapted from the package missForest (Stekhoven and Buhlmann 2012) as suggested by Rutkoski et al (2013).
}
}
\value{
Returns the genomic matrix without missing values, redundancy or low MAF markers.
}
\references{
Rutkoski, J. E., Poland, J., Jannink, J. L., & Sorrells, M. E. (2013). Imputation of unordered markers and the impact on genomic selection accuracy. G3: Genes| Genomes| Genetics, 3(3), 427-439.

Stekhoven, D. J. and Buhlmann, P. 2012. MissForest - nonparametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
}
\author{
Alencar Xavier, Katy Rainey and William Muir
}
\examples{
data(tpod)
gen=reference(gen)
gen=snpQC(gen=gen,psy=1,MAF=0.05,remove=TRUE,impute=FALSE)
}