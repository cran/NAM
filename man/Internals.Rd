\name{Internals}
\alias{calcSize}
\alias{funI}
\alias{funX}
\alias{inputRow}
\alias{RcppExports}
\alias{timesMatrix}
\alias{timesVec}
\alias{SAMP}
\alias{SAMP2}
\alias{gs}
\alias{KMUP}
\alias{Gdist}
\alias{markov}
\alias{eigX}
\alias{covar}
\alias{PedMat}
\alias{PedMat2}
\alias{LD}
\alias{G2A_Kernels}
\alias{Import_data}
\alias{.Random.seed}
\alias{NNsrc}
\alias{NNcov}
\title{
  Internal functions
}
\description{
Function under optimization, complimentary statistics and loops written in C++ to speed up \eqn{gwas}, \eqn{gibbs} and \eqn{wgr}.
}
\details{

Some of the functions available for users include:

01) \code{Import_data(file,type=c('GBS','HapMap','VCF')):}{  This function can be used to import genotypic data in the NAM format, providing a list with a genotypic matrix \code{gen} coded as \code{012} and a vector \code{chr} with count of markers per chromosome. Currently, it helps users to import three types of files: GBS text, HapMap and VCF.}

02) \code{markov(gen,chr):}{  Imputation method based forwards Markov model for SNP data coded as \code{012}. We recommend users to remove non-segregating markers before using this function.}

03) \code{LD(gen):}{  Computes the linkage disequilibrium in terms of r2 for SNP data coded as \code{012}. Missing data is not allowed.}

04) \code{PedMat(ped):}{  Builds a kinship from a pedigree. Input format is provided with \code{PedMat()}.}

05) \code{PedMat2(ped,gen=NULL,IgnoreInbr=FALSE,PureLines=FALSE):}{  Builds a kinship from a genomic data and pedigree. Useful when not all individuals are genotyped. Row names of \code{gen} must indicate the genotype id.}

06) \code{Gdist(gen, method = 1):}{  Computes genetic distance among individuals. Five methods are available: 1) Nei distance; 2) Edwards distance; 3) Reynolds distance; 4) Rogers distance; 5) Provesti's distance. 6) Modified Rogers distance}

07) \code{covar(sp=NULL,rho=3.5,type=1,dist=2.5):}{  Builds a spatial kernel from field plot information. Input format is provided with \code{covar()}. Parameter \code{rho} detemines the decay of relationship among neighbor plots. \code{type} defines if the kernel is exponential (1), Gaussian (2) or some intermediate. \code{dist} informs the distance ratio between range neighbors and row neighbors.}

08) \code{Gdist(gen,fam):}{  Computes the input of the argument \code{EIG} of the function \code{gwas2}.}

09) \code{G2A_Kernels(gen):}{  Computes a list of orthogonal kernels containing additive, dominant and first-order epistatic effects, in accordance to the G2A model from ZB Zeng et al. 2005. These kernels can be used for description of genetic architecture through variance components, for that we recommend packages \code{varComp} and \code{BGLR}.}

10) \code{NNsrc(sp=NULL,rho=1,dist=3):}{  Using the same field data input required by the function \code{covar}, this function provides a list of nearest neighbor plots for each entry.}

11) \code{NNcov(NN,y):}{  This function utilizes the output of \code{NNsrc} to generate a numeric vector, averageing the observed values of \code{y}. This function is useful to generate field covariates to control micro-environmental variance without krigging.}

}
\author{
  Alencar Xavier and Tiago Pimenta
}
\examples{
# Forward gen imputation
data(tpod)
fast.impute = markov(gen,chr)

# A matrix
PedMat()

# Pairwise LD
ld = LD(gen[,1:10])
heatmap(ld)

# Spatial correlation (kernel-based)
covar()

# Spatial correlation (NN-based)
NNsrc()

# Genetic distance
round(Gdist(gen[1:10,],method=1),2)

# PCs of a NAM kinship
eG = eigX(gen,fam)
plot(eG[[2]],col=fam)

# Polygenic kinship matrices
Ks = G2A_Kernels(gen)
ls(Ks)
}