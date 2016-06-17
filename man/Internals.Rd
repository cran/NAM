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
\alias{.Random.seed}
\title{
  Internal functions
}
\description{
  Complimentary statistics and functions written in C++ to speed up \eqn{gwas}, \eqn{gibbs} and \eqn{wgr}.
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

# Spatial correlation kernel
covar()

# Genetic distance
round(Gdist(gen[1:10,],method=1),2)

# PCs of a NAM kinship
eG = eigX(gen,fam)
plot(eG[[2]],col=fam)
}