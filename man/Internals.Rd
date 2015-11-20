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
\alias{KMUP}
\alias{Gdist}
\alias{markov}
\alias{covar}
\alias{PedMat}
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
fast.impute.gen = markov(gen,chr)

# A matrix
PedMat()

# Pairwise LD
ld = LD(gen[,1:3])
ld = LD(gen[,1:10])
heatmap(ld$r2)

# Spatial correlation kernel
covar()

# Genetic distance
round(Gdist(gen[1:10,],method=1),2)
}