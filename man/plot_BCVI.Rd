\name{plot_BCVI}
\alias{plot_BCVI}
\title{
Plots for visualizing BCVI
}
\description{
Plot Bayesian cluster validity index (BCVI) with and without standard deviation error bars and the underlying index.
}
\usage{
plot_BCVI(B.result, mult.err.bar = 2)
}
\arguments{
  \item{B.result}{a result from one of the functions \code{B_XB.IDX, B_Wvalid, B_WP.IDX, B_WL.IDX, B_TANG.IDX, B_STRPBM.IDX, B_SF.IDX, B_PBM.IDX, B_PB.IDX, B_KWON.IDX, B_KWON2.IDX, B_KPBM.IDX, B_HF.IDX, B_GC.IDX, B_DI.IDX, B_DB.IDX, B_CSL.IDX, B_CH.IDX, B_CCV.IDX and B_BayesCVIs.IDX}
}
  \item{mult.err.bar}{
a multiplier of the stadard deviations to be used for plotting error bars
}
}
\details{BCVI is defined as follows.

Let \cr
\deqn{r_k(\bf x) = \dfrac{\max_j CVI(j)- CVI(k)}{\sum_{i=2}^K (\max_j CVI(j) - CVI(i))}} for a cluster validity index (CVI) such that the smallest value indicates the optimal number of clusters and
\deqn{r_k(\bf x) = \dfrac{CVI(k)-\min_j CVI(j)}{\sum_{i=2}^K (CVI(i)-\min_j CVI(j))} }
for a CVI such that the largest indicates the optimal number of clusters.
Assume that \cr
\deqn{f({\bf x}|{\bf p}) = C({\bf p}) \prod_{k=2}^Kp_k^{nr_k(x)}}
represents the conditional probability density function of the dataset given \eqn{\bf p}, where \eqn{C({\bf p})} is the normalizing constant. Assume further that  \eqn{{\bf p}} follows a Dirichlet prior distribution with parameters \eqn{{\bm \alpha} = (\alpha_2,\ldots,\alpha_K)}. The posterior distribution of \eqn{\bf p} still remains a Dirichlet distribution with parameters \eqn{(\alpha_2+nr_2({\bf x}),\ldots,\alpha_K+nr_K({\bf x}))}.

The BCVI is then defined as  \cr
\deqn{BCVI(k) = E[p_k|{\bf x}] = \frac{\alpha_k + nr_k({\bf x})}{\alpha_0+n}}
where \eqn{\alpha_0 = \sum_{k=2}^K \alpha_k.} \cr

The variance of \eqn{p_k} can be computed as
\deqn{Var(p_k|{\bf x}) = \dfrac{(\alpha_k + nr_k(x))(\alpha_0 + n -\alpha_k-nr_k(x))}{(\alpha_0 + n)^2(\alpha_0 + n +1 )}.}
}
\value{
  \item{plot_index}{a plot of the underlying index for the number of groups from \eqn{2} to \eqn{kmax} according to \code{B.result}}
  \item{plot_BCVI}{a plot of BCVI for the number of groups from \eqn{2} to \eqn{kmax} according to \code{B.result}}
  \item{error_bar_plot}{a plot of BCVI with error bars for the number of groups from \eqn{2} to \eqn{kmax} according to \code{B.result}}
}
\references{
N. Wiroonsri, O. Preedasawakul, "A Bayesian cluster validity index", arXiv:2402.02162, 2024.}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{\code{\link{B_STRPBM.IDX}, \link{B_TANG.IDX}, \link{B_XB.IDX}, \link{B_Wvalid}, \link{B_WP.IDX}, \link{B_DB.IDX}}
}
\examples{

library(BayesCVI)
library(UniversalCVI)

##Soft clustering

# The data included in this package.
data = B7_data[,1:2]

# alpha
aalpha = c(5,5,5,20,20,20,0.5,0.5,0.5)

B.XB = B_XB.IDX(x = scale(data), kmax =10, method = "FCM", fzm = 2,
              nstart = 20, iter = 100, alpha = aalpha, mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.XB)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot


## Hard clustering

# The data included in this package.
data = B2_data[,1:2]

K.STR = STRPBM.IDX(scale(data), kmax = 10, kmin = 2, method = "kmeans",
  indexlist = "STR", nstart = 100)

# WP.IDX values
result = K.STR$STR$STR


aalpha = c(20,20,20,5,5,5,0.5,0.5,0.5)
B.STR = BayesCVIs(CVI = result,
          n = nrow(data),
          kmax = 10,
          opt.pt = "max",
          alpha = aalpha,
          mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.STR)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot

}
