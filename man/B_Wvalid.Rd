\name{B_Wvalid}
\alias{B_Wvalid}
\title{
BCVI-Wiroonsri (WI) index
}
\description{
Compute Bayesian cluster validity index (BCVI) from two to \code{kmax} groups using Wiroonsri (WI) as the underling cluster validity index (CVI) with the user's selected Dirichlet prior parameters. The full detail of BCVI can be found in the paper Wiroonsri and Preedasawakul (2024).}
\usage{
B_Wvalid(x, kmax, method = "kmeans", corr = "pearson", nstart = 100,
      sampling = 1, NCstart = TRUE, alpha = "default", mult.alpha = 1/2)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"kmeans"}, \code{"hclust_complete"}, \code{"hclust_average"}, \code{"hclust_single"}). The default is \code{"kmeans"}.
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed (\code{"pearson"}, \code{"kendall"} or \code{"spearman"}). The default is \code{"pearson"}.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
  \item{sampling}{a number greater than 0 and less than or equal to 1 indicating the undersampling proportion of data to be used. This argument is intended for handling a large dataset. The default is \code{1}.}
  \item{NCstart}{logical for \code{indexlist} includes the \code{"NC"}, \code{"NCI"}, \code{"NCI1"}, and \code{"NCI2"}), if \code{TRUE}, the NC correlation at \code{k=1} is defined as the ratio introduced in the reference. Otherwise, it is assigned as \code{0}.
}
  \item{alpha}{Dirichlet prior parameters \eqn{\alpha_2,...,\alpha_k} where \eqn{\alpha_k} is the parameter corresponding to "the probability of having k groups" (selecting each \eqn{\alpha_k} between 0 to 30 is recommended and using the other parameter \code{mult.alpha} to be its multiplier. The default is \code{"default"}.
}
  \item{mult.alpha}{the power \eqn{s} from \eqn{n^s} to be multiplied to the Dirichlet prior parameters \code{alpha} (selecting \code{mult.alpha} in \code{[0,1)} is recommended). The default is \eqn{\frac{1}{2}}.
}
}
\details{
BCVI-WI is defined as follows.
Let \cr
\deqn{r_k(\bf x) = \dfrac{WI(k)-\min_j WI(j)}{\sum_{i=2}^K (WI(i)-\min_j WI(j))} }
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
\item{BCVI}{the dataframe where the first and the second columns are the number of groups \code{k} and BCVI\eqn{(k)}, respectively, for \code{k} from \code{2} to \code{kmax}.}
  \item{VAR}{the data frame where the first and the second columns are the number of groups \code{k} and the variance of \eqn{p_k}, respectively, for \code{k} from \code{2} to \code{kmax}.}
  \item{CVI}{the data frame where the first and the second columns are the number of groups \code{k} and the original WI\eqn{(k)}, respectively, for \code{k} from \code{2} to \code{kmax}.}
}
\references{
N. Wiroonsri, "Clustering performance analysis using a new correlation based cluster validity index," Pattern Recognition, 145, 109910, 2024. \doi{10.1016/j.patcog.2023.109910}

N. Wiroonsri, O. Preedasawakul, "A Bayesian cluster validity index", arXiv:2402.02162, 2024.
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{
\code{\link{B2_data}, \link{B_TANG.IDX}, \link{B_WP.IDX}, \link{B_STRPBM.IDX}, \link{B_DB.IDX}}
}
\examples{
library(BayesCVI)

# The data included in this package.
data = B2_data[,1:2]

# alpha
aalpha = c(5,5,5,20,20,20,0.5,0.5,0.5)

B.WI = B_Wvalid(x = scale(data), kmax = 10, method = "kmeans", corr = "pearson",
              nstart = 100, sampling = 1, NCstart = TRUE, alpha = aalpha,
              mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.WI)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot

}
