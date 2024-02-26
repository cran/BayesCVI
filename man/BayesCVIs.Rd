\name{BayesCVIs}
\alias{BayesCVIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian cluster validity index
}
\description{
Compute Bayesian cluster validity index (BCVI) from two to \code{kmax} groups using an underlying cluster validity index (CVI) and Dirichlet prior parameters of the user's choice. The full detail of BCVI can be found in the paper Wiroonsri and Preedasawakul (2024).
}
\usage{
BayesCVIs(CVI, n, kmax, opt.pt, alpha = "default", mult.alpha = 1/2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CVI}{the CVI values for \code{k} from \code{2} to \code{kmax} to be used as the underlying index for computing BCVI.
}
  \item{n}{a number of data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{opt.pt}{a character string indicating whether the maximum or the minimum of \code{CVI} specifies the optimal number of groups (\code{"min"} or \code{"max"}).
}
  \item{alpha}{Dirichlet prior parameters \eqn{\alpha_2,...,\alpha_k} where \eqn{\alpha_k} is the parameter corresponding to "the probability of having k groups" (selecting each \eqn{\alpha_k} between 0 to 30 is recommended and using the other parameter \code{mult.alpha} to be its multiplier. The default is \code{"default"}.)
}
  \item{mult.alpha}{the power \eqn{s} from \eqn{n^s} to be multiplied to the Dirichlet prior parameters \code{alpha} (selecting \code{mult.alpha} in \code{[0,1)} is recommended). The default is \eqn{\frac{1}{2}}.
}
}
\details{
BCVI is defined as follows.
Let \cr
\deqn{r_k(\bf x) = \dfrac{\max_j CVI(j)- CVI(k)}{\sum_{i=2}^K (\max_j CVI(j) - CVI(i))}} for a CVI such that the smallest value indicates the optimal number of clusters and
\deqn{r_k(\bf x) = \dfrac{CVI(k)-\min_j CVI(j)}{\sum_{i=2}^K (CVI(i)-\min_j CVI(j))} }
for a CVI such that the largest value indicates the optimal number of clusters.
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
  \item{CVI}{the data frame where the first and the second columns are the number of groups \code{k} and the original CVI\eqn{(k)}, respectively, for \code{k} from \code{2} to \code{kmax}.}
  \item{opt.pt}{a character string indicating whether the maximum or the minimum of \code{CVI} specifies the optimal number of groups (\code{"min"} or \code{"max"}) that user select.
  }
}
\references{
N. Wiroonsri, O. Preedasawakul, "A Bayesian cluster validity index", arXiv:2402.02162, 2024.
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}


\seealso{
\code{\link{B2_data}, \link{B_TANG.IDX}, \link{B_WP.IDX}, \link{B_Wvalid}, \link{B_DB.IDX}}
}
\examples{

# install a package for computing an underlying CVI
# install.packages("UniversalCVI")

library(UniversalCVI)
library(BayesCVI)

data = R1_data[,-3]

# Compute WP index by WP.IDX using default gamma
FCM.WP = WP.IDX(scale(data), cmax = 10, cmin = 2, corr = 'pearson', method = 'FCM', fzm = 2,
                iter = 100, nstart = 20, NCstart = TRUE)


# WP.IDX values
result = FCM.WP$WP$WPI


aalpha = c(20,20,20,5,5,5,0.5,0.5,0.5)
B.WP = BayesCVIs(CVI = result,
          n = nrow(data),
          kmax = 10,
          opt.pt = "max",
          alpha = aalpha,
          mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.WP)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot

}
