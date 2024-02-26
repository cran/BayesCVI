\name{B_PBM.IDX}
\alias{B_PBM.IDX}
\title{
BCVI-Pakhira-Bandyopadhyay-Maulik (PBM) index
}
\description{
Compute Bayesian cluster validity index (BCVI) from two to \code{kmax} groups using Pakhira-Bandyopadhyay-Maulik (PBM) as the underling cluster validity index (CVI) with the user's selected Dirichlet prior parameters. The full detail of BCVI can be found in the paper Wiroonsri and Preedasawakul (2024).}
\usage{
B_PBM.IDX(x, kmax, method = "FCM", fzm = 2, nstart = 20,
        iter = 100, alpha = "default", mult.alpha = 1/2)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"FCM"} or \code{"EM"}). The default is \code{"FCM"}.
}
  \item{fzm}{a number greater than 1 giving the degree of fuzzification for \code{method = "FCM"}. The default is \code{2}.
}
  \item{nstart}{a maximum number of initial random sets for FCM for \code{method = "FCM"}. The default is \code{20}.
}
  \item{iter}{a maximum number of iterations for \code{method = "FCM"}. The default is \code{100}.
}
  \item{alpha}{Dirichlet prior parameters \eqn{\alpha_2,...,\alpha_k} where \eqn{\alpha_k} is the parameter corresponding to "the probability of having k groups" (selecting each \eqn{\alpha_k} between 0 to 30 is recommended and using the other parameter \code{mult.alpha} to be its multiplier. The default is \code{"default"}.
}
  \item{mult.alpha}{the power \eqn{s} from \eqn{n^s} to be multiplied to the Dirichlet prior parameters \code{alpha} (selecting \code{mult.alpha} in \code{[0,1)} is recommended). The default is \eqn{\frac{1}{2}}.
}
}
\details{
BCVI-PBM is defined as follows.
Let \cr
\deqn{r_k(\bf x) = \dfrac{PBM(k)-\min_j PBM(j)}{\sum_{i=2}^K (PBM(i)-\min_j PBM(j))} }
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
  \item{CVI}{the data frame where the first and the second columns are the number of groups \code{k} and the original PBM\eqn{(k)}, respectively, for \code{k} from \code{2} to \code{kmax}.}
}
\references{
M. K. Pakhira, S. Bandyopadhyay, and U. Maulik, “Validity index for crisp and fuzzy clusters,” Pattern recognition, vol. 37, no. 3, pp. 487–501, 2004.

N. Wiroonsri, O. Preedasawakul, "A Bayesian cluster validity index", arXiv:2402.02162, 2024.
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{
\code{\link{B7_data}, \link{B_TANG.IDX}, \link{B_WP.IDX}, \link{B_Wvalid}, \link{B_DB.IDX}}
}
\examples{
library(BayesCVI)

# The data included in this package.
data = B7_data[,1:2]

# alpha
aalpha = c(5,5,5,20,20,20,0.5,0.5,0.5)

B.PBM = B_PBM.IDX(x = scale(data), kmax =10, method = "FCM", fzm = 2, nstart = 20,
                iter = 100, alpha = aalpha, mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.PBM)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot
}
