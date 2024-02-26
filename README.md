
# BayesCVI

**BayesCVI** package is developed for computing and generating plots with and without error bars for Bayesian cluster validity index (BCVI), introduced in Wiroonsri and Preedasawakul(2024), based on several underlying cluster validity indices (CVIs) as listed below. It also allows users to input any other CVIs of their choices. The package is compatible with K-means, fuzzy C means, EM clustering, and hierarchical clustering (single, average, and complete linkage).
  
  
BayesCVI requires the use of the four packages:  [`e1071`](https://CRAN.R-project.org/package=e1071), [`mclust`](https://CRAN.R-project.org/package=mclust) for performing the fuzzy C-means (FCM) and EM algorithms, respectively, [`ggplot2`](https://CRAN.R-project.org/package=ggplot2) for plotting BCVI and existing CVIs, and [`UniversalCVI`](https://CRAN.R-project.org/package=UniversalCVI) for some required datasets.

In addition to the evaluation tools, the BayesCVI package also includes 7 simulated datasets intially used for testing BCVI in several perspectives written in Wiroonsri and Preedasawakul(2024). 

The underlying CVIs available in this package are listed as follows:

**Hard clustering:**

Dunn's index, Calinski–Harabasz index, Davies–Bouldin’s index, Point biserial correlation index, Chou-Su-Lai measure, Davies–Bouldin*’s index, Score function, Starczewski index, Pakhira–Bandyopadhyay–Maulik (for crisp clustering) index, and Wiroonsri index.

**Fuzzy clustering:**

Xie–Beni index, KWON index, KWON2 index, TANG index , HF index, Wu–Li index, Pakhira–Bandyopadhyay–Maulik (for fuzzy clustering) index, KPBM index, Correlation Cluster Validity index, Generalized C index, Wiroonsri and Preedasawakul index.

**Remark**

Though BCVI is compatible with any underlying existing CVIs, we recommend users to use either WI or WP as the underlying CVI. BCVI is only effective when underlying indices are present, providing meaningful options for ranking local peaks for the final number of clusters. This point has only been tested with either WI or WP indices.

## Installation

If you have not already installed `mclust`, `e1071`, `ggplot2` and `UniversalCVI` in your local system, install these package as follows: 

``` r
install.packages(c('e1071','mclust','ggplot2','UniversalCVI'))
```
Install `BayesCVI` package

``` r
install.packages('BayesCVI')
```

``` r
 suppressPackageStartupMessages({
library(BayesCVI) 
library(UniversalCVI)
library(e1071)
library(mclust)
library(ggplot2)
})
```

## Example

### Compute BCVI for hard clustering
Use B_Wvalid to compute BCVI with WI as the underlying CVI for a clustering results from 2 to 10 groups:

``` r
library(BayesCVI)

# The data included in this package.
data = B2_data[,1:2]

# alpha
aalpha = c(5,5,5,20,20,20,0.5,0.5,0.5)

B.WI = B_Wvalid(x = scale(data), kmax = 10, method = "kmeans", corr = "pearson", nstart = 100, sampling = 1, NCstart = TRUE, alpha = aalpha, mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.WI)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot

```

### Compute BCVI for soft clustering
Use B_WP.IDX to compute BCVI with WP as the underlying CVI for a clustering results from 2 to 10 groups:

``` r
library(BayesCVI)

# The data included in this package.
data = B7_data[,1:2]

# alpha
aalpha = c(20,20,20,5,5,5,0.5,0.5,0.5)

B.WP = B_WP.IDX(x = scale(data), kmax =10, corr = "pearson", method = "FCM",
                fzm = 2, sampling = 1, iter = 100, nstart = 20, NCstart = TRUE,
                alpha = aalpha, mult.alpha = 1/2)

# plot the BCVI

pplot = plot_BCVI(B.WP)
pplot$plot_index
pplot$plot_BCVI
pplot$error_bar_plot

```
### Compute BCVI 
Use BayesCVIS to compute BCVI with any selected underlying CVI for a clustering results from 2 to 10 groups:

``` r

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

```

### MRI brain tumor dataset
Use B_WP.IDX to compute BCVI with WP as the underlying CVI for a clustering results from 2 to 8 groups:


``` r

library(UniversalCVI)
library(BayesCVI)
library(imager)

# Download MRI data from https://www.kaggle.com/datasets/navoneel/brain-mri-images-for-brain-tumor-detection

x = "https://storage.googleapis.com/kagglesdsdata/datasets/165566/377107/yes/Y164.JPG?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=databundle-worker-v2%40kaggle-161607.iam.gserviceaccount.com%2F20240218%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20240218T124934Z&X-Goog-Expires=345600&X-Goog-SignedHeaders=host&X-Goog-Signature=269c3888a6cdc0cb4e9ea127d1e7bef2ecd798260c164acaec727c9cfa19a77428ac3ef792f0267129f20be3a2b8c8ff782f12701a7bd34b1fe7c228f517875906c2e5589c026ed89f2d474e0c3929743a644cdcccbc9567e32c8ee872d03cd77d9d38f4309dd2e5341dc32b04eaae63471d0763e85c4dab7104d0729495c15cc7b983406c4708b65ffc1ffff67ada77bab961cce25ffb4de4a349c81d6dbb35a5e495f8fad105ea3a2478826a70568f09a1cffa8935e29f90ae3be451bc3a2f53f4ac46d6510fc829c5db15d37ba1cb654ec3ab1544e95e451d35689252ee84096bfbd92afdd1afe7243d4555894bfcf7e5f382323f7052a7a98e1548c07955"
download.file(x,'y.jpg', mode = 'wb')

IMG1 <- load.image("y.jpg")

IMG.dat = data.frame()

IMG.dat[1,"NAME"] = paste0("IMG",1)
IMG.dat[1,"DIM1"] = dim(IMG1)[1]
IMG.dat[1,"DIM2"] = dim(IMG1)[2]
IMG.dat[1,"DIM3"] = dim(IMG1)[3]

# convert to RGB

img.rgb = data.frame(
  x = rep(1:IMG.dat[1,"DIM2"], each = IMG.dat[1,"DIM1"]),
  y = rep(IMG.dat[1,"DIM1"]:1, IMG.dat[1,"DIM2"]),
  R = as.vector(get(paste0(IMG.dat[1,"NAME"]))[,,1]),
  G = as.vector(get(paste0(IMG.dat[1,"NAME"]))[,,2]),
  B = as.vector(get(paste0(IMG.dat[1,"NAME"]))[,,3]))

IMG1.RGB = img.rgb

aalpha = c(25,25,2,2,0.5,0.5,0.5)

# use sampling in function to reduce MRI image size

WP.MRI = B_WP.IDX(x = IMG1.RGB[, c("R", "G", "B")], kmax = 8, corr = "pearson", method = "FCM", fzm = 2, sampling = 0.3, iter = 100,
             nstart = 20, NCstart = TRUE, alpha = aalpha, mult.alpha = 1/2)


pp = plot_BCVI(WP.MRI)
pp$plot_index
pp$plot_BCVI
pp$error_bar_plot

```



## License

The BayesCVI package as a whole is distributed under [GPL(>=3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
