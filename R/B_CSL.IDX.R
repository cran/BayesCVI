B_CSL.IDX <- function(x, kmax,
                    method = 'kmeans',
                    nstart = 100,
                    alpha = "default",
                    mult.alpha = 1/2){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(kmax))
    stop("Missing input argument. A maximum number of clusters is required")
  if(!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if(kmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if(!any(method  == c("kmeans","hclust_complete","hclust_average","hclust_single")))
    stop("Argument 'method' should be one of 'kmeans', 'hclust_complete', 'hclust_average', 'hclust_single'")
  if(method == "kmeans"){
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
  }
  if(startsWith(method,"hclust_")){
    H.model = hclust(dist(x),method = sub("hclust_", "", method))
  }
  if(!is.numeric(mult.alpha))
    stop("Argument 'mult.alpha' must be numeric")
  n = nrow(x)
  kmin = 2 #fix value
  if(any(alpha %in% "default")){
    alpha = rep(1,length(kmin:kmax))
  }
  if(length(kmin:kmax) != length(alpha)) # check
    stop("The length of kmin to kmax must be equal to the length of alpha")
  adj.alpha = alpha*(n)^mult.alpha
  # index part
  dm = dim(x)
  csl = rep(0,kmax-kmin+1)
  for(k in kmin:kmax){
    centroid = matrix(0,k,dm[2])
    if(method == "kmeans"){
      K.model = kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
      centroid = K.model$centers
      for(j in 1:k){
        dd = as.matrix(dist(x[cluss==j,]))
        csl[k-kmin+1] = csl[k-kmin+1] + sum(apply(dd,2,max))/sum(cluss==j)
      }
    } else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
      for (j in 1:k){
        if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
          centroid[j,] = as.numeric(x[cluss==j,])
        } else {
          centroid[j,] = colMeans(x[cluss==j,])
        }
        dd = as.matrix(dist(x[cluss==j,]))
        csl[k-kmin+1] = csl[k-kmin+1] + sum(apply(dd,2,max))/sum(cluss==j)
      }
    }
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    dd2 = as.matrix(dist(centroid))
    dd2 = matrix(dd2[dd2 > 0],k-1,k)
    csl[k-kmin+1] = csl[k-kmin+1]/sum(apply(dd2,2,min))
  }
  CVI.dframe = data.frame("C" = kmin:kmax,"Index" = csl)
  maxGI = max(CVI.dframe[,"Index"]) # The smallest value of the GI indicates the optimal number of cluster
  rk = (maxGI - CVI.dframe[,"Index"])/sum(maxGI - CVI.dframe[,"Index"])
  nrk = n*rk
  ex = (adj.alpha + nrk) / (sum(adj.alpha)+ n)
  var = ((adj.alpha+nrk)*(sum(adj.alpha)+n - adj.alpha - nrk))/((sum(adj.alpha)+n)^2*(sum(adj.alpha)+n+1))
  BCVI = data.frame("k" = kmin:kmax,"BCVI" = ex)
  VarBCVI = data.frame("k" = kmin:kmax,"Var" = var)
  colnames(CVI.dframe) = c("k","CSL")
  CSL.result = list("BCVI" = BCVI,"VAR" = VarBCVI,"Index" = CVI.dframe)
  return(CSL.result)
}
