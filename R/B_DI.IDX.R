B_DI.IDX <- function(x, kmax,
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
  di = vector()
  for(k in kmin:kmax){
    if(method == "kmeans"){
      K.model = kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
    } else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
    } # End check algorithm
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    size = table(cluss)
    dunn.dem = 0
    dunn.num = 1e10
    for (i in 1:(k-1)){
      dunn.dem = max(dunn.dem,max(dist(x[cluss==i,])))
      for(j in (i+1):k){
        dunn.num = min(dunn.num,min(as.matrix(dist(rbind(x[cluss==i,],x[cluss==j,])))[(size[i]+size[j]):(size[i]+1),1:size[i]]))
      }
    }
    dunn.dem = max(dunn.dem,max(dist(x[cluss==k,])))
    di[k-kmin+1] = dunn.num / dunn.dem
  }
  # Bayesian
  CVI.dframe = data.frame("C" = kmin:kmax,"Index" = di)
  minGI = min(CVI.dframe[,"Index"]) # The largest value of the GI indicates the optimal number of cluster
  rk = (CVI.dframe[,"Index"] - minGI)/sum(CVI.dframe[,"Index"] - minGI)
  nrk = n*rk
  ex = (adj.alpha + nrk) / (sum(adj.alpha)+ n)
  var = ((adj.alpha+nrk)*(sum(adj.alpha)+n - adj.alpha - nrk))/((sum(adj.alpha)+n)^2*(sum(adj.alpha)+n+1))
  BCVI = data.frame("k" = kmin:kmax,"BCVI" = ex)
  VarBCVI = data.frame("k" = kmin:kmax,"Var" = var)
  colnames(CVI.dframe) = c("k","DI")
  DI.result = list("BCVI" = BCVI,"VAR" = VarBCVI,"Index" = CVI.dframe)
  return(DI.result)
}
