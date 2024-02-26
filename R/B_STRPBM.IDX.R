B_STRPBM.IDX <- function(x, kmax,
                       method = 'kmeans',
                       indexlist = 'all',  #c(,"all","STR","PBM")
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
  if(!any(indexlist %in% c("all","STR","PBM")))
    stop("Argument 'indexlist' should be 'all', 'STR', 'PBM'")
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
  str = rep(0,kmax-kmin+1)
  pbm = rep(0,kmax-kmin+1)
  EK = rep(0,kmax-kmin+2)
  DK = rep(0,kmax-kmin+3)
  md = rep(0,kmax-kmin+3)

  if (kmin == 2){
    lb = 2
  } else {
    lb = kmin-1
  }

  for(k in lb:(kmax+1)){
    xnew = matrix(0,dm[1],dm[2])
    centroid = matrix(0,k,dm[2])
    if(method == "kmeans"){
      K.model = kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
      centroid = K.model$centers
      xnew = centroid[cluss,]
    } else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
      for (j in 1:k){
        if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
          centroid[j,] = as.numeric(x[cluss==j,])
        } else {
          centroid[j,] = colMeans(x[cluss==j,])
        }
      }
      xnew = centroid[cluss,]
    } # End check algorithm
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    EK[k-kmin+2] = sum(sqrt(rowSums((x - xnew)^2)))
    ddd = dist(centroid)
    md[k-kmin+2] = max(ddd)
    DK[k-kmin+2] = max(ddd)/min(ddd)
  }
  E0 = sum(sqrt(rowSums((x-colMeans(x))^2)))
  if (kmin == 2){
    EK[1] = E0
  }
  EKK = E0/EK
  str = (EKK[2:(length(EKK)-1)]-EKK[1:(length(EKK)-2)])*(DK[3:(length(DK))]-DK[2:(length(DK)-1)])
  pbm = EKK[2:(length(EKK)-1)]*md[2:(length(EKK)-1)]/(kmin:kmax)
  # Bayesian part
  if(any(indexlist %in% "all")){
    indexlist = c("STR","PBM")
  }
  STR.list = list()
  for (idx in seq(length(indexlist))) {
    CVI.dframe = data.frame("C" = kmin:kmax,"Index" = get(tolower(indexlist[idx])))
    minGI = min(CVI.dframe[,"Index"]) # The largest value of the GI indicates the optimal number of cluster
    rk = (CVI.dframe[,"Index"] - minGI)/sum(CVI.dframe[,"Index"] - minGI)
    nrk = n*rk
    ex = (adj.alpha + nrk) / (sum(adj.alpha)+ n)
    var = ((adj.alpha+nrk)*(sum(adj.alpha)+n - adj.alpha - nrk))/((sum(adj.alpha)+n)^2*(sum(adj.alpha)+n+1))
    BCVI = data.frame("k" = kmin:kmax,"BCVI" = ex)
    VarBCVI = data.frame("k" = kmin:kmax,"Var" = var)
    colnames(CVI.dframe) = c("k",paste0(indexlist[idx]))
    list.re = list("BCVI" = BCVI,"VAR" = VarBCVI,"Index" = CVI.dframe)
    assign(paste0(indexlist[idx],"_list"),list.re)
    STR.list[[paste0(indexlist[idx])]] = get(paste0(indexlist[idx],"_list"))
  }
  if (sum(indexlist == "all")==1){
    return(STR.list)
  } else {
    return(STR.list[indexlist])
  }
}
