B_DB.IDX <- function(x, kmax,
                   method = 'kmeans',
                   indexlist = 'all', #c(,"all","DB","DBs")
                   p = 2, q = 2,
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
  if(!any(indexlist %in% c("all","DB","DBs")))
    stop("Argument 'indexlist' should be 'all', 'DB', 'DBs'")
  if(!is.numeric(p))
    stop("Argument 'p' must be numeric")
  if(!is.numeric(q))
    stop("Argument 'q' must be numeric")
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
  db = vector()
  dbs = vector()
  for(k in kmin:kmax){
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
    # Si,q
    S = vector() #length = k
    sizecluss  = as.vector(table(cluss))
    for(i in 1:k){
      C = sizecluss[i]
      if(C>1){
        cenI = xnew[cluss ==i,]
        S[i] = (sum(sqrt(rowSums((x[cluss==i,]-cenI)^2))^q)/C)^(1/q)
      }else{
        S[i] = 0
      }
    }
    m = as.matrix(dist(centroid,method="minkowski",p=p))
    R = matrix(0,k,k)
    r = vector()
    rs = vector()
    wcdd = vector()
    for (i in 1:k){
      C = sizecluss[i]
      r[i] = max((S[i] + S[-i])/m[i,][m[i,]!=0])
      rs[i] = max(S[i] + S[-i])/min(m[i,][m[i,]!=0])
      # for SF index
      wcdd[i] = sum(dist(rbind(centroid[i,],x[cluss==i,]))[1:C])/C
    }
    db[k-kmin+1] = mean(r)
    dbs[k-kmin+1] = mean(rs)
  }
  if(any(indexlist %in% "all")){
    indexlist = c("DB","DBs")
  }
  DB.list = list()
  for (idx in seq(length(indexlist))) {
    CVI.dframe = data.frame("C" = kmin:kmax,"Index" = get(tolower(indexlist[idx])))
    maxGI = max(CVI.dframe[,"Index"]) # The smallest value of the GI indicates the optimal number of cluster
    rk = (maxGI - CVI.dframe[,"Index"])/sum(maxGI - CVI.dframe[,"Index"])
    nrk = n*rk
    ex = (adj.alpha + nrk) / (sum(adj.alpha)+ n)
    var = ((adj.alpha+nrk)*(sum(adj.alpha)+n - adj.alpha - nrk))/((sum(adj.alpha)+n)^2*(sum(adj.alpha)+n+1))
    BCVI = data.frame("k" = kmin:kmax,"BCVI" = ex)
    VarBCVI = data.frame("k" = kmin:kmax,"Var" = var)
    colnames(CVI.dframe) = c("k",paste0(indexlist[idx]))
    list.re = list("BCVI" = BCVI,"VAR" = VarBCVI,"Index" = CVI.dframe)
    assign(paste0(indexlist[idx],"_list"),list.re)
    DB.list[[paste0(indexlist[idx])]] = get(paste0(indexlist[idx],"_list"))
  }
  if (sum(indexlist == "all")==1){
    return(DB.list)
  } else {
    return(DB.list[indexlist])
  }
}
