B_CCV.IDX <- function(x, kmax, indexlist = "all", method = 'FCM', fzm = 2,
                    iter = 100, nstart = 20,
                    alpha = "default",
                    mult.alpha = 1/2){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(kmax))
    stop("Missing input argument. A maximum number of clusters  is required")
  if(!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if(kmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if(!any(indexlist %in% c("all","CCVP", "CCVS")))
    stop("Argument 'indexlist' is not in 'all', 'CCVP', 'CCVS'")
  if(!any(method  == c("FCM","EM")))
    stop("Argument 'method' should be one of 'FCM','EM' ")
  if(method == "FCM"){
    if(fzm <= 1)
      stop("Argument 'fcm' should be the number greater than 1",call. = FALSE)
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
    if(!is.numeric(iter))
      stop("Argument 'iter' must be numeric")
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
  # Defined vector
  ccvp = vector()
  ccvs = vector()
  distance =dist(x,diag = TRUE,upper= TRUE)
  # FOR CCVP CCVS
  distc = as.vector(as.matrix(distance))
  # start k loop
  for(k in kmin:kmax){
    if(method == "EM"){ # EM Algorithm
      EM.model <- Mclust(x,G=k,verbose=FALSE)
      assign("m",EM.model$z)
      assign("c",t(EM.model$parameters$mean))
    }else if(method == "FCM"){ # FCM Algorithm
      wd = Inf
      # cm.out = list()
      for (nr in 1:nstart){
        FCM.model = cmeans(x,k,iter,verbose=FALSE,method="cmeans",m=fzm)
        if (FCM.model$withinerror < wd){
          wd = FCM.model$withinerror
          FCM.model2 =FCM.model
        }
      }
      assign("m",FCM.model2$membership)
      assign("c",FCM.model2$centers)
    }
    uut = m%*%t(m)
    vnew = as.vector(1-(uut/max(uut)))

    if(sum(indexlist %in% c("all","CCVP"))>=1){
      ccvp[k-kmin+1] = cor(distc-mean(distc),vnew-mean(vnew),method = "pearson") #NW
    }
    if(sum(indexlist %in% c("all","CCVS"))>=1){
      ccvs[k-kmin+1] = cor(distc,vnew,method = "spearman") #NW
    }
  } # END CCVP CCVS index
  # Bayesian part
  if(any(indexlist %in% "all")){
    indexlist = c("CCVP","CCVS")
  }
  CCV.list = list()
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
    CCV.list[[paste0(indexlist[idx])]] = get(paste0(indexlist[idx],"_list"))
  }
  if (sum(indexlist %in% "all")>=1){
    return(CCV.list)
  } else {
    return(CCV.list[indexlist])
  }
}
