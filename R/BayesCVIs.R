BayesCVIs <- function(CVI,
                     n,
                      kmax,
                      opt.pt,
                      alpha = "default",
                      mult.alpha = 1/2){
  if(missing(CVI))
    stop("Missing input argument. A numeric result from cluster validity index is required")
  if(missing(n))
    stop("Missing input argument. A numeric data point value is required")
  if(!is.numeric(n))
    stop("Argument 'n' must be numeric")
  if(missing(kmax))
    stop("A maximum number of clusters from cluster validity index is required")
  if(!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if(!opt.pt %in% c("min","max"))
    stop("Argument 'opt.pt' should be one of 'min', 'max'")
  if(!is.numeric(mult.alpha))
    stop("Argument 'mult.alpha' must be numeric")
  kmin = 2
  if(any(alpha %in% "default")){
    alpha = rep(1,length(kmin:kmax))
  }
  adj.alpha = alpha*(n)^mult.alpha
  if(length(kmin:kmax) != length(adj.alpha))
    stop("The length of kmin to kmax must be equal to the length of alpha")
  # create new dataframe for index
  CVI.dframe = data.frame("C" = kmin:kmax,"Index" = CVI)
  # Bayesian CVI
  if(opt.pt == "max"){
    minGI = min(CVI.dframe[,"Index"])
    rk = (CVI.dframe[,"Index"] - minGI)/sum(CVI.dframe[,"Index"] - minGI)
  } else {
    maxGI = max(CVI.dframe[,"Index"])
    rk = (maxGI - CVI.dframe[,"Index"])/sum(maxGI - CVI.dframe[,"Index"])
  }
  nrk = n*rk
  ex = (adj.alpha + nrk) / (sum(adj.alpha)+ n)
  var = ((adj.alpha+nrk)*(sum(adj.alpha)+n - adj.alpha - nrk))/((sum(adj.alpha)+n)^2*(sum(adj.alpha)+n+1))
  BCVI = data.frame("k" = kmin:kmax,"BCVI" = ex)
  VarBCVI = data.frame("k" = kmin:kmax,"Var" = var)
  result = list("BCVI" = BCVI,"VAR" = VarBCVI,"Index" = CVI.dframe,"opt.pt" = opt.pt)
  return(result)
}
