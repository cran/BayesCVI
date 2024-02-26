plot_BCVI <- function(B.result,
                      mult.err.bar = 2){
  if(sum(names(B.result) %in% c("BCVI","VAR","Index","opt.pt")) <=2)
    stop("Bad input, 'B.result' is not a result from a function in the package.")
  if(sum(names(B.result) %in% c("BCVI","VAR","Index","opt.pt")) ==4){
    opt.pt = B.result[["opt.pt"]]
  } else if((any(colnames(B.result[["Index"]][2]) %in% c("WP","PBM","KPBM","CCVP","CCVS","CH","DI","PB","NCI","STR")))){
    opt.pt = "max"
  } else{
    opt.pt = "min"
  }
  if(!is.numeric(mult.err.bar))
    stop("Argument 'mult.err.bar' must be numeric")
  # BCVI plot
  pp.BCVI = ggplot(B.result[["BCVI"]], aes(k, BCVI))+
    geom_point(size = 2)+
    geom_point(data = B.result[["BCVI"]][which.max(B.result[["BCVI"]][,2]),],colour = "red")+
    geom_line()+
    theme_minimal()
  pp = data.frame(B.result[["BCVI"]],"VAR" = mult.err.bar*sqrt(B.result[["VAR"]][,2]))
  pp.VAR = ggplot(pp, aes(k, BCVI)) +
    geom_line() +
    geom_errorbar(aes(ymin = BCVI-VAR, ymax = BCVI+VAR,width = 0.2)) +
    geom_errorbar(data = pp[which.max(pp[,"BCVI"]),],aes(ymin = BCVI-VAR, ymax = BCVI+VAR,width = 0.2),
                  colour = "red",size=1)+
    geom_point(size = 1)+
    geom_text(aes(label= round(BCVI,4)), colour="black", size =3,hjust=-0.3)+
    theme_minimal()
  pp.Index = ggplot(B.result[["Index"]], aes_string(x = colnames(B.result[["Index"]])[1],y =colnames(B.result[["Index"]])[2]))+
    geom_point(size = 2)+
    geom_point(data = if(opt.pt=="max"){
      B.result[["Index"]][which.max(B.result[["Index"]][,2]),]
    }else{
      B.result[["Index"]][which.min(B.result[["Index"]][,2]),]
    },colour = 'red')+
    geom_line()+
    theme_minimal()
 result = list("plot_BCVI" = pp.BCVI,"error_bar_plot" = pp.VAR,"plot_index" = pp.Index)
}
