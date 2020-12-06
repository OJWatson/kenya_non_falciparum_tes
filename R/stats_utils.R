
#' @noRd
reorder_aic <- function(tab, x){
  tab$DeltaAIC <- tab$AIC-tab$AIC[which.min(tab$AIC)]
  tab$Weights <- round(exp(-0.5*tab$DeltaAIC) / sum(exp(-0.5*tab$DeltaAIC)),3)
  tab$AICc <- tab$AIC + ((2*tab$Params)*(tab$Params+1))/(x-tab$Params-1)
  tab <- tab[order(tab$AICc),]
  tab
}

#' @noRd
ll_mod_table <- function(res,models){

  ll_row <- function(x, mod){
    ll <- dmultinom(x$data, prob = x$params$multinom, log=TRUE)
    p <- length(x$params$params)
    aic = -2*ll + 2*p
    data.frame(Model = mod,AIC = aic,LogLik = ll,Params = p,stringsAsFactors = F)
  }


  tab <- purrr::map2_df(res,models,ll_row)
  tab <- reorder_aic(tab,sum(res[[1]]$data))

  return(tab)

}
#' @noRd
chi_val <- function(x,y){
  -2*x - -2*y
}
