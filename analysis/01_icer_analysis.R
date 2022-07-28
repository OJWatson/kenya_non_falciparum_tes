
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
icer_dat <- daily %>% filter(day == 0 & infx != "Uninfected" & infx != "Missing") %>% pull(infx) %>% table
icer_dat <- icer_dat[!names(icer_dat) %in% c("Uninfected", "Missing")]

# maximum models to run locally
icer_res_poisson <- icer::cooccurence_test(icer_dat, max_moi = 10, boot_iter = 10000,poisson = TRUE)
icer_res <- icer::cooccurence_test(icer_dat, max_moi = 10,boot_iter = 10000)
icer_res_pois_3p_interference <- icer::cooccurence_test(icer_dat,boot_iter = 10000,poisson=TRUE,
                                                 density_func = icer:::interference,max_moi = 10,
                                                 k_12=0.8, k_13=0.75, k_14=1.5, k_23=1, k_24=1.5, k_34=1.2)
icer_res_3p_interference <- icer::cooccurence_test(icer_dat,boot_iter = 10000,density_func = icer:::interference,max_moi = 10,
                                            k_12=0.8, k_13=0.75, k_14=1.5, k_23=1, k_24=1.5, k_34=1.2)

icers <- list("indpnt_pois" = icer_res_poisson,
     "indpnt_nb" = icer_res,
     "intctn_pois" = icer_res_pois_3p_interference,
     "intctn_nb" = icer_res_3p_interference)

saveRDS(icers, cp_path("analysis/data/derived/icer.rds"))
icers <- readRDS(cp_path("analysis/data/derived/icer.rds"))

# cluster work - this is specific for Imperial DIDE Cluster.
# To run locally convert qlapply below to simple lapply and remove obj, name and overwrite args
################################

obj <- roxer::contextualise(workdir = "M:/OJ/icer_results",new_dir = "results",packages="icer",use_workers = FALSE,cluster = "big",initialise = TRUE)


# hoseah's data
load("analysis/data/derived/data.rda")

library(icer)

pl_list <- list()
co <- 1
t <- tail(icers$intctn_nb$params$params,6)
for(j in 1:5){
  c <- combinat::combn(1:6,j)
  for(k in 1:ncol(c)){
    pl_list[[co]] <- t[c[,k]]
    co <- co+1
  }
}


grp_62_pois <- queuer::qlapply(pl_list,obj=obj,
                               FUN=function(x){
                                 hos <- c("1469","71","5","13","3","60","345","20","6","1","3","10","19","2")
                                 hos <- as.numeric(hos)
                                 names(hos) <- c("Pf","Pf/Pm","Pf/Pm/PoC","Pf/Pm/PoW","Pf/Pm/PoW/PoC","Pf/PoC",
                                                 "Pf/PoW","Pf/PoW/PoC","Pm","Pm/PoC","Pm/PoW","PoC","PoW","PoW/PoC")
                                 cooccurence_test(data = hos,  boot_iter = 50000, poisson=TRUE,
                                                  plot = TRUE,size =100,
                                                  quantiles = c(0.025,0.975),
                                                  lower = NULL,upper=NULL,
                                                  density_func = icer:::interference,
                                                  max_moi = 25,x)
                               },
                               name = "hoseah_62_pois_tes",overwrite = TRUE)

grp_62_nb <- queuer::qlapply(pl_list,obj=obj,
                             FUN=function(x){
                               hos <- c("1469","71","5","13","3","60","345","20","6","1","3","10","19","2")
                               hos <- as.numeric(hos)
                               names(hos) <- c("Pf","Pf/Pm","Pf/Pm/PoC","Pf/Pm/PoW","Pf/Pm/PoW/PoC","Pf/PoC",
                                               "Pf/PoW","Pf/PoW/PoC","Pm","Pm/PoC","Pm/PoW","PoC","PoW","PoW/PoC")
                               cooccurence_test(data = hos,  boot_iter = 50000, poisson=FALSE,
                                                plot = TRUE,size = 100,
                                                quantiles = c(0.025,0.975),
                                                lower = NULL,upper=NULL,
                                                density_func = icer:::interference,
                                                max_moi = 25,x)
                             },
                             name = "hoseah_62_nb_tes",overwrite = TRUE)

###

l <- lapply(grp_62_pois$tasks,function(x){x$result()})
pt <- icer:::ll_mod_table(l,lapply(pl_list,function(x){paste0(names(x),collapse="/")}))
pt$Dist <- "Poisson"

ln <- lapply(grp_62_nb$tasks,function(x){x$result()})
ptn <- ll_mod_table(ln,lapply(pl_list,function(x){paste0(names(x),collapse="/")}))
ptn$Dist <- "Negative Binomial"
save(l,ln,file="analysis/data/derived/complete_icer.rda")
load("analysis/data/derived/complete_icer.rda")

#############

reorder_aic <- function(tab, x){
  tab$DeltaAIC <- tab$AIC-tab$AIC[which.min(tab$AIC)]
  tab$Weights <- round(exp(-0.5*tab$DeltaAIC) / sum(exp(-0.5*tab$DeltaAIC)),3)
  tab$AICc <- tab$AIC + ((2*tab$Params)*(tab$Params+1))/(x-tab$Params-1)
  tab <- tab[order(tab$AICc),]
  tab
}

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

chi_val <- function(x,y){
  -2*x - -2*y
}

res <- list(icers[[1]],icers[[2]],icers[[3]],icers[[4]])
models <- c("Independent","Independent","Complete Interference","Complete Interference")
tab <- ll_mod_table(res,models)
tab$Dist <- c("Poisson","Poisson","Negative Binomial","Negative Binomial")

full <- rbind(pt,ptn,tab)
full <- full[order(full$AICc),c(1,8,2:7)]
full <- reorder_aic(full,2027)
full$Model[grep("k",full$Model)] <- paste(full$Model[grep("k",full$Model)],"Inteference")
write.csv(full[,1:8], "analysis/data/derived/complete_aic_table.csv")
