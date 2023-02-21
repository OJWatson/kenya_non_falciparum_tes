
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

# Setting Up Cluster From New

# Log in to didehpc
credentials = "C:/Users/ow813/.smbcredentials"
options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "ow813")

## ------------------------------------
## 1. Package Installations
## ------------------------------------
# drat:::add("mrc-ide")
# install.packages("pkgdepends")
# install.packages("didehpc")
# install.packages("orderly")

## ------------------------------------
## 2. Setting up a cluster configuration
## ------------------------------------

options(didehpc.cluster = "fi--didemrchnb")

# not if T is not mapped then map network drive
didehpc::didehpc_config_global(temp=didehpc::path_mapping("tmp",
                                                          "T:",
                                                          "//fi--didef3.dide.ic.ac.uk/tmp",
                                                          "T:"),
                               home=didehpc::path_mapping("OJ",
                                                          "L:",
                                                          "//fi--didenas5/malaria",
                                                          "L:"),
                               credentials=credentials,
                               cluster = "fi--didemrchnb")

# Creating a Context
context_name <- cp_path("context")

ctx <- context::context_save(
  path = context_name,
  package_sources = conan::conan_sources(
    packages = c(
      "OJWatson/icer"
    )
  )
)

# set up a specific config for here as we need to specify the large RAM nodes
config <- didehpc::didehpc_config(template = "24Core")
config$resource$parallel <- "FALSE"
config$resource$type <- "Cores"

# Configure the Queue
obj <- didehpc::queue_didehpc(ctx, config = config)

library(icer)

pl_list <- list()
co <- 1
t <- tail(icers$intctn_nb$params$params,6)
for(j in 1:5){
  c <- combinat::combn(1:6,j)
  for(k in 1:ncol(c)){
    pl_list[[co]] <- list()
    pl_list[[co]]$x <- t[c[,k]]
    pl_list[[co]]$data = icer_dat
    co <- co+1
  }
}

grp_62_pois <- obj$lapply(pl_list,
                               FUN=function(x){
                                 icer::cooccurence_test(data = x$data,  boot_iter = 50000, poisson=TRUE,
                                                  plot = TRUE,size =100,
                                                  quantiles = c(0.025,0.975),
                                                  lower = 0.01,upper=10,
                                                  density_func = icer:::interference,
                                                  max_moi = 25,x$x)
                               },
                               name = "hoseah_62_pois_nftes",overwrite = TRUE)

grp_62_nb <- obj$lapply(pl_list,
                        FUN=function(x){
                          icer::cooccurence_test(data = x$data,  boot_iter = 50000, poisson=FALSE,
                                           plot = TRUE,size =100,
                                           quantiles = c(0.025,0.975),
                                           lower = 0.01,upper=10,
                                           density_func = icer:::interference,
                                           max_moi = 25,x$x)
                        },
                        name = "hoseah_62_nb_nftes",overwrite = TRUE)

###

l <- lapply(grp_62_pois$tasks,function(x){x$result()})
pt <- ll_mod_table(l,lapply(pl_list,function(x){paste0(names(x$x),collapse="/")}))
pt$Dist <- "Poisson"
pt$run <- rownames(pt)

ln <- lapply(grp_62_nb$tasks,function(x){x$result()})
ptn <- ll_mod_table(ln,lapply(pl_list,function(x){paste0(names(x$x),collapse="/")}))
ptn$Dist <- "Negative Binomial"
ptn$run <- rownames(ptn)

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
tab$run <- 1:4

full <- rbind(pt,ptn,tab)
full <- full[order(full$AICc),c(1,8,2:7, 9)]
full <- reorder_aic(full,161)
full$Model[grep("k",full$Model)] <- paste(full$Model[grep("k",full$Model)],"Inteference")
write.csv(full[,1:8], "analysis/tables/complete_aic_table.csv")

if(full$Dist[1] == "Poisson"){
  best_fit <- l[[as.numeric(full$run[1])]]
} else {
  best_fit <- ln[[as.numeric(full$run[1])]]
}
saveRDS(best_fit, "analysis/data/derived/best_icer.rds")
best_fit <- readRDS("analysis/data/derived/best_icer.rds")
icers <- readRDS(cp_path("analysis/data/derived/icer.rds"))

# comparison between best and independent

gga <- icer:::dens_plot(x = best_fit$plot$vals,
                       real = best_fit$data,
                       levels = names(best_fit$params$multinom),
          density = FALSE)

ggb <- icer:::dens_plot(x = icers$indpnt_pois$plot$vals,
                        real = best_fit$data,
                        levels = names(icers$indpnt_pois$params$multinom),
                        density = FALSE)


gg <- cowplot::plot_grid(ggb, gga, ncol = 1, labels ="auto")
save_figs("icer_plots", gg, width = 10, height = 10)
