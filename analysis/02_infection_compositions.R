
# ------------------------------------------------------------------------------
# Running infection compositon models - N.B. Uses Imperial Cluster set up
# This can be run locally (see comments re lapply) but the resuts are saved
# in the data directory and read in in the next code chunk ("Loading icer runs")
# ------------------------------------------------------------------------------

# Running icer runs on Imperial cluster
contextualise <- function(workdir = "M:/OJ/icer_results",
                          packages = c("icer"),
                          linux = TRUE,
                          package_sources = provisionr::package_sources(
                            local=getwd()),
                          new_dir = "tests",
                          use_workers=TRUE,
                          sources=NULL,
                          cluster = "fi--dideclusthn",
                          cores=NULL,
                          rtools=TRUE,
                          template=NULL,
                          initialise=TRUE){

  if (is.null(template)){
    template <- "GeneralNodes"
  }

  dir.create(paste0(workdir,"/",new_dir),showWarnings = FALSE)
  if(linux) {

    workdir <- gsub("M:/OJ","/home/oj/net/Malaria/OJ",workdir)
    if(grepl("nas", workdir)){
      home <- "//fi--didenas1/Malaria"
      local_home <- "/home/oj/net/nas"
    } else{
      home <- "//fi--didef3.dide.ic.ac.uk/Malaria"
      local_home <- "/home/oj/net/Malaria"
    }
    didehpc::didehpc_config_global(workdir=workdir,
                                   credentials="/home/oj/.smbcredentials",
                                   temp=didehpc::path_mapping("tmp",
                                                              "/home/oj/net/temp",
                                                              "//fi--didef3.dide.ic.ac.uk/tmp",
                                                              "T:"),
                                   home=didehpc::path_mapping("OJ",
                                                              local_home,
                                                              home,
                                                              "M:"),
                                   cluster = cluster)
  } else {

    didehpc::didehpc_config_global(workdir=workdir,
                                   credentials="C:\\Users\\Oliver\\.smbcredentials",
                                   temp=didehpc::path_mapping("tmp",
                                                              "T:",
                                                              "//fi--didef3.dide.ic.ac.uk/tmp",
                                                              "T:"),
                                   home=didehpc::path_mapping("OJ",
                                                              "M:",
                                                              "//fi--didef3.dide.ic.ac.uk/Malaria",
                                                              "M:"),
                                   cluster = cluster)
  }

  didehpc::web_login()

  root <- file.path(workdir, new_dir)
  packages.vector <- packages

  context::context_log_start()
  ## set up context
  ctx <- context::context_save(root,sources = sources,
                               packages = packages.vector,
                               package_sources = package_sources)

  if(is.null(cores)){
    config <- didehpc::didehpc_config(use_workers = use_workers, rtools = rtools,
                                      template = template)
  } else {
    config <- didehpc::didehpc_config(use_workers = use_workers, rtools = rtools,
                                      template = template,cores = cores)
  }
  config$resource$parallel <- "FALSE"
  config$resource$type <- "Cores"

  obj <- didehpc::queue_didehpc(ctx, config = config,initialise = initialise)
  return(obj)
}
obj <- contextualise(workdir = "M:/OJ/icer_results",
                     new_dir = "results",
                     packages="icer",
                     use_workers = FALSE,
                     cluster = "big",
                     initialise = TRUE)

# Kenya data
df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))
spec <- table(df$species)

# Use the icer package written for the analysis (www.github.com/OJWatson/icer)
library(icer)

# create a parameter list for running all the models
pl_list <- list()

# initial starting conditions
t <- c(0.9, 0.5, 2, 1.2, 2, 1.8)
names(t) <- c("k_12", "k_13", "k_14", "k_23", "k_24", "k_34")
co <- 1
for(j in 1:5){
  c <- combinat::combn(1:6,j)
  for(k in 1:ncol(c)){
    pl_list[[co]] <- t[c[,k]]
    co <- co+1
  }
}

# Running on our cluster we pass in the species composition numbers into the following function calls

# Run all our Poisson Models on cluster (to Replicate locally simply call lapply(pl_list, FUN = ...))
grp_62_pois <- queuer::qlapply(pl_list,obj=obj,
                               FUN=function(x){
                                 spec <- c("1469","71","5","13","3","60","345","20","6","1","3","10","19","2")
                                 spec <- as.numeric(spec)
                                 names(spec) <- c("Pf","Pf/Pm","Pf/Pm/PoC","Pf/Pm/PoW","Pf/Pm/PoW/PoC","Pf/PoC",
                                                 "Pf/PoW","Pf/PoW/PoC","Pm","Pm/PoC","Pm/PoW","PoC","PoW","PoW/PoC")
                                 cooccurence_test(data = spec,  boot_iter = 50000, poisson=TRUE,
                                                  plot = TRUE,size =100,
                                                  quantiles = c(0.025,0.975),
                                                  lower = NULL,upper=NULL,
                                                  density_func = icer:::interference,
                                                  max_moi = 25,x)
                               },
                               name = "spec_62_pois",overwrite = TRUE)

# Run all our Negbinom Model on cluster (to Replicate locally simply call lapply(pl_list, FUN = ...))
grp_62_nb <- queuer::qlapply(pl_list,obj=obj,
                             FUN=function(x){
                               spec <- c("1469","71","5","13","3","60","345","20","6","1","3","10","19","2")
                               spec <- as.numeric(spec)
                               names(spec) <- c("Pf","Pf/Pm","Pf/Pm/PoC","Pf/Pm/PoW","Pf/Pm/PoW/PoC","Pf/PoC",
                                               "Pf/PoW","Pf/PoW/PoC","Pm","Pm/PoC","Pm/PoW","PoC","PoW","PoW/PoC")
                               cooccurence_test(data = spec,  boot_iter = 50000, poisson=FALSE,
                                                plot = TRUE,size = 100,
                                                quantiles = c(0.025,0.975),
                                                lower = NULL,upper=NULL,
                                                density_func = icer:::interference,
                                                max_moi = 25,x)
                             },
                             name = "spec_62_nb",overwrite = TRUE)

###
# Grab results from the cluster
l <- lapply(grp_62_pois$tasks,function(x){x$result()})
ln <- lapply(grp_62_nb$tasks,function(x){x$result()})

# remove plot for memory reasons
l <- lapply(l, function(x) {x$plot <- NULL; return(x)})
ln <- lapply(ln, function(x) {x$plot <- NULL; return(x)})

save(l,ln,file="analysis/data/derived/complete_icer.rda")

###
# Run the analysis for the complete independent and inteference models
spec <- table(df$species)
spec_res_poisson <- cooccurence_test(spec, boot_iter = 10000,poisson = TRUE)
spec_res <- cooccurence_test(spec, boot_iter = 10000)
spec_res_pois_3p_interference <- cooccurence_test(spec,boot_iter = 10000,poisson=TRUE,
                                                 density_func = icer:::interference,
                                                 k_12=0.8, k_13=0.75, k_14=1.5, k_23=1, k_24=1.5, k_34=1.2)
spec_res_3p_interference <- cooccurence_test(spec,boot_iter = 10000,density_func = icer:::interference,
                                            k_12=0.8, k_13=0.75, k_14=1.5, k_23=1, k_24=1.5, k_34=1.2)

# save these results to be read in easily
save(spec_res, spec_res_poisson, spec_res_pois_3p_interference, spec_res_3p_interference,
     file = file.path(here::here(),"analysis/data/derived/spec_icer.rda"))

# ------------------------------------------------------------------------------
# Reading in infection model results (Loading icer runs)
# ------------------------------------------------------------------------------

load(file.path(here::here(),"analysis/data/derived/complete_icer.rda"))
load(file.path(here::here(),"analysis/data/derived/spec_icer.rda"))

# Kenya data
df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))
spec <- table(df$species)

# Create our table of model results
res <- list(spec_res_poisson,spec_res,spec_res_pois_3p_interference,spec_res_3p_interference)
models <- c("Independent","Independent","Complete Interference","Complete Interference")
tab <- ll_mod_table(res,models)
tab$Dist <- c("Poisson","Poisson","Negative Binomial","Negative Binomial")

# create a parameter list for running all the models
pl_list <- list()

# pl list used in the full model createion
t <- c(0.9, 0.5, 2, 1.2, 2, 1.8)
names(t) <- c("k_12", "k_13", "k_14", "k_23", "k_24", "k_34")
co <- 1
for(j in 1:5){
  c <- combinat::combn(1:6,j)
  for(k in 1:ncol(c)){
    pl_list[[co]] <- t[c[,k]]
    co <- co+1
  }
}

# And table for the complete parameter scan
pt <- ll_mod_table(l,lapply(pl_list,function(x){paste0(names(x),collapse="/")}))
pt$Dist <- "Poisson"
ptn <- ll_mod_table(ln,lapply(pl_list,function(x){paste0(names(x),collapse="/")}))
ptn$Dist <- "Negative Binomial"

# Pull it all together
full <- rbind(pt,ptn,tab)
full <- full[order(full$AICc),c(1,8,2:7)]
full <- reorder_aic(full,2027)
full$Model[grep("k",full$Model)] <- paste(full$Model[grep("k",full$Model)],"Interference")
write.csv(full[,c("Model", "Dist", "AICc", "AIC", "DeltaAIC", "Params", "LogLik")],
          file.path(here::here(), "analysis/tables/complete_aic_table.csv"))

# And the smaller best tables
best <- full[c(1,which(full$AICc==min(full$AICc[full$Dist=="Negative Binomial"])),
               which(full$Model %in% c("Independent","Complete Interference"))),]
best <- reorder_aic(best, 2027)
best_mods <-  list(l[[2]],ln[[2]],  # lowest AIC model
                   spec_res_poisson,
                   spec_res_pois_3p_interference,
                   spec_res,
                   spec_res_3p_interference)
best_table <- cbind(best, rbind_list_base(lapply(best_mods,function(x){

  pfs <- as.data.frame(round(t(x$params$params[1:5]),3))
  if(length(x$params$params) == 5) {
    pfs$Extra <- ""
  } else {
    ex_p <- 6:length(x$params$params)
    extra <- paste0(mapply(paste,names(x$params$params[ex_p]),round(x$params$params[ex_p],3),MoreArgs = list(sep="=")),collapse=", ")
    pfs$Extra <- extra
  }
  return(pfs)
})))

best_table <- best_table[,c(c("Model", "Dist","AICc","AIC", "DeltaAIC", "pf", "pm", "poc", "pow", "mu","Extra"))]
best_table[,c(3,4,5)] <- round(best_table[,c(3,4,5)],3)
write.csv(best_table, file.path(here::here(), "analysis/tables/best_aic_table.csv"))

# ------------------------------------------------------------------------------
# Figure 2
# ------------------------------------------------------------------------------

# Create our plots from our best runs
a <- icer:::coinf_plot(reps = 50000,probs = spec_res_poisson$params$multinom ,
                levels = names(spec_res_poisson$params$multinom),
                total = sum(spec_res_poisson$data),
                plot = TRUE, real = spec_res$data,
                density = FALSE)

b <- icer:::coinf_plot(reps = 50000,probs = l[[2]]$params$multinom ,
                levels = names(l[[2]]$params$multinom),  # lowest AIC model
                total = sum(l[[2]]$data),  # lowest AIC model
                plot = TRUE, real = spec_res$data,
                density = FALSE)

# Format them nicely
scaleFUN <- function(x) sprintf("%.2f", x)

a_plot <- a$plot + ggplot2::ylab("Density") +
  scale_x_continuous(labels= scales::number_format(accuracy = 1)) +
  theme(text = element_text(size=16))

b_plot <- b$plot + ggplot2::ylab("Density") +
  scale_x_continuous(labels= scales::number_format(accuracy = 1)) +
  theme(text = element_text(size=16))

# Save figure 2
gp <- cowplot::plot_grid(a_plot,
                         NULL,
                         b_plot,
                         labels = c("a","","b"),
                         ncol=1, rel_heights = c(1,0.05,1))

save_figs("fig2", gp, width = 16, height = 16)

