library(icer)
library(tidyverse)
load(file.path(here::here(),"analysis/data/derived/spec_icer.rda"))

# Explore the model fits over the parameter space
brackets <- function(max=100,
                     num_brackets=20,
                     geometric_brackets=TRUE){

  if(geometric_brackets){
    ## Create the geometric age brackets
    ratio <- (max/0.1)^(1/num_brackets)
    vec <- 0.1 * ratio ** (1:num_brackets)
    vec[1] <- 0
  } else {
    vec <- seq(0,max,num_brackets)
  }
  return(vec)
}

# grid spacing
mu <- seq(1,5,length.out = 29)
size <- brackets(1000,30)[-1]

# parameter grid to scan over
sweep <- expand.grid(mu,size,KEEP.OUT.ATTRS = FALSE)
ll <- rep(0, nrow(sweep))
for(i in 1:nrow(sweep)) {
  message(i)
  p <- spec_res$params$params
  p[5:6] <- as.numeric(sweep[i,])
  m <- icer:::independent(p, spec_res$data, pci_list = spec_res$pci_list, max_moi = 25)
  ll[i] <- dmultinom(spec_res$data, prob = m, log = TRUE)

}
z <- ll %>% matrix(ncol=length(mu),byrow = TRUE) * -2

# turn into something that can be plotted
df <- sweep
df$LogLik <- as.numeric(ll)
names(df)[1:2] <- c("mu","size")
df <- as.data.frame(df)
closest <- function(x,y) { which.min(abs(x-y)) }

# heatmap for plotting
heatmap <- ggplot(df,mapping = aes(x=mu,y=size,fill=LogLik,z=LogLik)) + ggplot2::geom_tile() +
  scale_y_log10(expand=c(0,0)) +
  ylab("r") +
  xlab(expression(~italic(mu))) +
  geom_contour(color="black") +
  metR::geom_text_contour(stroke=0.4,stroke.color = viridis::viridis(841)[sapply(c(-93,-70,-50,-50,-36),closest,sort(df$LogLik))],size=5) +
  viridis::scale_fill_viridis(name="Log Likelihood") + coord_equal() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size=20,face = "italic"),
        text = element_text(size=14))


save_figs("r_mu_heatmap", heatmap, height = 8, width = 8)
