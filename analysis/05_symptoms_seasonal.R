# ------------------------------------------------------------------------------
## whether complaints differ by species
# ------------------------------------------------------------------------------

df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))
load(file.path(here::here(), "analysis/data/derived/spec_icer.rda"))
load(file.path(here::here(), "analysis/data/derived/complete_icer.rda"))

# sort our species by complaints grid
tb <- table(df$complaint,df$species)
colnames(tb) <- tolower(colnames(tb))
colnames(tb)[5] <- "pf/pm/poc/pow"
colnames(tb)[8] <- "pf/poc/pow"
colnames(tb)[14] <- "poc/pow"
huh <- match(names(spec_res_poisson$params$multinom),colnames(tb))
tb <- tb[,huh]
colnames(tb) <- names(spec_res_poisson$params$multinom)
tb[,14] <- 0

# what are the expected ratios
expected <- round(outer(rowSums(tb),l[[2]]$params$multinom ))
tbsp <- as.table(as.numeric(tb))
p <- as.numeric(expected)/sum(expected)
names(p) <- names(tbsp)

# bootstrap them
plot_l = list()
m_list <- names(which(table(df$complaint)>5))[-1]
boots <- list()
for(i in seq_along(m_list)){
  message(i)
  s <- match(m_list[i],rownames(tb))
  boots[[i]] <-
    icer:::coinf_plot(
      reps = 50000,
      probs = expected[s,]/sum(expected) ,
      levels = names(expected[s,]),
      total = sum(expected[s,]),
      plot = TRUE,
      real = tb[s,],
      density = FALSE
    )
  boots[[i]]$plot <- boots[[i]]$plot + ggtitle(m_list[i])
}

# we can see from this that for all the symptoms apart from fever they are
# well explained
cp <- cowplot::plot_grid(plotlist = lapply(boots,"[[","plot"))

# so save the fever plot out
save_figs("supp_fever_composition", boots[[3]]$plot, width = 10, height = 6)

### ---------------------------------
### by month difference
### ---------------------------------

df$month <- month.name[df$month]
tbs <- table(df$month,df$species)
colnames(tbs) <- tolower(colnames(tbs))
colnames(tbs)[5] <- "pf/pm/poc/pow"
colnames(tbs)[8] <- "pf/poc/pow"
colnames(tbs)[14] <- "poc/pow"
huh <- match(names(l[[2]]$params$multinom),colnames(tbs))
tbs <- tbs[,huh]
colnames(tbs) <- names(l[[2]]$params$multinom)
tbs[,14] <- 0

# expected diffs
expected <- round(outer(rowSums(tbs),l[[2]]$params$multinom ))
names(expected) <- month.name
tbsp <- as.table(as.numeric(tbs))
p <- as.numeric(expected)/sum(expected)
names(p) <- names(tbsp)

# bootstrap them
plot_l = list()
boots <- list()
for(i in seq_len(nrow(tbs))){
  message(i)
  s <- rownames(tbs)[i]
  boots[[i]] <-
    icer:::coinf_plot(
      reps = 50000,
      probs = expected[s,]/sum(expected[s,]) ,
      levels = names(expected[s,]),
      total = sum(expected[s,]),
      plot = TRUE,
      real = tbs[s,],
      density = FALSE
    )
  boots[[i]]$plot <- boots[[i]]$plot + ggtitle(s)
}

# we can see here that all species compositions apart from pf/pow are explained well
cp <- cowplot::plot_grid(plotlist = lapply(boots,"[[","plot"))

## there is possibly a monthly difference for pf/pow
for(i in 1:12) {boots[[i]]$vals$month <- rownames(expected)[i]}
x <- do.call(rbind, lapply(boots, function(x) {x$vals}))

# nice density plot using melted data.frame
month_plot <- function(x, tbs,
                       var = "variable", val = "value",
                       quantiles = c(0.025, 0.975),
                       scales = "free", ncol = 6,
                       density = TRUE, spec = "pf/pow",
                       spec_label = "Pf/PoW") {

  levels <- colnames(tbs)
  x <- filter(x, variable == spec)
  real <- data.frame("value" = as.numeric(tbs),
                     "variable" = as.character(sapply(colnames(tbs), rep, 12)),
                     "month" = rep(rownames(tbs), ncol(tbs))) %>%
    filter(variable == spec_label)

  sum_melt <- group_by(x, variable, month) %>%
    summarise(
      q_low = quantile(.data[[val]], quantiles[1]),
      q_high = quantile(.data[[val]], quantiles[2]),
      mdx = median(.data[[val]])
    )

  df_dens <- group_by(x, month) %>% do(dens = density(.[[val]]))
  df_dens <- rbind_list_base(by(df_dens, df_dens$month, function(x) {
    data.frame("x" = x$dens[[1]]$x,
               "y" = x$dens[[1]]$y,
               "month" = x$month[1])
  }))

  for (v in sum_melt$month) {
    df_dens <- df_dens[-which(c(df_dens$month == v &
                                  (df_dens$x > sum_melt$q_high[sum_melt$month == v] |
                                     df_dens$x < sum_melt$q_low[sum_melt$month == v]))), ]
  }

  area_mdx <- apply(sum_melt, 1, function(x) {
    up <- which.min(abs(df_dens$x[df_dens$month == x[2]] - as.numeric(x[4])))
    which(
      near(df_dens$x,df_dens$x[df_dens$month == x[2]][up], tol = 0.001) &
        df_dens$month == x[2]
    )
  })

  real[[var]] <- factor(real[[var]], levels = levels)
  real[["month"]] <- factor(real[["month"]], levels = month.name)


  new_x <- do.call(rbind,lapply(unique(x$month), function(m){
      vls <- table(x$value[x$month==m])/sum(table(x$value[x$month==m]))
      data.frame(x = as.numeric(names(vls)), y = as.numeric(vls), variable = spec, month = m)
    }))

  alt_x <- new_x
    for(j in unique(x$month)) {
      alt_x <- alt_x[!(alt_x$month == j & alt_x$x < sum_melt$q_low[sum_melt$month == j]),]
      alt_x <- alt_x[!(alt_x$month == j & alt_x$x > sum_melt$q_high[sum_melt$month == j]),]
    }

  quant_x <- new_x
    for(j in unique(x$month)) {
      quant_x <- quant_x[(quant_x$month == j& quant_x$x < sum_melt$q_low[sum_melt$month == j]) |
                           (quant_x$month == j& quant_x$x > sum_melt$q_high[sum_melt$month == j]) |
                           quant_x$month != j ,]

    }

  int_breaks <- function(x, n = 5) {
    l <- pretty(x, n)
    l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
  }

  new_x$month <- factor(new_x$month, levels = month.name)
  sum_melt$month <- factor(sum_melt$month, levels = month.name)
  real$month <- factor(real$month, levels = month.name)
  alt_x$month <- factor(alt_x$month, levels = month.name)
  quant_x$month <- factor(quant_x$month, levels = month.name)

  g1 <- ggplot(new_x) +
    geom_bar(aes(x=x, y = y), stat = "identity", fill = "white", color = "black", width = 1, size = 1) +
    geom_bar(aes(x=x, y = y), data = alt_x, stat = "identity", fill = "#002366", linetype = 0, width = 1) +
    geom_bar(aes(x=x, y = y), data = quant_x, stat = "identity", fill = "white", linetype = 0,width = 1) +
    # geom_step(aes(x=.data[[val]], y=..density..), data = x, colour = "black",
    #           stat="bin", binwidth=1, direction = "mid") +
    geom_vline(data = sum_melt, mapping = aes(xintercept = .data$mdx),
               col = "white", size = 1) +
    geom_vline(data = real, mapping = aes(xintercept = .data[[val]]),
               col = "#e62400", size = 1, linetype = "dashed") +
    xlab("Counts") + ylab("Density") +
    scale_x_continuous(breaks = int_breaks) +
    facet_wrap(~month, scales = "free", ncol = 6)  + theme_bw() +
    theme(axis.text.x = element_text(),
          panel.spacing = unit(10, units = "pt"),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = 12))


  g1

}
g1 <- month_plot(x, table(df$month,df$species))
save_figs("temporal_pf_pow", g1, 10,6)

