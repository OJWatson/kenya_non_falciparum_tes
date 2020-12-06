library(tidyverse)
our_theme <- ggplot2::theme(plot.title = element_text(hjust = 0.5),
                            panel.border = element_blank(),
                            axis.line = element_line(colour = "black"))

# ------------------------------------------------------------------------------
# Fig 3
# ------------------------------------------------------------------------------

# read in raw data
df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))

# Species plots over time
# Summarize values over sites and months and  years
df_summary <- df %>%
  group_by(.dots = c("site", "year", "month")) %>%
  summarize(endemicity_sites = unique(endemicity_sites),
            Pf_mean = mean(pf, na.rm=TRUE),
            Pm_mean = mean(pm, na.rm=TRUE),
            PoW_mean = mean(pow, na.rm=TRUE),
            PoC_mean = mean(poc, na.rm=TRUE),
            Pf_alone = mean(pf & !pm & !pow & !poc, na.rm=TRUE),
            No_Pf = mean(!pf, na.rm=TRUE),
            Pf = sum(pf, na.rm=TRUE),
            Pm = sum(pm, na.rm=TRUE),
            PoW = sum(pow, na.rm=TRUE),
            PoC = sum(poc, na.rm=TRUE),
            n = n()) %>%
  as.data.frame()

# add some useful columns
df_summary$month_cont <- df_summary$year*12 + df_summary$month
df_summary$date <- df_summary$year
df_summary$date[df_summary$date<10] <- paste0("0", df_summary$date[df_summary$date<10])
df_summary$date <- lubridate::ymd(paste0("20",df_summary$date,"-",df_summary$month,"-01"))

melt <- reshape2::melt(df_summary,measure.vars=c("Pm_mean","Pf_mean","PoC_mean","PoW_mean"))
levels(melt$variable) <- c("P. malariae","P. falciparum", "P. ovale curtisi", "P. ovale wallikeri")

species_over_time <- ggplot(melt) +
  geom_point(aes(x=date, y=value, size = n),col="#002366", alpha = 0.6) +
  theme_bw()  + facet_grid(variable~site) +
  theme(axis.text.x = element_text(), panel.spacing = unit(10,units = "pt"),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size=12)) +
  scale_x_date(labels = scales::date_format("'%y")) + scale_y_continuous(labels = scales::percent) +
  xlab("Year") + ylab("Percentage Species") + scale_size(name="Samples")

species_over_time_no_site <- ggplot(melt) +
  geom_point(aes(x=date, y=value, size = n),col="#002366", alpha = 0.6) +
  theme_bw()  + facet_wrap(~variable) +
  theme(axis.text.x = element_text(), panel.spacing = unit(10,units = "pt"),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size=12)) +
  scale_x_date(labels = scales::date_format("'%y")) + scale_y_continuous(labels = scales::percent) +
  xlab("Year") + ylab("Percentage Species") + scale_size(name="Samples")

# Mixed-effects models
df_summary$year_time <- (df_summary$month_cont - 96)/12

curt_mod2 <- lme4::lmer(PoC_mean ~ year_time + (0+year_time|endemicity_sites),
                        data = df_summary, weights = n,
                        REML = TRUE)

wali_mod2 <- lme4::lmer(PoW_mean ~ year_time + (0+year_time|endemicity_sites),
                        data = df_summary, weights = n,
                        REML = TRUE)

malariae_mod2 <- lme4::lmer(Pm_mean ~ year_time + (0+year_time|endemicity_sites),
                            data = df_summary, weights = n,
                            REML = TRUE)

falc_only_mod2 <- lme4::lmer(Pf_alone ~ year_time + (0+year_time|endemicity_sites),
                             data = df_summary, weights = n,
                             REML = TRUE)

# model predictions add to model
df_summary$predictions_malariae <- predict(malariae_mod2,type="response")
df_summary$predictions_falc_only_mod <- predict(falc_only_mod2,type="response")
df_summary$predictions_wali <- predict(wali_mod2,type="response")
df_summary$predictions_curt <- predict(curt_mod2,type="response")

# plots over time with predictions
wali_over_time <- ggplot(df_summary) +
  geom_point(aes(x=date, y=PoW_mean, size = n),col="#002366", alpha = 0.6) + theme_bw() +
  geom_line(aes(x=date, y=predictions_wali),col="#e62400",size=1) +
  xlab("Year") + ylab("% P. ovale wallikeri") + scale_size(name="Samples") +
  facet_wrap(~endemicity_sites,nrow=1) +
  theme(axis.text.x = element_text(), panel.spacing = unit(10,units = "pt"),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size=12)) +
  scale_x_date(labels = scales::date_format("'%y")) + scale_y_continuous(labels = scales::percent)

curt_over_time <- ggplot(df_summary) +
  geom_point(aes(x=date, y=PoC_mean, size = n),col="#002366", alpha = 0.6) + theme_bw() +
  geom_line(aes(x=date, y=predictions_curt),col="#e62400",size=1) +
  xlab("Year") + ylab("% P. ovale curtisi") + scale_size(name="Samples") +
  facet_wrap(~endemicity_sites,nrow=1) +
  theme(axis.text.x = element_text(), panel.spacing = unit(10,units = "pt"),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size=12)) +
  scale_x_date(labels = scales::date_format("'%y")) + scale_y_continuous(labels = scales::percent)

malariae_over_time <- ggplot(df_summary) +
  geom_point(aes(x=date, y=Pm_mean, size = n),col="#002366", alpha = 0.6) + theme_bw() +
  geom_line(aes(x=date, y=predictions_malariae),col="#e62400",size=1) +
  xlab("Year") + ylab("% P. malariae") + scale_size(name="Samples") +
  facet_wrap(~endemicity_sites,nrow=1) +
  theme(axis.text.x = element_text(), panel.spacing = unit(10,units = "pt"),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size=12)) +
  scale_x_date(labels = scales::date_format("'%y")) + scale_y_continuous(labels = scales::percent)

falciparum_only_over_time <- ggplot(df_summary) +
  geom_point(aes(x=date, y=Pf_alone, size = n),col="#002366", alpha = 0.6) + theme_bw() +
  geom_line(aes(x=date, y=predictions_falc_only_mod),col="#e62400",size=1) +
  xlab("Year") + ylab("% P. falciparum Only") + scale_size(name="Samples") +
  facet_wrap(~endemicity_sites,nrow=1) +
  theme(axis.text.x = element_text(), panel.spacing = unit(10,units = "pt"),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size=12)) +
  scale_x_date(labels = scales::date_format("'%y")) + scale_y_continuous(labels = scales::percent)

# combine them
combined_model_plots <- cowplot::plot_grid(
  curt_over_time,
  wali_over_time,
  malariae_over_time,
  falciparum_only_over_time,
  nrow=4
)

save_figs("fig3", combined_model_plots, width = 9.5, height = 10.8)

# ------------------------------------------------------------------------------
# Linear Model Coefficient Table
# ------------------------------------------------------------------------------

## Stargazer
writeLines(
  stargazer::stargazer(curt_mod2, wali_mod2, malariae_mod2, falc_only_mod2,
                       report=('vc*p'),
                       keep.stat = c("n","ll"),
                       column.labels = c("P. ovale curtisi",
                                         "P. ovale wallikeri",
                                         "P. malariae",
                                         "P. falciparum single infections"),
                       column.sep.width = "30pt",
                       digits = 4,
                       dep.var.labels.include = FALSE,
                       type="text",
                       covariate.labels = c("Time (years)"),
                       title = "Parameters for the hierchical analysis of Plasmodium species frequency changes over time",
                       omit = "Constant", dep.var.caption = ""),
  file.path(here::here(), "analysis/tables/linear_model_table.txt")
)

# ------------------------------------------------------------------------------
# Random effects plot
# ------------------------------------------------------------------------------

## Random effects plot
wali_re <- sjPlot::plot_model(wali_mod2, type="re",
                              vline.color="#A9A9A9", dot.size=1.5,
                              show.values=T, value.offset=.2, digits = 3) +
  theme_bw() + ylim(c(-0.05,0.05)) + ggtitle("P. ovale wallikeri")

curt_re <- sjPlot::plot_model(curt_mod2, type="re",
                              vline.color="#A9A9A9", dot.size=1.5,
                              show.values=T, value.offset=.2, digits = 3) +
  theme_bw() + ylim(c(-0.05,0.05)) + ggtitle("P. ovale curtisi")

malariae_re <- sjPlot::plot_model(malariae_mod2, type="re",
                                  vline.color="#A9A9A9", dot.size=1.5,
                                  show.values=T, value.offset=.2, digits = 3) +
  theme_bw() + ylim(c(-0.0125,0.0125)) + ggtitle("P. malariae")

falc_only_re <- sjPlot::plot_model(falc_only_mod2, type="re",
                                   vline.color="#A9A9A9", dot.size=1.5,
                                   show.values=T, value.offset=.2, digits = 3) +
  theme_bw() + ylim(c(-0.05,0.05)) + ggtitle("P. falciparum only")

re_all <- cowplot::plot_grid(curt_re, wali_re, malariae_re, falc_only_re, ncol = 2)
save_figs("random_effects", re_all, width = 7, height = 7)

