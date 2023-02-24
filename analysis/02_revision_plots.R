library(tidyverse)
devtools::load_all()

binconf_ci <- function(x,i) {
  Hmisc::binconf(sum(x), length(x))[i]
}

# ------------------------------------------------------------------------------
#  ACPR by species figure
# ------------------------------------------------------------------------------

# create table one from our data
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))

dailymod <- daily %>% filter(day == 0) %>%
  filter(day28 != "N/A") %>%
  filter(day42 != "N/A") %>%
  filter(infx != "Missing") %>%
  mutate(day28acpr = day28 == "ACPR") %>%
  mutate(day42acpr = day42 == "ACPR")

glm(day28acpr ~ pf, data = dailymod, family = "binomial") %>% summary
glm(day28acpr ~ pm, data = dailymod, family = "binomial") %>% summary
glm(day28acpr ~ poc, data = dailymod, family = "binomial") %>% summary
glm(day28acpr ~ pow, data = dailymod, family = "binomial") %>% summary
glm(day28acpr ~ pf+pm+poc+pow, data = dailymod, family = "binomial") %>% summary

acpr28_summs <- rbind(
  dailymod %>%
    filter(pf) %>%
    summarise(infx = "P. falciparum",
              n = n(),
              med = binconf_ci(day28acpr,1),
              min = binconf_ci(day28acpr,2),
              max = binconf_ci(day28acpr,3)),
  dailymod %>%
    filter(pm) %>%
    summarise(infx = "P. malariae",
              n = n(),
              med = binconf_ci(day28acpr,1),
              min = binconf_ci(day28acpr,2),
              max = binconf_ci(day28acpr,3)),
  dailymod %>%
    filter(poc) %>%
    summarise(infx = "P. ovale curtisi",
              n = n(),
              med = binconf_ci(day28acpr,1),
              min = binconf_ci(day28acpr,2),
              max = binconf_ci(day28acpr,3)),
  dailymod %>%
    filter(pow) %>%
    summarise(infx = "P. ovale wallikeri",
              n = n(),
              med = binconf_ci(day28acpr,1),
              min = binconf_ci(day28acpr,2),
              max = binconf_ci(day28acpr,3))
)



acpr42_summs <- rbind(
  dailymod %>%
    filter(pf) %>%
    summarise(infx = "P. falciparum",
              n = n(),
              med = binconf_ci(day42acpr,1),
              min = binconf_ci(day42acpr,2),
              max = binconf_ci(day42acpr,3)),
  dailymod %>%
    filter(pm) %>%
    summarise(infx = "P. malariae",
              n = n(),
              med = binconf_ci(day42acpr,1),
              min = binconf_ci(day42acpr,2),
              max = binconf_ci(day42acpr,3)),
  dailymod %>%
    filter(poc) %>%
    summarise(infx = "P. ovale curtisi",
              n = n(),
              med = binconf_ci(day28acpr,1),
              min = binconf_ci(day28acpr,2),
              max = binconf_ci(day28acpr,3)),
  dailymod %>%
    filter(pow) %>%
    summarise(infx = "P. ovale wallikeri",
              n = n(),
              med = binconf_ci(day42acpr,1),
              min = binconf_ci(day42acpr,2),
              max = binconf_ci(day42acpr,3))
)


acpr28 <- acpr28_summs %>%
  ggplot(aes(x = infx, y = med, ymin = min, ymax = max)) +
  geom_pointrange() +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        axis.text.x = element_text(face = "italic")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Day 28 ACPR") +
  xlab("Detected at Day 0")

acpr42 <- acpr42_summs %>%
  ggplot(aes(x = infx, y = med, ymin = min, ymax = max)) +
  geom_pointrange() +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        axis.text.x = element_text(face = "italic")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Day 42 ACPR") +
  xlab("Detected at Day 0")

spec_acpr <- cowplot::plot_grid(acpr28)
save_figs("species_acpr28", spec_acpr, height = 4, width = 6)

# ------------------------------------------------------------------------------
#  Microscopy vs PCR Clearance
# ------------------------------------------------------------------------------

# read in our data
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds")) %>%
  mutate(pf_only = infx == "Pf")

hourly <- readRDS(cp_path("analysis/data/derived/hourly_data.rds"))

# join this with our clearance data
combdf <- left_join(
  hourly,
  daily %>% select(subjectid, infx, pf_only, day),
  by = c("subjectid", "day")
)

# select just the time points that have both micro and PCR
combdf <- combdf %>% group_by(subjectid) %>% mutate(para = gmeanpfemia/max(gmeanpfemia)) %>%
  filter(!(time %in% c(4,8,12,16,18,20,24,30,36,42,48,54,60,66))) %>%
  filter(!is.na(infx))

# filter out missing PCR data and select only samples with all time points recorded
combdf <- combdf %>%
  group_by(day) %>%
  filter(infx != "Missing") %>%
  filter(subjectid %in% as.integer(names(which(table(combdf$subjectid)==7))))

# create inf plot
inf_prop <- combdf %>%
  summarise(n = n(),
            micro = sum(gmeanpfemia == 0)/n,
            pcr = sum(infx == "Uninfected", na.rm = TRUE)/n) %>%
  pivot_longer(micro:pcr) %>%
  ggplot(aes(day, 1-value, color = name)) +
  geom_step(size = 2) +
  theme_bw() +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(axis.line = element_line(),
        panel.grid.major.y = element_line(size = 0.25, color = "grey", linetype = "solid"),
        panel.grid.major.x = element_line(size = 0),
        axis.text.y = element_text(size = 6),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(color = "black", size = 1),
        legend.position = c(0.8, 0.84)) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_viridis_d(name = "Diagnostic", labels = c("Microscopy", "PCR"), end = 0.8) +
  ylab("% Infected by Diagnostic Method") +
  xlab("\nTime (days)")

pcr_non_clears <- combdf %>%
  filter(day %in% c(7,14) & infx != "Uninfected") %>%
  pull(subjectid) %>% unique()

daily <- daily %>%
  filter(infx != "Missing") %>%
  #filter(subjectid %in% names(which(table(subjectid) == 7))) %>%
  filter(infx != "Uninfected") %>%
  group_by(subjectid) %>%
  mutate(n = n())

non_clear_pcrs <- daily %>%
  ggplot(aes(x = day, y = fct_reorder(as.character(subjectid), n, .desc = TRUE), color = infx)) +
  geom_point(aes(shape = pf), size = 3, data = . %>% filter(subjectid %in% pcr_non_clears), color = "black") +
  geom_point(aes(shape = pf), size = 2, data = . %>% filter(subjectid %in% pcr_non_clears)) +
  scale_color_manual(values = colors15[seq_along(levels(daily$infx))],
                     labels = levels(daily$infx),
                     drop = TRUE,
                     name = "Infection \nComposition") +
  scale_x_continuous(breaks = seq(0, 42, 7)) +
  scale_shape(name = "P. falciparum \nInfection", labels = c("Pf-", "Pf+")) +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(legend.position = "right",
        axis.line = element_line(),
        panel.grid.major.y = element_line(size = 0.25, color = "grey", linetype = "solid"),
        panel.grid.major.x = element_line(size = 0),
        legend.key = element_rect(fill = "white"),
        axis.text.y = element_text(size = 6)) +
  ylab("\nIndividuals positive for infection by PCR \non days 7 or 14 after treatment\n") +
  xlab("\nDays after treatment")

micro_pcr_clear_fails <- cowplot::plot_grid(inf_prop, non_clear_pcrs, ncol = 2, align = "h", rel_widths = c(0.75,1), labels = "auto")
save_figs("micro_vs_pcr_clears", micro_pcr_clear_fails, height = 6, width = 12)


# ------------------------------------------------------------------------------
#  Reinfection work
# ------------------------------------------------------------------------------

# time to reinfection
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
slide <- readRDS(cp_path("analysis/data/derived/slide_data.rds"))

identify_cryptic <- function(pf, pm, poc, pow, cut_off = 4) {

  day_0 <- which(c(pf[1], pm[1], poc[1], pow[1]))
  not_day_0 <- setdiff(1:4, day_0)

  day7_onwards <- cbind(pf[-1], pm[-1], poc[-1], pow[-1])
  sum(day7_onwards[seq_len(cut_off),not_day_0])>0

}

cryptic_ids <- daily %>% group_by(subjectid) %>%
  summarise(cryptic = identify_cryptic(pf,pm,poc,pow)) %>%
  filter(cryptic) %>%
  pull(subjectid)


reinfection_data <- daily %>%
  filter(!subjectid %in% cryptic_ids) %>%
  group_by(subjectid) %>%
  mutate(all_inf = (sum(infx == "Uninfected") == n())) %>%
  filter(infx != "Missing") %>%
  filter(!all_inf) %>%
  filter(subjectid %in% names(which(table(subjectid) == 7))) %>%
  group_by(subjectid, arm) %>%
  summarise(dayuninf = day[which(infx == "Uninfected")[1]],
            dayreinf = day[day>dayuninf][which(infx[day>dayuninf] != "Uninfected")[1]]) %>%
  filter(!is.na(dayuninf)) %>%
  mutate(dayreinf = replace_na(dayreinf, 49)) %>%
  mutate(proph_duration = dayreinf - dayuninf) %>%
  na.omit()



pos_test <- function(x){any(as.character(x) == "POS", na.rm = TRUE)}
retreated <- which(apply(slide[slide$subjectid %in% unique(reinfection_data$subjectid),3:8], 1, pos_test))
retreated <- slide$subjectid[which(slide$subjectid %in% unique(reinfection_data$subjectid))][retreated]

# remove those that had been retreated
reinfection_data_surv <- reinfection_data %>% filter(!(subjectid %in% retreated))
reinfection_data_surv <- reinfection_data_surv %>%
  mutate(status = 1) %>%
  mutate(status = replace(status, dayreinf == 49, 0))

surv_gg <- ggsurvfit::survfit2(survival::Surv(proph_duration, status) ~ arm, data = reinfection_data_surv) %>%
  ggsurvfit::ggsurvfit(type = "survival", size = 1.5) +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) +
  ggsurvfit::add_confidence_interval( size = 1.5) +
  ggsurvfit::add_confidence_interval(type = "lines", size = 1) +
  geom_point(color = "black", size = 3) +
  geom_point(size = 2) +
  #geom_errorbar(width =1, size = 1) +
  scale_color_viridis_d(name = "ACT:", end = 0.8) +
  scale_fill_viridis_d(name = "ACT:", end = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  ylab("% Uninfected") +
  xlab("Days Since First Clearing All Infections") +
  theme(legend.position = "top", legend.background = element_rect(color = "white"))


surv_table <- survival::coxph(
  survival::Surv(proph_duration, status) ~ ACT,
  data = reinfection_data_surv %>% mutate(ACT = arm, ACT = replace(ACT, ACT == "AL", "AL (Ref)"))) %>%
  gtsummary::tbl_regression(exp = TRUE)

surv_fig <- cowplot::plot_grid(
  surv_gg, NA, as_ggplot(surv_table),
  ncol = 4, rel_widths = c(1,0.05,0.5,0.025),
  labels = c("a", "", "b", "")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))
save_figs("surv_fig", surv_fig, height = 6, width = 10)

# ------------------------------------------------------------------------------
#  hgb work
# ------------------------------------------------------------------------------

# time to reinfection
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
slide <- readRDS(cp_path("analysis/data/derived/slide_data.rds"))

identify_cryptic <- function(pf, pm, poc, pow, cut_off = 4) {

  day_0 <- which(c(pf[1], pm[1], poc[1], pow[1]))
  not_day_0 <- setdiff(1:4, day_0)

  day7_onwards <- cbind(pf[-1], pm[-1], poc[-1], pow[-1])
  sum(day7_onwards[seq_len(cut_off),not_day_0])>0

}

cryptic_ids <- daily %>% group_by(subjectid) %>%
  summarise(cryptic = identify_cryptic(pf,pm,poc,pow)) %>%
  filter(cryptic) %>%
  pull(subjectid)

# identify which were pf only at day 0
daily <- daily %>% filter(day == 0 & !(infx %in% c("Uninfected", "Missing"))) %>%
  mutate(pf_only = infx == "Pf")

daily %>% ggplot(aes(gmeanpfemia, hgb, color = arm)) +
  geom_point() +
  geom_smooth(method = "lm")

hgb_a <- ggplot(daily, aes(hgb, slope_half_life, color = arm, fill = arm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  xlab("Haemoglobin Levels on Enrollment") +
  ylab("Parasite Clearance Slope Half Life") +
  scale_color_viridis_d(name = "ACT", end = 1) +
  scale_fill_viridis_d(name = "ACT", end = 1)


hgb_b <- ggplot(daily, aes(hgb, log(gmeanpfemia), color = arm, fill = arm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  xlab("Haemoglobin Levels on Enrollment") +
  ylab("Log Parasitaemia") +
  scale_color_viridis_d(name = "ACT", end = 1) +
  scale_fill_viridis_d(name = "ACT", end = 1)

hgb_plot <- cowplot::plot_grid(
  hgb_a + theme(legend.position = "none"),
  NA + theme(plot.background = element_rect(fill = "white", color = "white")),
  hgb_b + theme(legend.position = "none"),
  cowplot::get_legend(hgb_a),
  ncol = 4, rel_widths = c(1,0.05,1,0.25),
  labels = c("a", "", "b", "")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

save_figs("hgb_effect_plot", hgb_plot, height = 4, width = 10)
