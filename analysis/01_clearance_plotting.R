library(tidyverse)

### -------------------------------------
### Parasite Clearance Curves
### -------------------------------------

# Fit models to explain parasite clearance by drug, pf and individual covariates

# read in our data
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
hourly <- readRDS(cp_path("analysis/data/derived/hourly_data.rds"))

# identify which were pf only at day 0
daily <- daily %>% filter(day == 0 & !(infx %in% c("Uninfected", "Missing"))) %>%
  mutate(pf_only = infx == "Pf")

# join this with our clearance data
arm_df <- left_join(hourly, daily %>% select(subjectid, arm, pf_only, age, gender, weight, slope_half_life, pc50, pc99, height)) %>%
  filter(time < 43 & gmeanpfemia < 1e7) %>%
  filter(!is.na(arm))
arm_df2 <- arm_df %>% group_by(subjectid) %>% mutate(para = gmeanpfemia/max(gmeanpfemia))

# create plot of clearance curvs with logistic decay
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),...)
}

# Incorrect transformation in pseudo log scale
# clearance_gg <- arm_df2 %>%
#   ggplot(aes(time, para, group = subjectid, color = arm)) +
#   geom_line(alpha = 0.1) +
#   binomial_smooth(aes(time, para, color = arm, fill = arm),
#                   inherit.aes = FALSE, lwd = 1, se = TRUE, span = 0.2, alpha = 0.3, level = 0.95, fullrange = TRUE) +
#   binomial_smooth(aes(time, para, color = arm, fill = arm),
#                   inherit.aes = FALSE, lwd = 1, se = FALSE, span = 0.9, alpha = 0.3, level = 0.95, fullrange = TRUE) +
#   scale_y_continuous(labels = scales::percent,
#                      #limits = c(0,1),
#                      trans=scales::pseudo_log_trans(sigma = 0.0001, base = exp(10)),
#                      breaks = c(0, 0.001,0.01,0.1,1)) +
#   ylab("Parasitaemia Relative to Max Parasitaemia (%)") +
#   xlab("Hours since treatment") +
#   scale_color_viridis_d(name = "ACT", end = 1) +
#   scale_fill_viridis_d(name = "ACT", end = 1) +
#   theme_bw()

clearance_gg <- arm_df2 %>%
  mutate(para = replace(para, para < 1e-4, 1e-4)) %>%
  ggplot(aes(time, para, p = para, color = arm)) +
  geom_line(alpha = 0.1, aes(group = subjectid)) +
  geom_smooth(method="glm",se=FALSE,  method.args=list(family="binomial"),
              formula = cbind(p) ~ x,
              data = arm_df2 %>%
                mutate(para = replace(para, para < 1e-4, 1e-4)) %>%
                filter((arm == "AL" & time > 8) | (arm != "AL" & time > 4)),
              aes(y = stage(para,y), color = arm), fullrange = TRUE) +
  #geom_smooth(se = FALSE)+
  scale_y_log10(labels = c("*0.01%", "0.1%","1%","10%","100%"),
                limits = c(1e-4,1),
                breaks = c(0.0001, 0.001,0.01,0.1,1)) +
  ylab("Parasitaemia Relative to Max Parasitaemia (%)") +
  xlab("Hours since treatment") +
  scale_color_viridis_d(name = "ACT", end = 1) +
  scale_fill_viridis_d(name = "ACT", end = 1) +
  theme_bw()

clearance_gg2 <- arm_df2 %>%
  mutate(para = replace(para, para < 1e-4, 1e-4)) %>%
  ggplot(aes(time, para, p = para, pf = pf_only, linetype = pf_only, color = arm)) +
  geom_line(aes(group = subjectid), alpha = 0.1) +
  geom_smooth(method="glm",se=FALSE,  method.args=list(family="binomial"),
              formula = cbind(p) ~ x,
              aes(y = stage(para,y), color = arm, linetype = pf_only), fullrange = TRUE,
              data = arm_df2 %>%
                mutate(para = replace(para, para < 1e-4, 1e-4)) %>%
                filter((arm == "AL" & time > 8) | (arm != "AL" & time > 4))) +
  scale_y_log10(labels = c("*0.01%", "0.1%","1%","10%","100%"),
                limits = c(1e-4,1),
                breaks = c(0.0001, 0.001,0.01,0.1,1)) +
  ylab("Parasitaemia Relative to Max Parasitaemia (%)") +
  xlab("Hours since treatment") +
  scale_linetype(name = "Single Pf Infection \nat Hour 0") +
  guides(linetype = guide_legend(override.aes = list(color = "black" ) ) ) +
  scale_color_viridis_d(name = "ACT", end = 1) +
  scale_fill_viridis_d(name = "ACT", end = 1) +
  theme_bw()

# Incorrect transformation in pseudo log scale
# clearance_gg2 <- arm_df2 %>%
#   ggplot(aes(time, para, group = subjectid, color = arm, linetype = pf_only)) +
#   geom_line(alpha = 0.1) +
#   binomial_smooth(aes(time, para, color = arm, fill = arm, linetype = pf_only),
#                   inherit.aes = FALSE, lwd = 1, se = FALSE, span = 0.9, alpha = 0.3, level = 0.95) +
#   scale_y_continuous(labels = scales::percent,
#                      #limits = c(0,1),
#                      trans=scales::pseudo_log_trans(sigma = 0.0001, base = exp(10)),
#                      breaks = c(0, 0.001,0.01,0.1,1)) +
#   ylab("Parasitaemia Relative to Max Parasitaemia (%)") +
#   xlab("Hours since treatment") +
#   scale_color_viridis_d(name = "ACT", end = 1) +
#   scale_fill_viridis_d(name = "ACT", end = 1) +
#   scale_linetype(name = "Single Pf Infection \nat Hour 0") +
#   guides(linetype = guide_legend(override.aes = list(color = "black" ) ) ) +
#   theme_bw()

clearance_plot <- cowplot::plot_grid(clearance_gg, clearance_gg2, labels = "auto")
save_figs("clearance_plot", clearance_plot, width = 12,height = 4.5)

### -------------------------------------
### Fit models to explain parasite clearance by drug, pf and individual covariates
### -------------------------------------

slope_glm <- glm(slope_half_life  ~  hgb + age + weight + gender + arm + pf_only + gmeanpfemia, data = daily)
pc50_glm <- glm(pc50  ~  hgb + age + weight + gender + arm + pf_only + gmeanpfemia, data = daily)
pc99_glm <- glm(pc99 ~  hgb + age + weight + gender + arm + pf_only + gmeanpfemia, data = daily)

covariate_labels <- c("Haemoglobin Levels g/dl", "Age (years)", "Weight (kg)",
                      "Gender = Male", "ACT = ASMQ", "ACT = DHA-PPQ", "Pf single infection",
                      "Parasitaemia (p/uL)")

writeLines(
  stargazer::stargazer(list(slope_glm, pc50_glm, pc99_glm),
                       report=('vc*p'),
                       keep.stat = c("n","ll"),
                       column.labels = c("PC Slope Half Life",
                                         "PC 50% Clearance",
                                         "PC 99% Clearance",
                                         "ACPR Day 28",
                                         "ACPR Day 42"),
                       column.sep.width = "90pt",
                       digits = 4,
                       dep.var.labels.include = FALSE,
                       type="html",
                       covariate.labels = covariate_labels,
                       title = "Generalised Linear Modelling of Parasite Clearance and Late Parasitological Failure",
                       omit = "Constant", dep.var.caption = ""),
  cp_path("analysis/tables/parasite_clearance_failure.html")
)

### -------------------------------------
### Effect sizes of covariates on slope half life
### -------------------------------------

# Tidy up mixed model data to be ORs from binom model
td <- broom.mixed::tidy(slope_glm, conf.int=TRUE)

# label our covariates more clearly
td$term[-1] <- covariate.labels

# format the labels and title
td$label <- paste0(round(td$estimate,2), " [",
                   round(td$conf.low,digits = 2), ", ",
                   round(td$conf.high,digits = 2), "]")
td <- td[-1,]
td <- rbind(td, td[1,])
td[nrow(td),2:ncol(td)] <- NA
td$term[nrow(td)] <- ""
td$label[nrow(td)] <- "Effect Size [95% CI]"

# order covariates sensible order
td <- td[c(2,3,4,1,5:7, 8),]
td$term <- factor(td$term, levels = td$term)

# create effect size plot
gg_effs <- ggplot(td, aes(estimate, term, color = p.value<0.05)) +
  geom_vline(xintercept = 0) + theme_bw() +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height=0.25)) +
  scale_color_manual(name="p-value < 0.05",values = c("#002366","#e62400")) +
  theme(text = element_text(size = 14), panel.spacing = unit(10,units = "pt"),) +
  xlab("Effect Size on Parasite Clearance Half Life") + ylab("") +
  geom_text(inherit.aes = FALSE, data = td,
            mapping = aes(y = term, x = 1, label = label)) +
  theme(legend.position = "none",
        panel.grid.major.y =  element_blank(),
        panel.grid.minor.y =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.4)) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-0.75,1.51))

# Add effect size desc to right
gg_effs <- gg_effs + geom_curve(data = data.frame(x = 0.5, y = nrow(td)-0.25,
                                  xend = 1.5, yend = nrow(td)-0.25),
                mapping = aes(x = x,y = y, xend = xend, yend = yend),
                angle = 0L, arrow = arrow(30L, unit(0L, "inches"), "last", "closed"),
                alpha = 1, inherit.aes = FALSE)

save_figs("slope_effect_sizes", gg_effs, height = 4, width = 6)


hgb_demo <- ggplot(daily, aes(hgb, slope_half_life, color = arm, fill = arm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  xlab("Haemoglobin Levels on Enrollment") +
  ylab("Parasite Clearance Slope Half Life") +
  scale_color_viridis_d(name = "ACT", end = 1) +
  scale_fill_viridis_d(name = "ACT", end = 1)
save_figs("hgb_effect_plot", hgb_demo, height = 4, width = 6)


### -------------------------------------
### Time to reinfection plots
### -------------------------------------

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

acpr28_glm <- glm(day28=="ACPR" ~  hgb + age + weight + gender + arm + pf_only + gmeanpfemia,
                  data = daily %>% filter(!subjectid %in% cryptic_ids), family = "binomial")
acpr42_glm <- glm(day42=="ACPR" ~  hgb + age + weight + gender + arm + pf_only + gmeanpfemia,
                  data = daily %>% filter(!subjectid %in% cryptic_ids), family = "binomial")

covariate_labels <- c("Haemoglobin Levels g/dl", "Age (years)", "Weight (kg)",
                      "Gender = Male", "ACT = ASMQ", "ACT = DHA-PPQ", "Pf single infection",
                      "Parasitaemia (p/uL)")

writeLines(
  stargazer::stargazer(list(acpr28_glm, acpr42_glm),
                       report=('vc*p'),
                       keep.stat = c("n","ll"),
                       column.labels = c("ACPR Day 28",
                                         "ACPR Day 42"),
                       column.sep.width = "90pt",
                       digits = 4,
                       dep.var.labels.include = FALSE,
                       type="html",
                       covariate.labels = covariate_labels,
                       title = "Generalised Linear Modelling of Parasite Clearance and Late Parasitological Failure",
                       omit = "Constant", dep.var.caption = ""),
  cp_path("analysis/tables/ACPR.html")
)


# focussing on those that aren't cryptic infections as we know these are not recrudescences
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
reinfection_data <- daily %>%
  filter(!subjectid %in% cryptic_ids) %>%
  group_by(subjectid) %>%
  mutate(all_inf = (sum(infx == "Uninfected") == n())) %>%
  filter(!all_inf) %>%
  filter(subjectid %in% names(which(table(subjectid) == 7))) %>%
  group_by(subjectid, arm) %>%
  summarise(dayuninf = day[which(infx == "Uninfected")[1]],
            dayreinf = day[day>dayuninf][which(infx[day>dayuninf] != "Uninfected")[1]]) %>%
  filter(!is.na(dayuninf)) %>%
  mutate(dayreinf = replace_na(dayreinf, 56)) %>%
  mutate(proph_duration = dayreinf - dayuninf) %>%
  na.omit()

reinfection_data <- reinfection_data %>%
  group_by(subjectid, arm) %>%
  summarise(reinfected = proph_duration) %>%
  mutate(day = reinfected) %>%
  mutate(reinfected = TRUE) %>%
  group_by(subjectid, arm) %>%
  complete(day = seq(0, 42, 7)) %>%
  group_by(subjectid, arm) %>%
  fill(reinfected, .direction = "down") %>%
  mutate(reinfected = replace_na(reinfected, FALSE)) %>%
  mutate(weight = 1) %>%
  mutate(weight = replace(weight, which(day == 0), 1))

pos_test <- function(x){any(as.character(x) == "POS", na.rm = TRUE)}

retreated <- which(apply(slide[slide$subjectid %in% unique(reinfection_data$subjectid),3:8], 1, pos_test))
retreated <- slide$subjectid[which(slide$subjectid %in% unique(reinfection_data$subjectid))][retreated]

# remove those that had been retreated
reinfection_data <- reinfection_data %>% filter(!(subjectid %in% retreated))

# create plot of clearance curvs with logistic decay
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"),...)
}

reinf_summary <- reinfection_data %>% filter(day <= 42) %>%
group_by(day, arm) %>%
  summarise(y = Hmisc::binconf(sum(reinfected), n())[1],
            ymin = Hmisc::binconf(sum(reinfected), n())[2],
            ymax = Hmisc::binconf(sum(reinfected), n())[3])

reinfection_gg <- ggplot(reinfection_data %>% filter(day <= 42),
       aes(day, as.integer(reinfected), color = arm, weight = weight, fill = arm)) +
  #geom_point() +
  binomial_smooth(se=TRUE) +
  binomial_smooth(se=FALSE) +
  theme_bw() +
  xlab("Days Since First Clearing All Infections") +
  ylab("% Reinfected")+
  scale_color_viridis_d(name = "ACT", end = 1) +
  scale_fill_viridis_d(name = "ACT", end = 1) +
  geom_errorbar(aes(day, y, ymin = ymin, ymax = ymax, color = arm), reinf_summary, inherit.aes = FALSE, width = 1, size = 1, alpha = 0.75) +
  geom_point(aes(day, y, ymin = ymin, ymax = ymax), color = "black", reinf_summary, inherit.aes = FALSE, size = 3) +
  geom_point(aes(day, y, ymin = ymin, ymax = ymax, color = arm), reinf_summary, inherit.aes = FALSE, size = 2) +
  scale_y_continuous(labels = scales::percent)

save_figs("reinfection", reinfection_gg, height = 4, width = 6)
