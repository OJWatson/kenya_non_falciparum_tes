library(tidyverse)

# read in our daily data
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))

# What are the colors and species labels we will be using

colors15 <- c(
  "dodgerblue2",
  "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "grey",
  "#FB9A99", # lt pink
  "palegreen2",
  "deeppink1", "blue1",  "green1", "yellow4",
  "darkorange4", "brown"
)

colors <- c("#90c3d5", "#c35e35", "#8bd757", "#c155c5")
specs <- c("Pf","Pm","Poc","Pow")
spec_names <- names(gen_color_combinations(specs, colors))

# ------------------------------------------------------------------------------
# . Read in our data and format
# ------------------------------------------------------------------------------

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

daily <- daily %>%
  filter(infx != "Missing") %>%
  filter(subjectid %in% names(which(table(subjectid) == 7))) %>%
  filter(infx != "Uninfected") %>%
  group_by(subjectid) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  mutate(cryptic = FALSE) %>%
  mutate(cryptic = replace(cryptic, which(subjectid %in% cryptic_ids), TRUE)) %>%
  mutate(infx = factor(infx)) %>%
  mutate(infx = factor(infx, levels = subset(spec_names, spec_names %in% levels(infx))))


repeat_infection_pf_gg <- daily %>%
  ggplot(aes(x = day, y = fct_reorder(as.character(subjectid), n, .desc = TRUE), color = infx)) +
  geom_point(aes(shape = pf), size = 2) +
  geom_point(aes(shape = pf), size = 3, data = . %>% filter(cryptic), color = "black") +
  geom_point(aes(shape = pf), size = 2) +
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
  ylab("\nIndividuals postive for infection within 42 days after treatment\n") +
  xlab("\nDays after treatment")

repeat_infection_cryptic_pf_gg <- daily %>%
  ggplot(aes(x = day, y = fct_reorder(as.character(subjectid), n, .desc = TRUE), color = infx)) +
  geom_point(aes(shape = pf), size = 3, data = . %>% filter(cryptic), color = "black") +
  geom_point(aes(shape = pf), size = 2, data = . %>% filter(cryptic)) +
  scale_color_manual(values = colors15[seq_along(levels(daily$infx))],
                     labels = levels(daily$infx),
                     drop = FALSE,
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
  ylab("Individuals positive for new species \nwithin 28 days after treatment\n") +
  xlab("\nDays after treatment")


repeat_infection_act_gg <- daily %>%
  ggplot(aes(x = day, y = fct_reorder(as.character(subjectid), n, .desc = TRUE), color = infx)) +
  geom_point(aes(shape = arm), size = 2) +
  geom_point(aes(shape = arm), size = 3, data = . %>% filter(cryptic), color = "black") +
  geom_point(aes(shape = arm), size = 2) +
  scale_color_manual(values = colors15[seq_along(levels(daily$infx))],
                     labels = levels(daily$infx),
                     drop = TRUE,
                     name = "Infection \nComposition") +
  scale_x_continuous(breaks = seq(0, 42, 7)) +
  scale_shape(name = "Drug Treatment") +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(legend.position = "right",
        axis.line = element_line(),
        panel.grid.major.y = element_line(size = 0.25, color = "grey", linetype = "solid"),
        panel.grid.major.x = element_line(size = 0),
        legend.key = element_rect(fill = "white"),
        axis.text.y = element_text(size = 6)) +
  ylab("\nIndividuals postive for infection within 42 days after treatment\n") +
  xlab("\nDays after treatment")

repeat_infection_cryptic_act_gg <- daily %>%
  ggplot(aes(x = day, y = fct_reorder(as.character(subjectid), n, .desc = TRUE), color = infx)) +
  geom_point(aes(shape = arm), size = 3, data = . %>% filter(cryptic), color = "black") +
  geom_point(aes(shape = arm), size = 2, data = . %>% filter(cryptic)) +
  scale_color_manual(values = colors15[seq_along(levels(daily$infx))],
                     labels = levels(daily$infx),
                     drop = FALSE,
                     name = "Infection \nComposition") +
  scale_x_continuous(breaks = seq(0, 42, 7)) +
  scale_shape(name = "Drug Treatment") +
  ggpubr::theme_pubclean(base_size = 14) +
  theme(legend.position = "right",
        axis.line = element_line(),
        panel.grid.major.y = element_line(size = 0.25, color = "grey", linetype = "solid"),
        panel.grid.major.x = element_line(size = 0),
        legend.key = element_rect(fill = "white"),
        axis.text.y = element_text(size = 6)) +
  ylab("Individuals positive for new species \nwithin 28 days after treatment\n") +
  xlab("\nDays after treatment")


repeat_plot <- cowplot::plot_grid(
  repeat_infection_pf_gg, NA, repeat_infection_act_gg,
  repeat_infection_cryptic_pf_gg, NA, repeat_infection_cryptic_act_gg,
  ncol = 3,
  rel_heights = c(1,1/3),
  rel_widths = c(1,0,1),
  labels = c("a)", "", "b)", "c)", "", "d)")
)

save_figs("repeat_infections", repeat_plot, width = 16, height = 20)
save_figs("repeat_infection_act", repeat_infection_act_gg + xlab("Days after treatment"), width = 9, height = 10)

repeat_infection_cryptic_act_gg$layers[[1]]$aes_params$size <- 5
repeat_infection_cryptic_act_gg$layers[[2]]$aes_params$size <- 4
save_figs("repeat_infection_cryptic_act", repeat_infection_cryptic_act_gg + xlab("Days after treatment"), width = 9, height = 6)

