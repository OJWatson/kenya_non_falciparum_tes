devtools::load_all()
library(tidyverse)

# create table one from our data
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
daily <- daily %>% filter(day == 0) %>%
  filter(!infx %in% c("Uninfected", "Missing"))

countCharOccurrences <- function( s,char) {
  s2 <- gsub(char,"",s)
  return (nchar(as.character(s)) - nchar(as.character(s2)))
}

slash_to_text <- function(x) {
  x <- gsub("/",", ",x)
  gsub(",([^,]*)$"," and\\1", x)
}

# write these results to table one
fig1 <- daily %>%
  group_by(infx) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(tot = sum(n)) %>%
  group_by(infx) %>%
  mutate(ymin = Hmisc::binconf(n,tot)[2]) %>%
  mutate(ymax = Hmisc::binconf(n,tot)[3]) %>%
  mutate(type = as.factor(sapply(infx, countCharOccurrences, "/") + 1)) %>%
  mutate(type = c("Single","Double","Triple","Quadruple")[type]) %>%
  mutate(type = factor(type ,levels=c("Single","Double","Triple","Quadruple"))) %>%
  mutate(infx = slash_to_text(infx)) %>%
  ggplot(aes(x = prop, y = forcats::fct_reorder(infx, n), fill = type, color = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = ymax, label=paste0(round(prop*100,1), "% (n=",n,")")) , color="black",hjust = -0.1) +
  geom_errorbar(aes(xmin = ymin, xmax = ymax, y = infx), width = 0.5, color = "black") +
  xlab("Percentage of Infections") + ylab("Infection Composition") +
  theme_bw() +
  ggpubr::theme_pubclean() +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  MetBrewer::scale_color_met_d("Egypt", name = "Infection Type") +
  MetBrewer::scale_fill_met_d("Egypt", name = "Infection Type") +
  theme(axis.line = element_line(), legend.position = c(0.8,0.3), panel.grid.major.y = element_blank())

save_figs("species_composition", fig1, width = 10, height = 5)
