library(tidyverse)

# create table one from our data
daily <- readRDS(cp_path("analysis/data/derived/daily_data.rds"))
daily <- daily %>% filter(day == 0) %>%
  mutate(pf_only = infx == "Pf")

# write these results to table one
write.csv(
daily %>%
  select(gender, age, weight, hgb, plt, gmeanpfemia, pf_only, day28, day42, slope_half_life, pc50, pc99, arm) %>%
  tableone::CreateTableOne(strata = "arm" , data = .) %>%
  print(., nonnormal = c("age", "weight", "gmeanpfemia"), formatOptions = list(big.mark = ","), noSpaces = TRUE),
cp_path("analysis/tables/table_one.csv")
)

# evidence plot for why normality chosen
normality_plot <- daily %>%
  select(subjectid, gender, age, weight, hgb, plt, gmeanpfemia, day28, day42, slope_half_life, pc50, pc99, arm) %>%
  rename(Age = age, Parasitaemia= gmeanpfemia, Haemoglobin = hgb, Platelets = plt,
         PCE_Slope_Half_Life = slope_half_life, PCE_PC50 = pc50, PCE_PC99 = pc99) %>%
  pivot_longer(-c("subjectid","arm","gender", "day28","day42")) %>%
  ggplot(aes(x=value, fill = arm)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~name, scales = "free") +
  scale_fill_viridis_d(name = "ACT") +
  theme_bw()
save_figs("normality_plot", normality_plot, 10, 8)
