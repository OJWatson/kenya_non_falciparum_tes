library(tidyverse)

df2 <- readxl::read_xlsx(cp_path("analysis/data/raw/vivo_scheduled_visits.xlsx"), sheet = 1)
df2 <- df2 %>%
  rename(subjectid = `SUBJECT ID`) %>%
  rename(D0 = `D0-D3`)

df3 <- readxl::read_xlsx(cp_path("analysis/data/raw/other_scheduled_visits.xlsx"), sheet = 1)
df3 <- df3 %>%
  rename(subjectid = `SUBJECT ID`) %>%
  rename(D0 = `DO-D3`) %>%
  select(-OUTCOME)

dfnew <- rbind(df2, df3)
dfnew <- dfnew %>% mutate(subjectid = gsub("-", "", subjectid))
saveRDS(dfnew, cp_path("analysis/data/derived/slide_data.rds"))
