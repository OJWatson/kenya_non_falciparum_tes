library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Read in our daily data and format
# ------------------------------------------------------------------------------

# Read in our raw data
df <- readxl::read_xlsx(cp_path("analysis/data/raw/raw_data.xlsx"), sheet = 3)
df <- select(df, subjectid:gender)

# Species labels
specs <- c("Pf","Pm","Poc","Pow")

pos_neg_to_int <- function(x) {
  x[x == "Pos"] <- "1"
  x[x == "Neg"] <- "0"
  x <- as.logical(as.integer(x))
  x
}


# pivot and format for plotting
df2 <- left_join(
  df %>%
  mutate(across(starts_with(specs), .fns = pos_neg_to_int)) %>%
  select(-starts_with(c("hgb","plt","ne"))) %>%
  pivot_longer(starts_with(specs),
               names_to = c("species", "day"),
               names_pattern = "([A-Za-z]+)(\\d+)",
               values_to = "infected"),
  df %>%
    select(-starts_with(specs)) %>%
    pivot_longer(starts_with(c("hgb","plt","ne")),
                 names_to = c("vitals", "day"),
                 names_pattern = "([A-Za-z]+)(\\d+)",
                 values_to = "value")
) %>%
  pivot_wider(names_from = vitals, values_from = value) %>%
  pivot_wider(names_from = species, values_from = infected) %>%
  mutate(
    infx = purrr::pmap_chr(list(pf, pm, poc, pow), ~ paste0(specs[c(..1,..2,..3,..4)],collapse="/")),
    infx = replace(infx, which(infx == ""), "Uninfected"),
    infx = replace(infx, which(infx == "NA/NA/NA/NA"), "Missing"),
    infx = factor(infx, levels = c(
      gen_combinations(specs),
      "Uninfected", "Missing"
      )),
    day = as.integer(day)
  )
saveRDS(df2, cp_path("analysis/data/derived/daily_data.rds"))

# ------------------------------------------------------------------------------
# 2. Read in our hourl data and format
# ------------------------------------------------------------------------------

dat <- haven::read_dta("analysis/data/raw/cbcmbfvivo1n2_19_Nov_2021.dta")

dat <- dat %>%
  mutate(time = as.numeric(substr(sample_id, 5,6))) %>%
  mutate(day = as.numeric(substr(sample_id, 2, 3 ))) %>%
  mutate(time = replace(time, which(is.na(time)), day[is.na(time)]*24))

saveRDS(dat, cp_path("analysis/data/derived/hourly_data.rds"))


