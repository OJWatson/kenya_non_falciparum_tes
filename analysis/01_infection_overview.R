library(tidyverse)
our_theme <- ggplot2::theme(plot.title = element_text(hjust = 0.5),
                            panel.border = element_blank(),
                            axis.line = element_line(colour = "black"))


# ------------------------------------------------------------------------------
# Fig 1
# ------------------------------------------------------------------------------

# read in raw data
df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))

# summarise by species
df_summary_species_mixed <- df %>%
  group_by(species) %>%
  summarise(mean = length(sample_id)/nrow(df),
            n = length(sample_id)) %>%
  as.data.frame()

# mixed infection composition
countCharOccurrences <- function( s,char) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

# summarise the data for plotting
df_summary_species_mixed$type <- as.factor(sapply(df_summary_species_mixed$species,countCharOccurrences,"/")+1)
df_summary_species_mixed <- rbind(df_summary_species_mixed,c("Pm/Poc/Pow",0,0,3))
df_summary_species_mixed$mean <- as.numeric(df_summary_species_mixed$mean)
df_summary_species_mixed$type <- c("Single","Double","Triple","Quadruple")[df_summary_species_mixed$type]
df_summary_species_mixed$type <- factor(df_summary_species_mixed$type ,levels=c("Single","Double","Triple","Quadruple"))
df_summary_species_mixed$ci025 <- vapply(as.numeric(df_summary_species_mixed$n), function(x) Hmisc::binconf(x,2027)[2], numeric(1))
df_summary_species_mixed$ci975 <- vapply(as.numeric(df_summary_species_mixed$n), function(x) Hmisc::binconf(x,2027)[3], numeric(1))

# fig 1 plot
mixed_infection <- ggplot(df_summary_species_mixed[rev(order(df_summary_species_mixed$mean)),],
                          aes(x=reorder(species,mean),y=mean,fill=type,color=type)) + geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(mean,4)*100,"% (",n,"/2027)")) , color="black",hjust = -0.4) +
  geom_errorbar(aes(ymin = ci025, ymax = ci975, x = species), color = "black", width = 0.25) +
  ylab("Percentage Infections") + xlab("Infection Composition") +
  theme_bw() + our_theme + coord_flip() +  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_color_brewer(palette = "Paired",name="Infection Type",) +
  scale_fill_brewer(palette = "Paired",name="Infection Type")

save_figs("fig1", mixed_infection, 9.5, 6)

df_summary_species_mixed$mean <- scales::percent(df_summary_species_mixed$mean, accuracy = 0.01)
df_summary_species_mixed$ci025 <- scales::percent(df_summary_species_mixed$ci025, accuracy = 0.01)
df_summary_species_mixed$ci975 <- scales::percent(df_summary_species_mixed$ci975, accuracy = 0.01)
df_summary_species_mixed <- df_summary_species_mixed[order(as.numeric(df_summary_species_mixed$n), decreasing = TRUE),]
df_summary_species_mixed$n <- paste0(df_summary_species_mixed$n, "/", sum(as.numeric(df_summary_species_mixed$n)))
df_summary_species_mixed$mean <- paste0(df_summary_species_mixed$mean, " [", df_summary_species_mixed$ci025, " - ", df_summary_species_mixed$ci975, "]")
df_summary_species_mixed <- df_summary_species_mixed[,c("species", "type", "mean", "n")]
write.csv(df_summary_species_mixed, file.path(here::here(), "analysis/tables/figure1_table.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Supp Figure By Sites
# ------------------------------------------------------------------------------

# summarise by sites
df_summary_sites <- df %>%
  group_by(site) %>%
  mutate(n=n()) %>%
  group_by(site,species) %>%
  summarise(mean = length(sample_id)/unique(n),
            n=n()) %>%
  ungroup() %>%
  arrange(site, (mean)) %>%
  mutate(.r = row_number())

df_summary_sites$type <- as.factor(sapply(df_summary_sites$species,countCharOccurrences,"/")+1)
df_summary_sites$type <- c("Single","Double","Triple","Quadruple")[df_summary_sites$type]
df_summary_sites$type <- factor(df_summary_sites$type ,levels=c("Single","Double","Triple","Quadruple"))

# create our plot by sites
mixed_infection_sites <- ggplot(df_summary_sites[rev(order(df_summary_sites$mean)),],
                                aes(icer:::reorder_within(species, mean, site),y=mean,fill=type,color=type)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(mean,4)*100,"% (n=",n,")")) , color="black",hjust = -0.1) +
  ylab("Percentage Infections") + xlab("Infection Composition") +
  theme_bw() + our_theme + coord_flip() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14), strip.text = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size=12),
        panel.background = element_rect(colour="black")) +
  facet_wrap(~site, scales="free",ncol = 2) +
  icer:::scale_x_reordered() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_color_brewer(palette = "Paired",name="Infection Type",) +
  scale_fill_brewer(palette = "Paired",name="Infection Type")

# save to results
save_figs("species_sites",mixed_infection_sites, 16, 16)
