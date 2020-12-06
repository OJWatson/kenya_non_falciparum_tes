library(malariaAtlas)

# Get the Kenya prevalence map
KEN_shp <- getShp(ISO = "KEN", admin_level = "admin0")
KEN_PfPR2_10 <- getRaster(surface = "Plasmodium falciparum PR2-10", shp = KEN_shp, year = 2015)
KEN_PfPR2_10_df <- as.MAPraster(KEN_PfPR2_10)
KEN_shp_df <- as.MAPshp(KEN_shp)
p <- autoplot(KEN_PfPR2_10_df, shp_df = KEN_shp_df, printed = FALSE)
pr <- getPR(country = c("Kenya"), species = "Pf")

# Adding in prevalence
coords <- data.frame(rbind(
  "KOM" = c(-0.103434, 34.518362),
  "KSI" = c(-0.6823459, 34.7080638),
  "KDH" = c(-0.131087, 34.8044322),
  "MDH" = c(-3.2266136, 40.1226615),
  "KCH" = c(-0.3727785, 35.2757856),
  "MGT" = c(0.4703683, 35.9795809)
))

hospitals <- c(
  "Kombewa" = "KOM",
  "Kisii" = "KSI",
  "Kisumu" = "KDH",
  "Malindi" = "MDH",
  "Kericho" = "KCH",
  "Marigat" = "MGT"
)
coords$name <- names(hospitals)

# fill in map
map <- p[[1]] +
  geom_point(data = coords, aes(X2, X1),fill="black", shape = 21)+
  ggrepel::geom_label_repel(data = coords, aes(X2, X1, label=name), min.segment.length = 0) +
  theme(plot.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_distiller(name = "PfPR", palette = "RdYlBu")

save_figs("pfpr_map", map, height = 6, width = 6)
