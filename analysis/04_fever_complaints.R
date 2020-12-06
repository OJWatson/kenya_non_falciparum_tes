library(tidyverse)

# ------------------------------------------------------------------------------
# Read In Raw Data
# ------------------------------------------------------------------------------

# Read in our raw data
df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))

# remove missing fields
df <- df[!is.na(df$age),]
df <- df[!is.na(df$fever),]

# Helper and Plotter Functions
# ------------------------------------------------------------------------------

# Mean and CI boot function
mean_and_ci <- function(x,na.rm=TRUE, perc=FALSE,nl=2){
  if(na.rm) {
    x <- x[!is.na(x)]
  }

  bm <- function(x, i){
    mean(x[i])
  }

  boot.mean = function(x,B=2000,binwidth=NULL) {
    n = length(x)
    boot <- boot::boot(x, bm, R = B)
    boot.ci <- boot::boot.ci(boot, type = "all")
    #print( interval )
    return( list(boot.statistics = boot.ci, interval=boot.ci$bca[4:5]) )
  }

  rsp <- function(z,d=2){
    out <- sprintf("%.2f",round(z,d))
    if(nchar(out) < (nl*2 + 1)){
      out <- paste0(rep(" ",(nl*2 + 1)-nchar(out)),out)
    }
    return(out)
  }

  y <- boot.mean(x,length(x+500))
  if(perc) {
    return(paste0(rsp(mean(x*100),2),"%",
                  " [",
                  paste0(rsp(y$interval[1]*100,2),"-",rsp(y$interval[2]*100,2)),
                  "]"))
  } else {
    return(paste0(rsp(mean(x),2),
                  " [",
                  paste0(rsp(y$interval[1],2),"-",rsp(y$interval[2],2)),
                  "]"))
  }
}

# ------------------------------------------------------------------------------
# Fig 4
# ------------------------------------------------------------------------------

# Construct model for Fever presentation in individuals with Pf
# controlling for variation between endimicity sites
tru_fev <- lme4::glmer(fever ~ 0 + pow + poc + pm + malaria_attacks_in_last_year + year_time + age + female + (1|endemicity_sites),
                       data = df[df$pf==1,], family = "binomial")

# Odds Ratio plot
or_plot <- function(glm_fev){

  # Tidy up mixed model data to be ORs from binom model
  td <- broom.mixed::tidy(glm_fev,conf.int=TRUE,) %>%
    mutate_at(vars(matches("est|conf")),.funs = exp)

  # label our covariates more clearly
  td$term[td$term == "malaria_attacks_in_last_year"] <- "Number of malaria\nattacks in last year"
  td$term[td$term == "treated_6"] <- "Treated for malaria\nin last 6 months"
  td$term[td$term == "year_time"] <- "Year"
  td$term[td$term == "multiply_infected"] <- "Multiply Infected"
  td$term[td$term == "female"] <- "Female"
  td$term[td$term == "age"] <- "Age (years)"
  td$term[td$term == "poc"] <- "Coinfection with\nP. ovale curtisi"
  td$term[td$term == "pm"] <- "Coinfection with\nP. malariae"
  td$term[td$term == "pow"] <- "Coinfection with\nP. ovale wallikeri"

  # format the labels and title
  td$label <- paste0(round(td$estimate,2), " (",
                     round(td$conf.low,digits = 2), "-",
                     round(td$conf.high,digits = 2), ")")
  td <- rbind(td, td[1,])
  td[nrow(td),2:ncol(td)] <- NA
  td$term[nrow(td)] <- ""
  td$label[nrow(td)] <- "Odds Ratio (95% CI)"

  # order covariates by length of label for easy plotting
  td <- td[order(nchar(td$term),decreasing = TRUE),]
  td$term <- factor(td$term, levels = td$term)
  if("effect" %in% names(td)) {
    td <- td[td$effect == "fixed",]
  }

  # create OR plot
  gg <- ggplot(td, aes(estimate, term, color = p.value<0.05)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height=0.25)) +
    geom_vline(xintercept = 1) + theme_bw() +
    scale_color_manual(name="p-value < 0.05",values = c("#002366","#e62400")) +
    theme(text = element_text(size = 14), panel.spacing = unit(10,units = "pt"),) +
    xlab("Odds Ratio") + ylab("") +
    geom_text(inherit.aes = FALSE, data = td,
              mapping = aes(y = term, x = 2, label = label)) +
    theme(legend.position = "none",
          panel.grid.major.y =  element_blank(),
          panel.grid.minor.y =  element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(hjust = 0.4)) +
    scale_x_continuous(breaks = c(0.5, 1, 1.5), limits = c(0,2.51))

  # Add OR desc to right
  gg + geom_curve(data = data.frame(x = 1.5, y = nrow(td)-0.25,
                                    xend = 2.5, yend = nrow(td)-0.25),
                  mapping = aes(x = x,y = y, xend = xend, yend = yend),
                  angle = 0L, arrow = arrow(30L, unit(0L, "inches"), "last", "closed"),
                  alpha = 1, inherit.aes = FALSE)

}

# Create our plot
or <- or_plot(tru_fev)

save_figs("fig4", or, width = 10, height = 5)

# ------------------------------------------------------------------------------
# Supp Fever
# ------------------------------------------------------------------------------

supp_fever <- df %>%
  mutate(age_bin = cut(age,c(0,4,8,12,16,20,30,40,80), include.lowest = TRUE)) %>%
  group_by(age_bin) %>%
  summarise(
    ych = stringr::str_trim(mean_and_ci(fever,na.rm = TRUE)),
    xch = stringr::str_trim(mean_and_ci(age,na.rm = TRUE)),
    n = n()
  ) %>%
  mutate(
    y = vapply(str_split(ych,"-|]|\\["), function(x){as.numeric(x[1])}, numeric(1)),
    ymin = vapply(str_split(ych,"-|]|\\["), function(x){as.numeric(x[2])}, numeric(1)),
    ymax = vapply(str_split(ych,"-|]|\\["), function(x){as.numeric(x[3])}, numeric(1)),
    x = vapply(str_split(xch,"-|]|\\["), function(x){as.numeric(x[1])}, numeric(1)),
    xmin = vapply(str_split(xch,"-|]|\\["), function(x){as.numeric(x[2])}, numeric(1)),
    xmax = vapply(str_split(xch,"-|]|\\["), function(x){as.numeric(x[3])}, numeric(1))
  ) %>%
  ggplot(aes(x=x,y=y,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax,weight=n)) +
  geom_errorbar(width = 0) +
  geom_errorbarh(height = 0) +
  geom_point(aes(size = n)) +
  geom_smooth(span=2, method="lm") +
  theme_bw() +
  xlab("Age (years)") + ylab("Proportion of clinical presentation of fever")

save_figs("supp_fev", supp_fever, width = 10, height = 6)

# ------------------------------------------------------------------------------
# Table 1
# ------------------------------------------------------------------------------

df <- data.table::fread(file.path(here::here(), "analysis/data/raw/data.txt"))
df_sum <- df %>%
  group_by(site) %>%
  summarise(age=mean_and_ci(age,na.rm=TRUE),
            sex=mean_and_ci(female,perc = TRUE),
            true_fev=mean_and_ci(fever,na.rm=TRUE,perc=TRUE),
            mal_12=mean_and_ci(malaria_attacks_in_last_year,na.rm=TRUE),
            n=n())
df_sum_all <- df %>%
  summarise(site="All",
            age=mean_and_ci(age,na.rm=TRUE),
            sex=mean_and_ci(female,perc = TRUE),
            true_fev=mean_and_ci(fever,na.rm=TRUE,perc=TRUE),
            mal_12=mean_and_ci(malaria_attacks_in_last_year,na.rm=TRUE),
            n=n())
df_grouped <- rbind(df_sum, df_sum_all)
names(df_grouped) <- c("Site","Age","Sex (% female)","Fever",
                        "Malaria episodes in last 12 months","n")

sites <- c("Kisumu", "Kombewa", "Malindi", "Kericho", "Kisii", "Marigat")
regions <- c("Lake Endemic","Lake Endemic","Coastal Endemic","Highland Epidemic", "Highland Epidemic", "Semi-Arid Zone")
df_grouped$Region <- regions[match(df_grouped$Site, sites)]
df_grouped <- df_grouped[,c(1,ncol(df_grouped),2:(ncol(df_grouped)-1))]
write.table(df_grouped, file = file.path(here::here(), "analysis/tables/table1.txt"), sep = "\t",row.names = FALSE)

# ----------------------------
## COMPLAINT BY SPECIES
# ----------------------------

library(UpSetR)
spr <- spread(df,complaint,complaint,fill = 0) %>%
  mutate_at(na.omit(unique(df$complaint)[nchar(unique(df$complaint))>0]),function(x){as.numeric(nchar(x)>1)})

plot_l = list()
m_list <- names(which(table(df$complaint)>5))
for(i in seq_along(m_list)) {

  p <- upset(spr[spr[[m_list[i]]]==1,c("pf","poc","pow","pm")],sets.x.label = NULL,
             order.by = "freq",keep.order = TRUE,sets=c("pf","pm","poc","pow"),
             text.scale=1.25)
  t <- ggplot(data.frame(x=1,y=1,text=m_list[i]),aes(x=x,y=y,label=text)) +
    geom_label(size=4) + theme_void()

  r <- cowplot::plot_grid(p$Main_bar,NULL, p$Matrix,NULL,ncol=1,rel_heights = c(2,0.1,1,0.15),align="v")
  l <- cowplot::plot_grid(t, p$Sizes,ncol=1,rel_heights = c(2,1.1),align="v")
  plot_l[[i]] <-cowplot::plot_grid(l,r,ncol=2,align="v",rel_widths = c(1.7,3))
}

complaints <- cowplot::plot_grid(plotlist = plot_l)
save_figs("complaints", complaints, width = 10, height = 6)

