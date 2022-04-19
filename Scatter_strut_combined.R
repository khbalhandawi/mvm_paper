## IMPORTANT NOTE: Code to visualize the scatter pairs of impact and absorption 
## matrices of FEA strut problem

########################################
## 1. Headers
########################################
library(ggthemes)
library(GGally)
library(ggpubr)
library(patchwork)

#HOME_WD <- "~/"
HOME_WD <- Sys.getenv("GITDIR")
## Where the MCMC chains are stored
main_wd <- paste0(HOME_WD,"/mvm_paper/")

## CHANGE TO MAIN WD
setwd(main_wd)
source(paste0(main_wd,"/plot_funcs.R"))

export_theme <- export_theme + theme(plot.margin=unit(c(0,0,0,0),units="cm"),plot.tag=element_text(size=10,face="bold"))

data_wd <- paste0(main_wd,"/strut")

########################################
## 1. Scatter plot matrix
########################################

df_matrix <- read.csv(paste0(data_wd,"/df_matrix.csv"))
df_matrix <- df_matrix %>% 
  mutate(node = as.factor(node)) # convert to categorical type

p_pairs <- ggpairs(df_matrix, columns=3:ncol(df_matrix), ggplot2::aes(colour=node), legend=1) +
  export_theme + theme(legend.position = "bottom")

########################################
## 2. Scatter plot of mean
########################################

df_mean <- read.csv(paste0(data_wd,"/df_mean.csv"))
df_mean <- df_mean %>% 
  mutate(node = as.factor(node)) # convert to categorical type

plot1 <- ggplot(df_mean, aes(x = Absorption, y = Impact, color = node)) + 
  geom_point(aes(color = node), size = 3) + 
  geom_point(shape = 1, color = "black", size = 3) + 
  stat_smooth(method = "lm", fullrange = TRUE) +
  geom_rug() + 
  scale_y_continuous(name = "Change absorption capability") + 
  scale_x_continuous(name = "Impact on performance") + 
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9)) +
  export_theme

dens1 <- ggplot(df_mean, aes(x = Absorption, fill = node)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") +
  export_theme

dens2 <- ggplot(df_mean, aes(x = Impact, fill = node)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip() +
  export_theme

p_mean <- dens1 + plot_spacer() + plot1 + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))


########################################
## 3. Export plots
########################################
ggsave(paste0(main_wd,"/images/scatter_pairs.pdf"),
       plot = p_pairs,
       device = cairo_pdf,
       dpi = 320,height=6,width=8)

ggsave(paste0(main_wd,"/images/scatter_mean.pdf"),
       plot = p_mean,
       device = cairo_pdf,
       dpi = 320,height=6,width=8)