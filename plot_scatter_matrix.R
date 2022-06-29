## IMPORTANT NOTE: Code to visualize the scatter pairs of impact and absorption 
## matrices of FEA strut problem

########################################
## 1. Headers
########################################
library(ggthemes)
library(GGally)
library(ggpubr)
library(patchwork)
library(tibble)

#HOME_WD <- "~/"
HOME_WD <- Sys.getenv("GITDIR")
## Where the MCMC chains are stored
main_wd <- paste0(HOME_WD,"/mvm_paper/")

# For fea example
data_dir <- paste0("data/strut_fea/C1/")
img_dir <- paste0("images/strut_fea/C1/")

# # For simple example
# data_dir <- paste0("data/strut/C1/")
# img_dir <- paste0("images/strut/C1/")

## CHANGE TO MAIN WD
setwd(main_wd)
source(paste0(main_wd,"/plot_funcs.R"))

export_theme <- theme_tufte() +
  theme(
    axis.text.x=element_text(size=8,family=""),
    axis.text.y=element_text(size=8,family=""),
    axis.title.x=element_text(size=12,family="",vjust=-1),
    axis.title.y=element_text(size=12,family=""),
    
    ## Axis lines
    axis.line = element_line(colour="black", size = 0.1),
    axis.ticks = element_line(size=0.1),
    
    ## Title
    plot.title = element_text(family="",size=10,face="bold",hjust=0.5),
    plot.tag = element_text(family="",size=12,face="bold"),
    
    ## Legends
    legend.title=element_text(size=10,family="",face="italic"),
    legend.text=element_text(size=10,family=""),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    
    ## Strips for facet_wrap
    strip.text=element_text(size=10,family="",face="bold"),
    strip.background=element_rect(fill="#f0f0f0",size=0.1)
    # strip.background=element_blank()
  )

########################################
## 1. Scatter plot matrix
########################################

# https://stackoverflow.com/a/34988306
my_dens <- function(data, mapping, ...){
  ggplot(data = data, mapping=mapping) +
    geom_density(...,alpha = 0.5, color = "black")
}

my_hist <- function(data, mapping, ...){
  ggplot(data = data, mapping=mapping) +
    geom_histogram(...,alpha = 0.5, color = "black", size=0.1)
}

# Defines function to color according to correlation
# https://stackoverflow.com/a/57455317
cor_func <- function(data, mapping, method, symbol, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  colFn <- colorRampPalette(c("brown1", "white", "dodgerblue"), 
                            interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]

  ggally_text(
    label = paste(symbol, as.character(round(corr, 2))), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  ) + #removed theme_void()
    theme(panel.background = element_rect(fill = fill))
}
  
df_matrix <- read.csv(paste0(main_wd,data_dir,"/df_matrix.csv"))
df_matrix <- df_matrix %>% 
  mutate(node = as.factor(node)) # convert to categorical type

p_pairs <- ggpairs(df_matrix, mapping = ggplot2::aes(colour=node, fill = node), 
                   columns=3:ncol(df_matrix), 
                   diag = list(continuous = my_hist),
                   upper = list(continuous = GGally::wrap(ggally_cor, stars = F)),
                   # upper = list(continuous = wrap(cor_func,method = 'spearman', symbol = expression('\u03C1 ='))),
                   lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5)),
                   legend=1) +
  export_theme +
  theme(legend.position = "bottom")


########################################
## 2. Export plots
########################################
ggsave(paste0(main_wd,img_dir,"scatter_pairs.pdf"),
       plot = p_pairs, height=6, width=9)