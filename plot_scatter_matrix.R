## IMPORTANT NOTE: Code to visualize the scatter pairs of impact and absorption
## matrices of FEA strut problem

########################################
## 1. Headers
########################################library(ggthemes)
library(GGally)
library(ggpubr)
library(patchwork)
library(tibble)
library(ggthemes)
library(dplyr)
library(latex2exp)

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
    legend.title=element_text(size=14,family="",face="italic"),
    legend.text=element_text(size=14,family=""),
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
my_dens <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
  geom_density(...,alpha = 0.5, color = "black")
}

my_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
  geom_histogram(...,alpha = 0.5, color = "black", size=0.1)
}

# Defines function to color according to correlation
# https://stackoverflow.com/a/57455317
cor_func <- function(data, mapping, method, symbol, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  corr <- cor(x, y, method=method, use="complete.obs")
  colFn <- colorRampPalette(c("brown1", "white", "dodgerblue"),
                            interpolate ="spline")
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]

  ggally_text(
    label = paste(symbol, as.character(round(corr, 2))),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = "black",
    ...
  ) + #removed theme_void()
  theme(panel.background = element_rect(fill = fill))
}

df_matrix <- read.csv(paste0(main_wd,data_dir,"/df_matrix.csv"))
df_matrix <- df_matrix %>%
  mutate(node = as.factor(node)) # convert to categorical type

p_pairs <- ggpairs(df_matrix,
                   mapping = ggplot2::aes(color=node, fill=factor(node, labels = c(expression(e1), expression(e2), expression(e3)))),
                   columns=3:ncol(df_matrix),
                   diag = list(continuous = my_hist),
                   upper = list(continuous = GGally::wrap(ggally_cor, stars = F)),
                   # upper = list(continuous = wrap(cor_func,method = "spearman", symbol = expression("\u03C1 ="))),
                   lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5)),
                   legend=1,
                   columnLabels = c("$W$","$c_{raw}$","$T_1$","$T_2$","$B_x$","$B_y$"),
                   labeller = label_parse_label
                  ) +
  labs(fill = "margin node") +
  export_theme +
  theme(legend.position = "bottom")
p_pairs

########################################
## 2. Scatter plot matrix of a node
########################################

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

df_total <- read.csv(paste0(main_wd,data_dir,"/df_total.csv"))

# get limits
df_maxmin <- bind_rows(apply(df_total,2,min),apply(df_total,2,max))

# normalize impact and absorption
df_total <- df_total %>% 
  # mutate(across(8:ncol(df_total), normalize)) %>%
  mutate(across(6:8, as.factor))

df_total$decision <- as.factor(with(df_total,paste(D1,D2,D3,sep=",")))

# get legend labels
labels <- levels(df_total$decision)
labels_plot <- process_labels(labels)

ID_levels = 1:length(levels(df_total$decision))
ID <- data.frame(ID_levels,levels(df_total$decision))
colnames(ID) <- c("ID","decision")
df_total <- merge(df_total, ID, by = "decision", all.x = TRUE) %>%
  mutate(decision = as.factor(ID))

# df_total %>% expand(nesting(D1, D2, D3))

my_cdf <- function(data, mapping, ...) {
  col = as.character(mapping$x)[2]
  fcol = as.character(mapping$colour)[2]
  data = data %>% arrange(col)
  
  df_gmin = data %>% 
    group_by_at(fcol) %>% 
    summarise_at(.vars=col, .funs=min)
  
  df_gmax = data %>% 
    group_by_at(fcol) %>% 
    summarise_at(.vars=col, .funs=max)
  
  P = ecdf(data[,col])
  x_cdf = data[,col]
  y_cdf = P(x_cdf)
  f_cdf = data[,fcol]
  df_cdf = data.frame(x_cdf,y_cdf,as.factor(f_cdf))
  
  ggplot(data = df_cdf, mapping=aes(x=x_cdf,y=y_cdf,color=f_cdf,fill=f_cdf)) +
    geom_line() +
    geom_area(position="identity") +
    scale_x_continuous(limits=c(as.numeric(df_maxmin[1,col]),
                                as.numeric(df_maxmin[2,col])))
}

my_cdf_stacked <- function(data, mapping, ...) {
  col = as.character(mapping$x)[2]
  fcol = as.character(mapping$colour)[2]
  
  if (grepl("S", col, fixed = TRUE)) {
    plot = ggplot(data = data, mapping=mapping) +
      # geom_density(...,alpha = 0.5, color = "black")
      geom_histogram(...,alpha = 0.5, color = "black", size=0.1)
  } else {
    dfs = process_cdf(data,df_maxmin,col,fcol)
    marks = dfs[[1]]
    df_cdf = dfs[[2]]
    # plot everything
    plot = ggplot(data = df_cdf, mapping=aes(x=x_cdf,y=y_cdf,color=decision,fill=decision)) +
      geom_line() +
      geom_segment(data=marks, aes(x=interval_vector, xend=interval_vector,y=0.0, yend=interval_limits,color=NA)) +
      geom_area(position="identity",alpha=0.5) +
      scale_x_continuous(limits=c(as.numeric(df_maxmin[1,col]),
                                  as.numeric(df_maxmin[2,col])))
  }
  return(plot)
}

my_corr <- function(data, mapping, ...) {
  xcol = as.character(mapping$x)[2]
  ycol = as.character(mapping$y)[2]
  fcol = as.character(mapping$colour)[2]
  
  ggally_cor(data, mapping = aes_string(x = xcol, y = ycol),color="black",size = 3,stars=F)
}

# make the plot
p_node <- ggpairs(df_total,
                  mapping = ggplot2::aes(color=decision, fill=factor(decision)),
                  columns=c("decision","S1","S2","S3","S4","I3","A3"),
                  diag = list(continuous = my_cdf_stacked),
                  upper = list(continuous = my_corr, combo = "box_no_facet"),
                  # upper = list(continuous = wrap(cor_func,method = "spearman", symbol = expression("\u03C1 ="))),
                  lower = list(continuous = wrap("points", alpha = 0.2, size = 0.8), 
                               combo = wrap("dot_no_facet", alpha = 0.4, size = 0.8)),
                  columnLabels = c("decision","$T_1$","$T_2$","$B_x$","$B_y$","$I_3$","$A_3$"),
                  labeller = label_parse_label,
                  legend=1
                  ) +
  export_theme +
  scale_colour_discrete(guide = "none") +
  scale_fill_discrete(labels = labels_plot) +
  labs(fill = unname(TeX("$w$, material, $n_{struts}$"))) +
  export_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10))
p_node

########################################
## 3. Export plots
########################################
ggsave(paste0(main_wd,img_dir,"scatter_pairs.pdf"),
       plot = p_pairs, height=6, width=9)

ggsave(paste0(main_wd,img_dir,"scatter_pairs_E3.pdf"),
       plot = p_node, height=6, width=9)
