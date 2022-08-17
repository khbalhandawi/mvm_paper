library(ggplot2)
library(dplyr)
library(ggthemes)

#HOME_WD <- "~/"
HOME_WD <- Sys.getenv("GITDIR")
main_wd <- paste0(HOME_WD,"/mvm_paper/")

# For Icord paper example
data_dir <- paste0("strut_manufacturability/")
img_dir <- paste0("images/strut_manufacturability/")

## CHANGE TO MAIN WD
setwd(main_wd)
source(paste0(main_wd,"/plot_funcs.R"))

export_theme <- theme_tufte() +
  theme(
    axis.text.x=element_text(size=12,family=""),
    axis.text.y=element_text(size=12,family=""),
    axis.title.x=element_text(size=14,family="",vjust=-1),
    axis.title.y=element_text(size=14,family=""),
    
    ## Axis lines
    axis.line = element_line(colour="black", size = 0.1),
    axis.ticks = element_line(size=0.1),
    
    ## Title
    plot.title = element_text(family="",size=10,face="bold",hjust=0.5),
    plot.tag = element_text(family="",size=12,face="bold"),
    
    ## Legends
    legend.title=element_text(size=14,family="",face="italic"),
    legend.text=element_text(size=12,family=""),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    
    ## Strips for facet_wrap
    strip.text=element_text(size=14,family="",face="bold"),
    strip.background=element_rect(fill="#f0f0f0",size=0.1)
    # strip.background=element_blank()
  )


df_matrix <- read.csv(paste0(main_wd,data_dir,"/impact_data.csv"))
df_matrix <- df_matrix %>% 
  mutate(Margin = as.factor(Margin)) # convert to categorical type

# Draw barplot with grouping & stacking
plot <- ggplot(df_matrix, aes(x = Concept, y = Impact, fill = Margin)) + 
  geom_bar(stat = "identity", width = 0.9, position = "stack") +
  coord_flip() +
  facet_wrap(~ facet, ncol = 1) +
  export_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position="none")
plot

# Draw barplot with grouping & stacking
df_inset <- df_matrix[df_matrix$facet=='manuf.',]
plot_inset <- ggplot(df_inset, aes(x = Concept, y = Impact, fill = Margin)) + 
  geom_bar(stat = "identity", width = 0.9, position = "stack") +
  export_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
plot_inset
########################################
## 2. Export plots
########################################
ggsave(paste0(main_wd,img_dir,"impact_bars.pdf"),
       plot = plot, height=4, width=6.5)

ggsave(paste0(main_wd,img_dir,"impact_inset.pdf"),
       plot = plot_inset, height=4, width=3)