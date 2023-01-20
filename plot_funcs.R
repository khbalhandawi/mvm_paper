# Run once
# library(remotes)
# install.packages("extrafont")
# install.packages("cairo_pdf")
# remotes::install_version("Rttf2pt1", version = "1.3.8")
library(extrafont)
# # Install **TTF** Latin Modern Roman fonts from www.fontsquirrel.com/fonts/latin-modern-roman
# # Import the newly installed LModern fonts, change the pattern according to the
# # filename of the lmodern ttf files in your fonts folder
# font_import(pattern = "lmroman*")

# Run each time
library(extrafont)
loadfonts(device = "pdf")

# get default color palette
library(scales)
hex_codes1 <- scales::hue_pal()(3) # Identify hex codes
# "#F8766D" "#7CAE00" "#00BFC4" "#C77CFF"
# "#F8766D" "#00BA38" "#619CFF" "#C77CFF" (pastel)

# get legend function
library(gridExtra)
get_legend<-function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

export_theme <- theme_tufte() +
  theme(
    axis.text.x=element_text(size=14,family=""),
    axis.text.y=element_text(size=14,family=""),
    axis.title.x=element_text(size=14,family="",vjust=-1),
    axis.title.y=element_text(size=14,family=""),

    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),

    ## Title
    plot.title = element_text(family="",size=10,face="bold",hjust=0.5),
    plot.tag = element_text(family="",size=10,face="bold"),

    ## Legends
    legend.title=element_text(size=16,family="",face="italic"),
    legend.text=element_text(size=15,family=""),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="",face="bold"),
    #strip.background=element_rect(fill="#f0f0f0")
    strip.background=element_blank(),
  )

# create color palette for design concepts
n_concepts <- c(4, 4, 4)
node <- c(1,2,3)
df_concepts <- data.frame(node, n_concepts)

colors_list <- list()
colors_vector <- c()

# palette <- scales::hue_pal(h = c(0, 360))(nrow(df_concepts))
palette <- c("#F8766D","#00BA38","#619CFF")
# palette <- c("#F8766D","#00BFC4") for simple example
for (i in 1:nrow(df_concepts)) {
  alphas <- seq(1.0, 0.2, length.out = df_concepts$n_concepts[i])
  region_c <- c()
  for (j in 1:df_concepts$n_concepts[i]) {
    sub_color <- adjustcolor(palette[i], alpha.f = alphas[j])
    region_c <- append(region_c,sub_color)
  }
  colors_list <- append(colors,list(region_c))
  colors_vector <- append(colors_vector,region_c)
  # scales::show_col(region_c)
}