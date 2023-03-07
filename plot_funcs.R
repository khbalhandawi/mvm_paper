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

# for ggpairs labeling
charcount <- function(x, w, fx = TRUE) {
  nchar(x) - nchar(gsub(w, "", x, fixed = fx))
}


look_for_tex <- function(x) {
  if (getOption("pmplots_TeX_labels", FALSE)) {
    return(TRUE)
  }
  charcount(x, "$") >= 2
}


parse_label <- function(x) {
  if (substr(x, 1, 2) == "!!") {
    x <- parse(text = substr(x, 3, nchar(x)))
    return(x)
  }
  if (look_for_tex(x)) {
    if (requireNamespace("latex2exp")) {
      return(latex2exp::TeX(x))
    }
  }
  x
}

label_parse_label <- function(x) {
  x <- lapply(x, as.character)
  lapply(x, function(values) {
    lapply(values, parse_label)
  })
}

# process pdf data into cdf with intervals
process_cdf <- function(data, df_maxmin, col, fcol) {
  
  data = data %>% arrange(col)
  
  df_gmin = data %>% 
    group_by_at(fcol) %>% 
    summarise_at(.vars=col, .funs=min)
  
  df_gmax = data %>% 
    group_by_at(fcol) %>% 
    summarise_at(.vars=col, .funs=max)
  
  df_interval = merge(df_gmin, df_gmax, by="decision", all.x=TRUE) %>% 
    arrange_at(paste0(col,".x")) %>% 
    mutate(start=0.0) %>% 
    mutate(end=0.0)
  df_interval$design_ID = 1:nrow(df_interval)
  
  for (row in 1:nrow(df_interval)) {
    if (row == 1) {
      df_interval[row, "start"] = df_maxmin[col][1]
    } else {
      df_interval[row, "start"] = min(c(df_interval[row-1, paste0(col,".y")],
                                        df_interval[row, paste0(col,".x")]))
    }
    if (row == nrow(df_interval)) {
      df_interval[row, "end"] = df_interval[row, paste0(col,".y")]
    }
    else {
      df_interval[row, "end"] = min(c(df_interval[row, paste0(col,".y")],
                                      df_interval[row+1, paste0(col,".x")]))
    }
  }
  interval_vector = c(df_interval$start,df_interval$end[nrow(df_interval)])
  
  P = ecdf(data[,col])
  
  # Create interval df
  interval_limits = P(interval_vector)
  marks = data.frame(interval_vector,interval_limits)
  marks$design_ID = 1:nrow(marks)
  marks = merge(marks, df_interval[,c("design_ID","decision")], by = "design_ID", all.x = TRUE)
  marks$decision = factor(marks$decision,levels=levels(data[,fcol]))
  
  # Create cdf df
  x_cdf = seq(as.numeric(df_maxmin[1,col]), as.numeric(df_maxmin[2,col]), length.out = 1000)
  y_cdf = P(x_cdf)
  df_cdf = data.frame(x_cdf,y_cdf) %>% 
    mutate(f_cdf = findInterval(x_cdf, interval_vector, all.inside=FALSE, rightmost.closed=TRUE))  
  df_cdf = merge(df_cdf, df_interval[,c("design_ID","decision")], by.x = "f_cdf", by.y = "design_ID",all.x = TRUE)
  df_cdf$decision = factor(df_cdf$decision,levels=levels(data[,fcol]))
  
  return(list(marks,df_cdf))
}

# get legend labels
process_labels <- function(labels,inconel=6) {
  labels_plot <- c()
  for (i in 1:length(labels)) {
    label <- c()
    label_og <- unlist(strsplit(labels[i], ","))
    for (j in 1:length(label_og)) {
      if (j!=2) {
        if (j==1) {
          label <- c(label,paste0(as.character(i),": ",label_og[j]))
        } else {
          label <- c(label,label_og[j])
        }
      } else {
        if (label_og[j] < inconel) {
          label <- c(label,"steel")
        } else if (label_og[j] == inconel) {
          label <- c(label,"Inconel")
        } else {
          label <- c(label,"titanium")
        }
      }
    }
    labels_plot <- c(labels_plot,paste(label, collapse = ','))
  }
  return(labels_plot)
}

# The paper's plot theme and color palette
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
extend_palette <- function(n_nodes,n_concepts,palette,show=FALSE) {
  
  n_concepts = rep(n_concepts,n_nodes)
  n_nodes = seq(1,n_nodes,by=1)
  df_concepts = data.frame(n_nodes, n_concepts)
  
  colors_list = list()
  colors_vector = c()
  
  for (i in 1:nrow(df_concepts)) {
    alphas = seq(1.0, 0.2, length.out = df_concepts$n_concepts[i])
    region_c = c()
    for (j in 1:df_concepts$n_concepts[i]) {
      sub_color = adjustcolor(palette[i], alpha.f = alphas[j])
      region_c = append(region_c,sub_color)
    }
    colors_list = append(colors,list(region_c))
    colors_vector = append(colors_vector,region_c)
    
  }
  if (show) {scales::show_col(colors_vector)}
  return(colors_vector)
}