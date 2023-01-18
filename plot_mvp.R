## IMPORTANT NOTE: Code to visualize the scatter pairs of impact and absorption
## matrices of FEA strut problem

########################################
## 1. Headers
########################################
library(ggthemes)
# library(GGally)
library(ggpubr)
# library(patchwork)
library(tibble)
library(gridExtra)
library(grid)
library(tidyverse)

#HOME_WD <- "~/"
HOME_WD <- Sys.getenv("GITDIR")
## Where the working directory is
main_wd <- paste0(HOME_WD,"/mvm_paper/")

# For fea example
data_dir <- paste0("data/strut_fea/")
img_dir <- paste0("images/strut_fea/")
dir.create(file.path(img_dir), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(img_dir,"C1"), showWarnings = FALSE, recursive = TRUE)

legend_aes <- list(linetype=c(0,0,0,2), shape=c(16,16,16,NA))
legend_dist_aes <- list(linetype=c(2,2,2,2), shape=c(16,16,16,NA))
palette <- c("#F8766D","#00BA38","#619CFF","#000000")
thetas <- c(0.00, 30.0, 0.00, 30.0)
heights <- c(15.0,15.0,20.0,20.0)
# fea limits
x0 <- array(c(-0.009630757, 1.152203562))
x1 <- array(c(0.3329746, 4.3151886))

# # For simple example
# data_dir <- paste0("data/strut/")
# img_dir <- paste0("images/strut/")
# legend_dist_aes <- list(linetype=c(2,2,2), shape=c(16,16,NA))
# legend_aes <- list(linetype=c(0,0,2), shape=c(16,16,NA))
# palette <- c("#F8766D","#00BFC4","#000000")
# thetas = c(0.00, 30.0, 0.00, 30.0)
# heights = c(15.0,15.0,20.0,20.0)
# # simple limits
# x0 <- array(c(0.00, 2.5))
# x1 <- array(c(0.075, 3.5))

## CHANGE TO MAIN WD
setwd(main_wd)
source(paste0(main_wd,"/plot_funcs.R"))

export_theme <- export_theme + theme(plot.margin=unit(c(0,0,0,0),units="cm"),
                                     plot.tag=element_text(size=14,face="bold"),
                                     plot.tag.position = c(0.15, 0.95))

#######################################
# Scaling function
scaling <- function(x, lb, ub, operation) {
  if (operation == 1) {
    # scale
    x_out <- (x - lb) / (ub - lb)
    return(x_out)
  } else if (operation == 2) {
    # unscale
    x_out <- lb + x * (ub - lb)
    return(x_out)
  }
}

#######################################
# Get nearest point and length
nearest <- function(p1, p2, s) {
  x1 <- p1[1]
  y1 <- p1[2]
  x2 <- p2[1]
  y2 <- p2[2]
  xs <- s[1]
  ys <- s[2]
  dx <- x2 - x1
  dy <- y2 - y1
  n <- (p2 - p1) / norm((p2 - p1),"2")
  n_o <- array(c(-n[2], n[1]))
  qp <- s - p1
  d <- (qp %*% n_o)[1] # distance
  np <- s - (d * n_o) # nearest point
  # det <- dx * dx + dy * dy
  # a <- (dy * (ys - y1) + dx * (xs - x1)) / det
  # # calculate distance
  # d <- (crossprod(p2 - p1, s - p1)) / norm((p2 - p1),"2")
  # np <- array(c(x1 + a * dx, y1 + a * dy))
  return(list(np,d))
}

#######################################
# Get shortest path data frame
get_shortest_df <- function(df_mean) {
  n_nodes <- length(unique(df_mean$node))
  mean_nodes <- aggregate(df_mean[, 3:4], list(df_mean$node), mean)
  mean_nodes_n <- aggregate(df_mean[, 5:6], list(df_mean$node), mean)

  datalist <- list()
  for (node_i in 1:n_nodes) {
    s <- array(c(mean_nodes[node_i,2],mean_nodes[node_i,3]))
    nearest_data <- nearest(x0,x1,s)

    s_n <- array(c(mean_nodes_n[node_i,2],mean_nodes_n[node_i,3]))
    nearest_data_n <- nearest(x0_n,x1_n,s_n)

    dat <- data.frame(X = c(s[1],unlist(nearest_data[1])[1]),
                      Y = c(s[2],unlist(nearest_data[1])[2]),
                      X_n = c(s_n[1],unlist(nearest_data_n[1])[1]),
                      Y_n = c(s_n[2],unlist(nearest_data_n[1])[2]),
                      dist = as.numeric(nearest_data[2]),
                      dist_n = as.numeric(nearest_data_n[2]),
                      node = node_i)

    datalist[[node_i]] <- dat # add it to your list
  }

  shortest_line <- do.call(rbind, datalist) %>%
    mutate(node = as.factor(node)) # convert to categorical type

  return(shortest_line)
}

########################################
## 1. get max min limits and 45 deg line
########################################

data_wd_1 <- paste0(main_wd,data_dir,"C1/")
data_wd_2 <- paste0(main_wd,data_dir,"C2/")
data_wd_3 <- paste0(main_wd,data_dir,"C3/")
data_wd_4 <- paste0(main_wd,data_dir,"C4/")

df_mean_1 <- read.csv(paste0(data_wd_1,"df_mean.csv"))
df_mean_1 <- df_mean_1 %>%
  mutate(node = as.factor(node)) # convert to categorical type

df_mean_2 <- read.csv(paste0(data_wd_2,"df_mean.csv"))
df_mean_2 <- df_mean_2 %>%
  mutate(node = as.factor(node)) # convert to categorical type

df_mean_3 <- read.csv(paste0(data_wd_3,"df_mean.csv"))
df_mean_3 <- df_mean_3 %>%
  mutate(node = as.factor(node)) # convert to categorical type

df_mean_4 <- read.csv(paste0(data_wd_4,"df_mean.csv"))
df_mean_4 <- df_mean_4 %>%
  mutate(node = as.factor(node)) # convert to categorical type


i_min <- min(c(df_mean_1$Impact,df_mean_2$Impact,df_mean_3$Impact))
a_min <- min(c(df_mean_1$Absorption,df_mean_2$Absorption,df_mean_3$Absorption))
i_max <- max(c(df_mean_1$Impact,df_mean_2$Impact,df_mean_3$Impact))
a_max <- max(c(df_mean_1$Absorption,df_mean_2$Absorption,df_mean_3$Absorption))

# # automatic limits
# x0 <- array(c(i_min,a_min))
# x1 <- array(c(i_max,a_max))

x0_n <- array(c(0,0))
x1_n <- array(c(1,1))

#######################################
# scale if necessary
df_mean_1 <- df_mean_1 %>%
  mutate(Impact, Impact_n = (Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n = (Absorption - x0[2]) / (x1[2] - x0[2]))
df_mean_2 <- df_mean_2 %>%
  mutate(Impact, Impact_n = (Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n = (Absorption - x0[2]) / (x1[2] - x0[2]))
df_mean_3 <- df_mean_3 %>%
  mutate(Impact, Impact_n = (Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n = (Absorption - x0[2]) / (x1[2] - x0[2]))
df_mean_4 <- df_mean_4 %>%
  mutate(Impact, Impact_n = (Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n = (Absorption - x0[2]) / (x1[2] - x0[2]))

ref_line_data <- data.frame(X = c(i_min, i_max),
                            neutral = c(a_min, a_max),
                            X_n = c(0, 1),
                            neutral_n = c(0, 1))

########################################
## 2. Scatter plot of three concepts
########################################
p1 <- ggplot(df_mean_1, aes(x = Impact_n, y = Absorption_n, color = node)) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  # geom_point(shape = 1, color = "black", size = 2, stroke = 0.1, alpha = 0.1) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  labs(color = "margin node\n") +
  guides(colour = guide_legend(override.aes = legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="A")

p2 <- ggplot(df_mean_2, aes(x = Impact_n, y = Absorption_n, color = node)) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  guides(colour = guide_legend(override.aes = legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="B")

p3 <- ggplot(df_mean_3, aes(x = Impact_n, y = Absorption_n, color = node)) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  # geom_point(shape = 1, color = "black", size = 2, stroke = 0.1, alpha = 0.1) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  guides(colour = guide_legend(override.aes = legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="C")

p4 <- ggplot(df_mean_4, aes(x = Impact_n, y = Absorption_n, color = node)) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  # geom_point(shape = 1, color = "black", size = 2, stroke = 0.1, alpha = 0.1) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  guides(colour = guide_legend(override.aes = legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="D")

p1 <- p1 + theme(legend.position = "top")
legend <- get_legend(p1)

# Remove the legends from the plots
# t,r,b,l
p1 <- p1 + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank())
p2 <- p2 + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())
p3 <- p3 + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())
p4 <- p4 + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())

ylab <- textGrob("Change absorption capability", rot = 90, gp = gpar(fontfamily = "", size = 10, cex = 1.5))
xlab <- textGrob("Impact on performance", gp = gpar(fontfamily = "", size = 10, cex = 1.5))
p1_label <- textGrob(substitute(paste("(a) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[1])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))
p2_label <- textGrob(substitute(paste("(b) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[2])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))
p3_label <- textGrob(substitute(paste("(c) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[3])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))
p4_label <- textGrob(substitute(paste("(d) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[4])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))

p_concept <- grid.arrange(ylab, legend, p1, p2, p3, p4, xlab, ncol=3, nrow=4,
             layout_matrix = rbind(c(1,2,2), c(1,3,4), c(1,5,6), c(1,7,7)),
             widths = c(0.2, 2.7, 2.7), heights = c(0.2, 2.5, 2.5, 0.2))

# Aggregate plots


########################################
## 3. Export distance plots
########################################
shortest_line <- get_shortest_df(df_mean_1)
distance <- aggregate(shortest_line[, 5], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
distance_n <- aggregate(shortest_line[, 6], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
dist_1 <- sum(distance$dist)
dist_1_n <- sum(distance_n$dist)

p1_dist <- ggplot(shortest_line, aes(x = X_n, y = Y_n, color = node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size = 3.5) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  labs(color = "margin node\n") +
  guides(colour = guide_legend(override.aes = legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="A")

shortest_line <- get_shortest_df(df_mean_2)
distance <- aggregate(shortest_line[, 5], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
distance_n <- aggregate(shortest_line[, 6], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
dist_2 <- sum(distance$dist)
dist_2_n <- sum(distance_n$dist)

p2_dist <- ggplot(shortest_line, aes(x = X_n, y = Y_n, color = node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size = 3.5) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  guides(colour = guide_legend(override.aes = legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="B")

shortest_line <- get_shortest_df(df_mean_3)
distance <- aggregate(shortest_line[, 5], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
distance_n <- aggregate(shortest_line[, 6], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
dist_3 <- sum(distance$dist)
dist_3_n <- sum(distance_n$dist)

p3_dist <- ggplot(shortest_line, aes(x = X_n, y = Y_n, color = node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size = 3.5) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  guides(colour = guide_legend(override.aes = legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="C")

shortest_line <- get_shortest_df(df_mean_4)
distance <- aggregate(shortest_line[, 5], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
distance_n <- aggregate(shortest_line[, 6], list(shortest_line$node), mean) %>%
  rename(node = Group.1) %>%
  rename(dist = x)
dist_4 <- sum(distance$dist)
dist_4_n <- sum(distance_n$dist)

p4_dist <- ggplot(shortest_line, aes(x = X_n, y = Y_n, color = node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size = 3.5) +
  geom_line(data = ref_line_data, aes(x = X_n, y = neutral_n, color = "neutral line"), linetype="dashed") +
  geom_point(aes(color = node), size = 1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=palette) +
  guides(colour = guide_legend(override.aes = legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="D")

p1_dist <- p1_dist + theme(legend.position = "top")
legend <- get_legend(p1_dist)

# Remove the legends from the plots
# t,r,b,l
p1_dist <- p1_dist + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank())
p2_dist <- p2_dist + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())
p3_dist <- p3_dist + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())
p4_dist <- p4_dist + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())

ylab <- textGrob("Change absorption capability", rot = 90, gp = gpar(fontfamily = "", size = 10, cex = 1.5))
xlab <- textGrob("Impact on performance", gp = gpar(fontfamily = "", size = 10, cex = 1.5))
p1_label <- textGrob(substitute(paste("(a) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[1])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))
p2_label <- textGrob(substitute(paste("(b) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[2])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))
p3_label <- textGrob(substitute(paste("(c) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[3])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))
p4_label <- textGrob(substitute(paste("(d) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[4])), gp = gpar(fontfamily = "LM Roman 10", size = 8, cex = 1.5))

p_dist <- grid.arrange(ylab, legend, p1_dist, p2_dist, p3_dist, p4_dist, xlab, ncol=3, nrow=4,
                          layout_matrix = rbind(c(1,2,2), c(1,3,4), c(1,5,6), c(1,7,7)),
                          widths = c(0.2, 2.7, 2.7), heights = c(0.2, 2.5, 2.5, 0.2))

# Aggregate plots


########################################
## 4. Export comparison plots
########################################
ggsave(paste0(main_wd,img_dir,"scatter_mean_comp.pdf"),
       plot = p_concept,
       dpi = 320, height = 10, width = 13)

ggsave(paste0(main_wd,img_dir,"scatter_mean_dist.pdf"),
       plot = p_dist,
       dpi = 320, height = 10, width = 13)

########################################
## 5. mean plot of single concept
########################################

df_mean <- df_mean_1

plot <- ggplot(df_mean, aes(x = Impact, y = Absorption, color = node)) +
  geom_line(data = ref_line_data, aes(x = X, y = neutral), color = "black", linetype="dashed") +
    geom_point(aes(color = node), size = 2.5) +
    geom_point(shape = 1, color = "black", size = 3, stroke = 0.1, alpha = 0.2) +
  geom_rug() +
  labs(color = "margin node\n") +
  scale_y_continuous(limits=c(a_min,a_max), name = "Change absorption capability") +
  scale_x_continuous(limits=c(i_min,i_max), name = "Impact on performance") +
  theme_pubr() +
  export_theme +
  theme(legend.position = c(0.15, 0.9),
        plot.margin=unit(c(-0.1,-0.1,0.0,0.0), "cm"))

dens1 <- ggplot(df_mean, aes(x = Impact, fill = node)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits=c(i_min,i_max), name = "") +
  labs(color = "margin node\n") +
  theme_void() +
  theme(legend.position = "none") +
  export_theme

dens2 <- ggplot(df_mean, aes(x = Absorption, fill = node)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits=c(a_min,a_max), name = "") +
  labs(color = "margin node\n") +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip() +
  export_theme

# Remove the x tick labels from the plots
# t,r,b,l
dens1 <- dens1 + theme(axis.ticks.x = element_blank(),
                       axis.text.x = element_blank(),
                       plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))
dens2 <- dens2 + theme(axis.ticks.y = element_blank(),
                       axis.text.y = element_blank(),
                       plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))

p_mean <- dens1 + plot_spacer() + plot + dens2 +
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = "collect") & theme(legend.position = "top")

########################################
## 5. CDF plot of single concept
########################################

plot <- ggplot(df_mean, aes(x = Impact, y = Absorption, color = node)) +
  geom_line(data = ref_line_data, aes(x = X, y = neutral), color = "black", linetype="dashed") +
  geom_point(aes(color = node), size = 2.5) +
  geom_point(shape = 1, color = "black", size = 3, stroke = 0.1, alpha = 0.2) +
  geom_rug() +
  labs(color = "margin node\n") +
  scale_y_continuous(limits=c(a_min,a_max), name = "Change absorption capability") +
  scale_x_continuous(limits=c(i_min,i_max), name = "Impact on performance") +
  theme_pubr() +
  export_theme +
  theme(legend.position = c(0.15, 0.9),
        plot.margin=unit(c(-0.1,-0.1,0.0,0.0), "cm"))

cdf1 <- ggplot(df_mean, aes(x = Impact, fill = node, color = node)) +
  stat_ecdf(geom = "step") +
  scale_x_continuous(limits=c(i_min,i_max), name = "") +
  ylab("Probability") +
  labs(color = "margin node\n") +
  theme_void() +
  theme(legend.position = "none") +
  export_theme

# get function values
tensor_comp <- df_mean$Impact[df_mean$node==2]
dat_ecdf <-
  data.frame(x=unique(tensor_comp),
             y=ecdf(tensor_comp)(unique(tensor_comp))*length(tensor_comp))
#rescale y to 0,1 range
dat_ecdf$y <-
  scale(dat_ecdf$y,center=min(dat_ecdf$y),scale=diff(range(dat_ecdf$y)))

# plot using computed CDF
ggplot(dat_ecdf[order(dat_ecdf$x),],aes(x,y)) +
  geom_step()

cdf2 <- ggplot(df_mean, aes(x = Absorption, fill = node, color = node)) +
  stat_ecdf(geom = "step") +
  scale_x_continuous(limits=c(a_min,a_max), name = "") +
  ylab("Probability") +
  labs(color = "margin node\n") +
  scale_y_continuous(breaks = c(0.0,0.5,1.0)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip() +
  export_theme

# Remove the x tick labels from the plots
# t,r,b,l
cdf1 <- cdf1 + theme(axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))
cdf2 <- cdf2 + theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))

p_cdf <- cdf1 + plot_spacer() + plot + cdf2 +
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = "collect") & theme(legend.position = "top")


########################################
## 5. PDF + CDF plot of single concept
########################################
pad <- plot_spacer() + theme_void()
#
# p_cdf_pdf <- grid.arrange(cdf1, dens1, plot, dens2, cdf2, pad, ncol=3, nrow=3,
#                        layout_matrix = rbind(c(1,6,6), c(2,6,6), c(3,4,5)),
#                         widths = c(4.0, 1.0, 1.0), heights = c(1.0, 1.0, 4.0))

p_cdf_pdf <- cdf1 + pad + pad + dens1 + pad + pad + plot + dens2 + cdf2 +
  plot_layout(ncol = 3, nrow = 3, widths = c(4, 1, 1), heights = c(1, 1, 4), guides = "collect") & theme(legend.position = "top")


########################################
## 6. Export plots
########################################
ggsave(paste0(main_wd,img_dir,"C1/scatter_mean.pdf"),
       plot = p_mean,
       dpi = 320, width = 8, height=8)

ggsave(paste0(main_wd,img_dir,"C1/scatter_cdf.pdf"),
       plot = p_cdf,
       dpi = 320, width = 8, height=8)

ggsave(paste0(main_wd,img_dir,"C1/scatter_cdf_pdf.pdf"),
       plot = p_cdf_pdf,
       dpi = 320, width = 8, height=8)