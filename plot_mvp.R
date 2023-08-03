## IMPORTANT NOTE: Code to visualize the scatter pairs of impact and absorption
## matrices of FEA strut problem

########################################
## 1. Headers
########################################
library(ggthemes)
library(latex2exp)
library(ggpubr)
library(patchwork)
library(tibble)
library(gridExtra)
library(grid)
library(tidyverse)
library(scales)

#HOME_WD <- "~/"
HOME_WD <- Sys.getenv("GITDIR")
## Where the working directory is
main_wd <- paste0(HOME_WD,"/mvm_paper/")

# # For fea example (polygonal)
# folder <- "strut_fea_poly"
# legend_aes <- list(linetype=c(0,0,0,2), shape=c(16,16,16,NA))
# legend_dist_aes <- list(linetype=c(2,2,2,2), shape=c(16,16,16,NA))
# node_labels <- c("$e_1$", "$e_2$", "$e_3$")
# shapes <- c(15,18,16,17)
# concept_labels <- c("1A", "1B", "1C", "1D")
# thetas <- c(0.00, 30.0, 0.00, 30.0)
# heights <- c(15.0,15.0,20.0,20.0)
# palette <- hue_pal()(3) # "#F8766D","#00BA38","#619CFF"
# # fea limits
# x0 <- array(c(-0.009630757, 1.152203562))
# x1 <- array(c(0.3329746, 4.3151886))

# For fea example (circumferential)
folder <- "strut_fea_circ"
legend_aes <- list(linetype=c(0,0,0,0,2), shape=c(16,16,16,16,NA))
legend_dist_aes <- list(linetype=c(2,2,2,2,2), shape=c(16,16,16,16,NA))
node_labels <- c("$e_1$", "$e_2$", "$e_3$", "$e_4$")
shapes <- c(3,4,8,6)
concept_labels <- c("2A", "2B", "2C", "2D")
thetas <- c(0.00, 30.0, 0.00, 30.0)
heights <- c(15.0,15.0,20.0,20.0)
palette <- c(hue_pal()(3),"#C77CFF") # "#F8766D" "#00BA38" "#619CFF" "#C77CFF"
# fea limits
x0 <- array(c(-0.009630757, 1.152203562))
x1 <- array(c(0.3329746, 4.3151886))

# # For simple example
# folder <- "strut"
# legend_dist_aes <- list(linetype=c(2,2,2), shape=c(16,16,NA))
# legend_aes <- list(linetype=c(0,0,2), shape=c(16,16,NA))
# thetas=c(0.00, 30.0, 0.00, 30.0)
# heights=c(15.0,15.0,20.0,20.0)
# palette <- hue_pal()(2) # "#F8766D","#00BFC4"
# # simple limits
# x0 <- array(c(0.00, 2.5))
# x1 <- array(c(0.075, 3.5))

## CHANGE TO MAIN WD
setwd(main_wd)
source(paste0(main_wd,"/plot_funcs.R"))
data_dir <- paste0("data/",folder,"/")
img_dir <- paste0("images/",folder,"/")
dir.create(file.path(img_dir), showWarnings=FALSE, recursive=TRUE)
colors_vector <- extend_palette(length(node_labels),length(thetas),palette)

export_theme <- export_theme + theme(plot.margin=unit(c(0,0,0,0),units="cm"),
                                     plot.tag=element_text(size=14,face="bold"),
                                     plot.tag.position=c(0.15, 0.95))

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

    dat <- data.frame(X=c(s[1],unlist(nearest_data[1])[1]),
                      Y=c(s[2],unlist(nearest_data[1])[2]),
                      X_n=c(s_n[1],unlist(nearest_data_n[1])[1]),
                      Y_n=c(s_n[2],unlist(nearest_data_n[1])[2]),
                      dist=as.numeric(nearest_data[2]),
                      dist_n=as.numeric(nearest_data_n[2]),
                      node=node_i)

    datalist[[node_i]] <- dat # add it to your list
  }

  shortest_line <- do.call(rbind, datalist) %>%
    mutate(node=as.factor(node)) # convert to categorical type

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
  mutate(node=as.factor(node)) # convert to categorical type

df_mean_2 <- read.csv(paste0(data_wd_2,"df_mean.csv"))
df_mean_2 <- df_mean_2 %>%
  mutate(node=as.factor(node)) # convert to categorical type

df_mean_3 <- read.csv(paste0(data_wd_3,"df_mean.csv"))
df_mean_3 <- df_mean_3 %>%
  mutate(node=as.factor(node)) # convert to categorical type

df_mean_4 <- read.csv(paste0(data_wd_4,"df_mean.csv"))
df_mean_4 <- df_mean_4 %>%
  mutate(node=as.factor(node)) # convert to categorical type


i_min <- min(c(df_mean_1$Impact,df_mean_2$Impact,df_mean_3$Impact,df_mean_4$Impact))
a_min <- min(c(df_mean_1$Absorption,df_mean_2$Absorption,df_mean_3$Absorption,df_mean_4$Absorption))
i_max <- max(c(df_mean_1$Impact,df_mean_2$Impact,df_mean_3$Impact,df_mean_4$Impact))
a_max <- max(c(df_mean_1$Absorption,df_mean_2$Absorption,df_mean_3$Absorption,df_mean_4$Absorption))

# # automatic limits
# x0 <- array(c(i_min,a_min))
# x1 <- array(c(i_max,a_max))

x0_n <- array(c(0,0))
x1_n <- array(c(1,1))

#######################################
# scale if necessary
df_mean_1 <- df_mean_1 %>%
  mutate(Impact, Impact_n=(Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n=(Absorption - x0[2]) / (x1[2] - x0[2])) %>%
  mutate(concept=concept_labels[1])
df_mean_2 <- df_mean_2 %>%
  mutate(Impact, Impact_n=(Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n=(Absorption - x0[2]) / (x1[2] - x0[2])) %>%
  mutate(concept=concept_labels[2])
df_mean_3 <- df_mean_3 %>%
  mutate(Impact, Impact_n=(Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n=(Absorption - x0[2]) / (x1[2] - x0[2])) %>%
  mutate(concept=concept_labels[3])
df_mean_4 <- df_mean_4 %>%
  mutate(Impact, Impact_n=(Impact - x0[1]) / (x1[1] - x0[1])) %>%
  mutate(Absorption, Absorption_n=(Absorption - x0[2]) / (x1[2] - x0[2])) %>%
  mutate(concept=concept_labels[4])

ref_line_data <- data.frame(X=c(i_min, i_max),
                            neutral=c(a_min, a_max),
                            X_n=c(0, 1),
                            neutral_n=c(0, 1))

########################################
## 2. Scatter plot of three concepts
########################################
p1 <- ggplot(df_mean_1, aes(x=Impact_n, y=Absorption_n, color=node)) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  # geom_point(shape=1, color="black", size=2, stroke=0.1, alpha=0.1) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  labs(color="margin node\n") +
  guides(colour=guide_legend(override.aes=legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(a)")

p2 <- ggplot(df_mean_2, aes(x=Impact_n, y=Absorption_n, color=node)) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  guides(colour=guide_legend(override.aes=legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(b)")

p3 <- ggplot(df_mean_3, aes(x=Impact_n, y=Absorption_n, color=node)) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  # geom_point(shape=1, color="black", size=2, stroke=0.1, alpha=0.1) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  guides(colour=guide_legend(override.aes=legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(c)")

p4 <- ggplot(df_mean_4, aes(x=Impact_n, y=Absorption_n, color=node)) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  # geom_point(shape=1, color="black", size=2, stroke=0.1, alpha=0.1) +
  geom_rug() +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  guides(colour=guide_legend(override.aes=legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(d)")

p1 <- p1 + theme(legend.position="top")
legend <- get_legend(p1)

# Remove the legends from the plots
# t,r,b,l
p1 <- p1 + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank())
p2 <- p2 + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())
p3 <- p3 + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())
p4 <- p4 + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())

ylab <- textGrob("Change absorption capability", rot=90, gp=gpar(fontfamily="", size=10, cex=1.5))
xlab <- textGrob("Impact on performance", gp=gpar(fontfamily="", size=10, cex=1.5))
p1_label <- textGrob(substitute(paste("(a) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[1])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))
p2_label <- textGrob(substitute(paste("(b) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[2])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))
p3_label <- textGrob(substitute(paste("(c) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[3])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))
p4_label <- textGrob(substitute(paste("(d) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[4])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))

p_concept <- grid.arrange(ylab, legend, p1, p2, p3, p4, xlab, ncol=3, nrow=4,
             layout_matrix=rbind(c(1,2,2), c(1,3,4), c(1,5,6), c(1,7,7)),
             widths=c(0.2, 2.7, 2.7), heights=c(0.2, 2.5, 2.5, 0.2))

# Aggregate plots
df_mean_agg <- do.call("rbind", list(df_mean_1,df_mean_2,df_mean_3,df_mean_4)) %>%
  mutate(concept=factor(concept,concept_labels))

# set up shapes
names(shapes) <- concept_labels
sizes <- c(2.0,2.0,2.0,2.0)
names(sizes) <- concept_labels
strokes <- seq(from=0.8,to=0.2,by=-0.2)
names(strokes) <- concept_labels
alphas <- seq(from=0.8,to=0.2,by=-0.2)
names(alphas) <- concept_labels

p_agg <- ggplot(df_mean_agg, aes(x=Impact_n, y=Absorption_n, color=node)) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n), color="black", linetype="dashed") +
  # geom_point(data=ref_line_data, shape=1, color="black", size=2, stroke=0.1, alpha=0.1) +
  geom_point(aes(shape=concept,size=concept,alpha=concept)) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c(palette),
                      labels = unname(TeX(node_labels)),
                      guide = guide_legend(override.aes=aes(size=4.0))) +
  scale_shape_manual(values=shapes) +
  scale_size_manual(values=sizes) +
  scale_alpha_manual(values=alphas,
                     guide = guide_legend(override.aes=aes(alpha=NA, size=4.0))) +
  labs(color="margin node\n") +
  xlab("Impact on performance") +
  ylab("Change absorption capability") +
  # guides(color=guide_legend(override.aes=legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(a)")
p_agg

if (folder == "strut_fea_circ") {
  df_mean_agg_circ <- df_mean_agg
} else {
  df_mean_agg_poly <- df_mean_agg
}

########################################
## 3. Export distance plots
########################################
shortest_line_1 <- get_shortest_df(df_mean_1) %>%
  mutate(concept=concept_labels[1])
distance <- aggregate(shortest_line_1[, 5], list(shortest_line_1$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
distance_n <- aggregate(shortest_line_1[, 6], list(shortest_line_1$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
dist_1 <- sum(distance$dist)
dist_1_n <- sum(distance_n$dist)

p1_dist <- ggplot(shortest_line_1, aes(x=X_n, y=Y_n, color=node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size=3.5) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  labs(color="margin node\n") +
  guides(colour=guide_legend(override.aes=legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(a)")

shortest_line_2 <- get_shortest_df(df_mean_2) %>%
  mutate(concept=concept_labels[2])
distance <- aggregate(shortest_line_2[, 5], list(shortest_line_2$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
distance_n <- aggregate(shortest_line_2[, 6], list(shortest_line_2$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
dist_2 <- sum(distance$dist)
dist_2_n <- sum(distance_n$dist)

p2_dist <- ggplot(shortest_line_2, aes(x=X_n, y=Y_n, color=node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size=3.5) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  guides(colour=guide_legend(override.aes=legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(b)")

shortest_line_3 <- get_shortest_df(df_mean_3) %>%
  mutate(concept=concept_labels[3])
distance <- aggregate(shortest_line_3[, 5], list(shortest_line_3$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
distance_n <- aggregate(shortest_line_3[, 6], list(shortest_line_3$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
dist_3 <- sum(distance$dist)
dist_3_n <- sum(distance_n$dist)

p3_dist <- ggplot(shortest_line_3, aes(x=X_n, y=Y_n, color=node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size=3.5) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  guides(colour=guide_legend(override.aes=legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(c)")

shortest_line_4 <- get_shortest_df(df_mean_4) %>%
  mutate(concept=concept_labels[4])
distance <- aggregate(shortest_line_4[, 5], list(shortest_line_4$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
distance_n <- aggregate(shortest_line_4[, 6], list(shortest_line_4$node), mean) %>%
  rename(node=Group.1) %>%
  rename(dist=x)
dist_4 <- sum(distance$dist)
dist_4_n <- sum(distance_n$dist)

p4_dist <- ggplot(shortest_line_4, aes(x=X_n, y=Y_n, color=node)) +
  geom_line(linetype="dashed",size=1) +
  geom_point(size=3.5) +
  geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n, color="neutral line"), linetype="dashed") +
  geom_point(aes(color=node), size=1.5) +
  scale_y_continuous(name="Absorption", limits=c(0,1)) +
  scale_x_continuous(name="Impact", limits=c(0,1)) +
  scale_colour_manual(values=c(palette,"black")) +
  guides(colour=guide_legend(override.aes=legend_dist_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(d)")

p1_dist <- p1_dist + theme(legend.position="top")
legend <- get_legend(p1_dist)

# Remove the legends from the plots
# t,r,b,l
p1_dist <- p1_dist + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank())
p2_dist <- p2_dist + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())
p3_dist <- p3_dist + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())
p4_dist <- p4_dist + theme(legend.position="none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank())

ylab <- textGrob("Change absorption capability", rot=90, gp=gpar(fontfamily="", size=10, cex=1.5))
xlab <- textGrob("Impact on performance", gp=gpar(fontfamily="", size=10, cex=1.5))
p1_label <- textGrob(substitute(paste("(a) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[1])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))
p2_label <- textGrob(substitute(paste("(b) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[2])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))
p3_label <- textGrob(substitute(paste("(c) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[3])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))
p4_label <- textGrob(substitute(paste("(d) ",theta==th*degree, h==height," mm"), list(th=thetas[1], height=heights[4])), gp=gpar(fontfamily="LM Roman 10", size=8, cex=1.5))

p_dist <- grid.arrange(ylab, legend, p1_dist, p2_dist, p3_dist, p4_dist, xlab, ncol=3, nrow=4,
                          layout_matrix=rbind(c(1,2,2), c(1,3,4), c(1,5,6), c(1,7,7)),
                          widths=c(0.2, 2.7, 2.7), heights=c(0.2, 2.5, 2.5, 0.2))

# Aggregate plots
df_line_agg <- do.call("rbind", list(shortest_line_1,shortest_line_2,shortest_line_3,shortest_line_4)) %>%
  mutate(concept=as.factor(concept))

p_dist_agg <- ggplot(data=df_line_agg, aes(x=X_n, y=Y_n, color=node, group=interaction(node,concept))) +
  geom_line(color="grey",linetype="solid",size=0.5) +
  # geom_point(size=3.5) +
  geom_point(data=df_line_agg %>% filter(row_number() %% 2 == 1), ## Select odd rows
             aes(color=node,shape=concept), size=2.5) +
  geom_line(data=ref_line_data %>% mutate(concept=NA) %>% mutate(node=NA), aes(x=X_n, y=neutral_n), color="black", linetype="dashed") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c(palette)) +
  scale_shape_manual(values=shapes) +
  labs(color="margin node\n") +
  xlab("Impact on performance") +
  ylab("Change absorption capability") +
  # guides(color=guide_legend(override.aes=legend_aes)) +
  theme_pubr() +
  export_theme +
  labs(tag="(b)")
p_dist_agg

if (folder == "strut_fea_circ") {
  df_line_agg_circ <- df_line_agg
} else {
  df_line_agg_poly <- df_line_agg
}

# Combine dist and scatter plots together
p_agg_export <- p_agg + theme(legend.position="top")
legend <- get_legend(p_agg_export)

# Remove the legends from the plots
# t,r,b,l
p_agg_export <- p_agg_export + theme(legend.position="none",
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank())
p_dist_agg_export <- p_dist_agg + theme(legend.position="none",
                                        axis.title.x=element_blank(),
                                        axis.title.y=element_blank(),
                                        axis.text.y=element_blank())

ylab <- textGrob("Change absorption capability", rot=90, gp=gpar(fontfamily="", size=10, cex=1.5))
xlab <- textGrob("Impact on performance", gp=gpar(fontfamily="", size=10, cex=1.5))
blank <- grid.rect(gp=gpar(col="white",alpha=0.0))

p_agg_combined <- grid.arrange(ylab, legend, p_agg_export, p_dist_agg_export, xlab, blank, ncol=5, nrow=4,
                               layout_matrix=rbind(c(1,6,6,6,6), c(1,2,2,2,6), c(1,3,6,4,6), c(1,5,5,5,6)),
                               widths=c(0.2, 2.7, 1.0, 2.7, 0.1), heights=c(0.01, 0.2, 2.5, 0.2),
                               vp=viewport(width=1.0, height=0.9))

########################################
## 4. Export comparison plots
########################################
ggsave(paste0(main_wd,img_dir,"scatter_mean_comp.pdf"),
       plot=p_concept,
       dpi=320, height=10, width=13)

ggsave(paste0(main_wd,img_dir,"scatter_mean_dist.pdf"),
       plot=p_dist,
       dpi=320, height=10, width=13)

ggsave(paste0(main_wd,img_dir,"scatter_all.pdf"),
       plot=p_agg_combined,
       dpi=320, height=6.5, width=13)

# report reliability
reliability <- list(A=nrow(df_mean_1)/3000,B=nrow(df_mean_2)/3000,
                    C=nrow(df_mean_3)/3000,D=nrow(df_mean_4)/3000)

# report value
value <- list(A=sum(shortest_line_1$dist),B=sum(shortest_line_2$dist),
              C=sum(shortest_line_3$dist),D=sum(shortest_line_4$dist))
signif(unlist(value),3)

# node values of each concept
signif(shortest_line_1$dist_n[seq(1,nrow(shortest_line_1),2)],3)
signif(shortest_line_2$dist_n[seq(1,nrow(shortest_line_1),2)],3)
signif(shortest_line_3$dist_n[seq(1,nrow(shortest_line_1),2)],3)
signif(shortest_line_4$dist_n[seq(1,nrow(shortest_line_1),2)],3)

########################################
## 5. Export combined plots for poly and circ
########################################
  
if (exists("df_mean_agg_poly") & exists("df_mean_agg_circ")) {
  
  df_mean_comb <- do.call("rbind", list(df_mean_agg_poly,df_mean_agg_circ)) %>%
    mutate(concept=as.factor(concept))
  
  # set up shapes
  shapes <- c(15,18,16,17,3,4,8,6)
  names(shapes) <- unique(df_mean_comb$concept)
  sizes <- c(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0)
  names(sizes) <- unique(df_mean_comb$concept)
  strokes <- c(seq(from=0.8,to=0.2,by=-0.2),seq(from=0.8,to=0.2,by=-0.2))
  names(strokes) <- unique(df_mean_comb$concept)
  alphas <- c(seq(from=0.8,to=0.2,by=-0.2),seq(from=0.8,to=0.2,by=-0.2))
  names(alphas) <- unique(df_mean_comb$concept)
  
  p_comb <- ggplot(df_mean_comb, aes(x=Impact_n, y=Absorption_n, color=node)) +
    geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n), color="black", linetype="dashed") +
    # geom_point(data=ref_line_data, shape=1, color="black", size=2, stroke=0.1, alpha=0.1) +
    geom_point(aes(shape=concept,size=concept,alpha=concept)) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(limits=c(0,1)) +
    scale_colour_manual(values = c(hue_pal()(3),"#C77CFF"), # "#F8766D" "#00BA38" "#619CFF" "#C77CFF"
                        labels = unname(TeX(node_labels)),
                        guide = guide_legend(override.aes=aes(size=4.0))) +
    scale_shape_manual(values=shapes) +
    scale_size_manual(values=sizes) +
    scale_alpha_manual(values=alphas,
                       guide = guide_legend(override.aes=aes(alpha=NA, size=4.0))) +
    guides(shape=guide_legend(nrow=1,byrow=TRUE)) +
    labs(color="margin node\n") +
    xlab("Impact on performance") +
    ylab("Change absorption capability") +
    # guides(color=guide_legend(override.aes=legend_aes)) +
    theme_pubr() +
    export_theme +
    labs(tag="(a)")
  p_comb
  
}

presentation <- TRUE
if (exists("df_mean_agg_poly") & exists("df_mean_agg_circ") & presentation) {
  
  chosen <- c('')
  # chosen <- c('1A')
  # chosen <- c('1A','1B')
  # chosen <- c('1A','1B','2A')
  # chosen <- c('1A','1B','1C','1D','2A','2B','2C','2D')

  df_mean_comb_pres <- do.call("rbind", list(df_mean_agg_poly,df_mean_agg_circ)) %>%
    mutate(concept=as.factor(concept)) %>% filter(concept %in% chosen)
  
  # set up shapes
  shapes <- c(15,18,16,17,3,4,8,6)
  names(shapes) <- unique(df_mean_comb$concept)
  sizes <- c(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0)
  names(sizes) <- unique(df_mean_comb$concept)
  strokes <- c(seq(from=0.8,to=0.2,by=-0.2),seq(from=0.8,to=0.2,by=-0.2))
  names(strokes) <- unique(df_mean_comb$concept)
  alphas <- c(seq(from=0.8,to=0.2,by=-0.2),seq(from=0.8,to=0.2,by=-0.2))
  names(alphas) <- unique(df_mean_comb$concept)
  
  node_labels_pres <- c("$w_{vane}$", "material", "$n_{struts}$", "$t_{shroud}$")
  
  p_pres <- ggplot(df_mean_comb_pres, aes(x=Impact_n, y=Absorption_n, color=node)) +
    geom_line(data=ref_line_data, aes(x=X_n, y=neutral_n), color="black", linetype="dashed") +
    # geom_point(data=ref_line_data, shape=1, color="black", size=2, stroke=0.1, alpha=0.1) +
    geom_point(aes(shape=concept,size=concept,alpha=concept)) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(limits=c(0,1)) +
    scale_colour_manual(values = c(hue_pal()(3),"#C77CFF"), # "#F8766D" "#00BA38" "#619CFF" "#C77CFF"
                        labels = unname(TeX(node_labels_pres)),
                        guide = guide_legend(override.aes=aes(size=4.0))) +
    scale_shape_manual(values=shapes) +
    scale_size_manual(values=sizes) +
    scale_alpha_manual(values=alphas,
                       guide = guide_legend(override.aes=aes(alpha=NA, size=4.0))) +
    guides(shape=guide_legend(nrow=1,byrow=TRUE)) +
    labs(color="margin node\n") +
    xlab("Impact on performance") +
    ylab("Change absorption capability") +
    # guides(color=guide_legend(override.aes=legend_aes)) +
    theme_pubr() +
    theme_tufte() +
    theme(
      axis.text.x=element_text(size=18,family=""),
      axis.text.y=element_text(size=18,family=""),
      axis.title.x=element_text(size=24,family="",vjust=-1),
      axis.title.y=element_text(size=24,family=""),
      
      ## Axis lines
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      
      ## Title
      plot.title = element_text(family="",size=10,face="bold",hjust=0.5),
      plot.tag = element_text(family="",size=10,face="bold"),
      
      ## Legends
      legend.title=element_text(size=18,family="",face="italic"),
      legend.text=element_text(size=15,family=""),
      legend.key.size= unit(0.5, "cm"),
      legend.margin = margin(0,0,0,0, "cm"),
      ## Strips for facet_wrap
      strip.text=element_text(size=8,family="",face="bold"),
      #strip.background=element_rect(fill="#f0f0f0")
      strip.background=element_blank(),
    ) +
    theme(legend.position="top", legend.box ="vertical")
  p_pres
  
  ggsave(paste0(main_wd,"images/","scatter_all_pres_",paste(chosen,collapse="_"),".pdf"),
         plot=p_pres,
         dpi=320, height=7.5, width=7.5)
  
}

if (exists("df_line_agg_poly") & exists("df_line_agg_circ")) {
  
  df_line_comb <- do.call("rbind", list(df_line_agg_poly,df_line_agg_circ)) %>%
    mutate(concept=as.factor(concept))
  
  
  p_dist_comb <- ggplot(data=df_line_comb, aes(x=X_n, y=Y_n, color=node, group=interaction(node,concept))) +
    geom_line(color="grey",linetype="solid",size=0.5) +
    # geom_point(size=3.5) +
    geom_point(data=df_line_comb %>% filter(row_number() %% 2 == 1), ## Select odd rows
               aes(color=node,shape=concept), size=2.5) +
    geom_line(data=ref_line_data %>% mutate(concept=NA) %>% mutate(node=NA), aes(x=X_n, y=neutral_n), color="black", linetype="dashed") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(limits=c(0,1)) +
    scale_colour_manual(values=c(hue_pal()(3),"#C77CFF")) + # "#F8766D" "#00BA38" "#619CFF" "#C77CFF"
    scale_shape_manual(values=shapes) +
    labs(color="margin node\n") +
    xlab("Impact on performance") +
    ylab("Change absorption capability") +
    # guides(color=guide_legend(override.aes=legend_aes)) +
    theme_pubr() +
    export_theme +
    labs(tag="(b)")
  p_dist_comb
  
}

if (exists("p_comb") & exists("p_dist_comb")) {
  
  # Combine dist and scatter plots together
  p_comb_export <- p_comb + theme(legend.position="top", legend.box ="horizontal")
  p_comb_export
  legend <- get_legend(p_comb_export)
  
  # Remove the legends from the plots
  # t,r,b,l
  p_comb_export <- p_comb_export + theme(legend.position="none",
                                         axis.title.x=element_blank(),
                                         axis.title.y=element_blank())
  p_dist_comb_export <- p_dist_comb + theme(legend.position="none",
                                            axis.title.x=element_blank(),
                                            axis.title.y=element_blank(),
                                            axis.text.y=element_blank())
  
  ylab <- textGrob("Change absorption capability", rot=90, gp=gpar(fontfamily="", size=10, cex=1.5))
  xlab <- textGrob("Impact on performance", gp=gpar(fontfamily="", size=10, cex=1.5))
  blank <- grid.rect(gp=gpar(col="white",alpha=0.0))
  
  p_comb_combined <- grid.arrange(ylab, legend, p_comb_export, p_dist_comb_export, xlab, blank, ncol=5, nrow=4,
                                 layout_matrix=rbind(c(1,6,6,6,6), c(1,2,2,2,6), c(1,3,6,4,6), c(1,5,5,5,6)),
                                 widths=c(0.2, 2.7, 1.0, 2.7, 0.1), heights=c(0.01, 0.2, 2.5, 0.2),
                                 vp=viewport(width=1.0, height=0.9))
  
  ggsave(paste0(main_wd,"images/","scatter_all.pdf"),
         plot=p_comb_combined,
         dpi=320, height=6.5, width=13)
  
}

########################################
## 6. mean plot of single concept
########################################
i_min <- min(c(df_mean_1$Impact))
a_min <- min(c(df_mean_1$Absorption))
i_max <- max(c(df_mean_1$Impact))
a_max <- max(c(df_mean_1$Absorption))

ref_line_int <- data.frame(approx(ref_line_data$X, ref_line_data$neutral, n=100)) %>%
  filter(x > i_min) %>%
  filter(y > a_min)

df_mean <- df_mean_1

plot <- ggplot(df_mean, aes(x=Impact, y=Absorption, color=node)) +
  geom_line(data=ref_line_int, aes(x=x, y=y), color="black", linetype="dashed") +
  geom_point(aes(color=node), size=2.5) +
  geom_point(shape=1, color="black", size=3, stroke=0.1, alpha=0.2) +
  geom_rug() +
  labs(color="margin node\n") +
  scale_y_continuous(limits=c(a_min,a_max), name="Change absorption capability") +
  scale_x_continuous(limits=c(i_min,i_max), name="Impact on performance") +
  scale_colour_discrete(guide = "none") +
  theme_pubr() +
  export_theme +
  theme(legend.position=c(0.15, 0.9),
        plot.margin=unit(c(-0.1,-0.1,0.0,0.0), "cm"))
plot

dens1 <- ggplot(df_mean, aes(x=Impact, fill=node)) +
  geom_density(alpha=0.4) +
  scale_x_continuous(limits=c(i_min,i_max), name="") +
  labs(fill="margin node\n") +
  scale_fill_discrete(labels = unname(TeX(node_labels))) +
  theme_void() +
  theme(legend.position="none") +
  export_theme

dens2 <- ggplot(df_mean, aes(x=Absorption, fill=node)) +
  geom_density(alpha=0.4) +
  scale_x_continuous(limits=c(a_min,a_max), name="") +
  labs(color="margin node\n") +
  scale_fill_discrete(guide = "none") +
  theme_void() +
  theme(legend.position="none") +
  coord_flip() +
  export_theme

# Remove the x tick labels from the plots
# t,r,b,l
dens1 <- dens1 + theme(axis.ticks.x=element_blank(),
                       axis.text.x=element_blank(),
                       plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))
dens2 <- dens2 + theme(axis.ticks.y=element_blank(),
                       axis.text.y=element_blank(),
                       plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))

p_mean <- dens1 + plot_spacer() + plot + dens2 +
  plot_layout(ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4), guides="collect") & theme(legend.position="top")

########################################
## 7. CDF plot of single concept
########################################

plot <- ggplot(df_mean, aes(x=Impact, y=Absorption, color=node)) +
  geom_line(data=ref_line_int, aes(x=x, y=y), color="black", linetype="dashed") +
  geom_point(aes(color=node), size=2.5) +
  geom_point(shape=1, color="black", size=3, stroke=0.1, alpha=0.2) +
  geom_rug() +
  labs(color="margin node\n") +
  scale_y_continuous(limits=c(a_min,a_max), name="Change absorption capability") +
  scale_x_continuous(limits=c(i_min,i_max), name="Impact on performance") +
  scale_colour_discrete(labels = unname(TeX(node_labels))) +
  theme_pubr() +
  export_theme +
  theme(legend.position=c(0.15, 0.9),
        plot.margin=unit(c(-0.1,-0.1,0.0,0.0), "cm"))

cdf1 <- ggplot(df_mean, aes(x=Impact, fill=node, color=node)) +
  stat_ecdf(geom="step") +
  scale_x_continuous(limits=c(i_min,i_max), name="") +
  ylab("probability") +
  labs(color="margin node\n") +
  scale_colour_discrete(guide = "none") +
  theme_void() +
  theme(legend.position="none") +
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

cdf2 <- ggplot(df_mean, aes(x=Absorption, fill=node, color=node)) +
  stat_ecdf(geom="step") +
  scale_x_continuous(limits=c(a_min,a_max), name="") +
  ylab("probability") +
  labs(color="margin node\n") +
  scale_colour_discrete(guide = "none") +
  scale_y_continuous(breaks=c(0.0,0.5,1.0)) +
  theme_void() +
  theme(legend.position="none") +
  coord_flip() +
  export_theme

# Remove the x tick labels from the plots
# t,r,b,l
cdf1 <- cdf1 + theme(axis.ticks.x=element_blank(),
                     axis.text.x=element_blank(),
                     plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))
cdf2 <- cdf2 + theme(axis.ticks.y=element_blank(),
                     axis.text.y=element_blank(),
                     plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))

p_cdf <- cdf1 + plot_spacer() + plot + cdf2 +
  plot_layout(ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4), guides="collect") & theme(legend.position="top")

########################################
## 8. CDF plot of single node
########################################
df_total <- read.csv(paste0(main_wd,data_dir,"C1/df_total.csv")) %>%
  mutate(across(6:8, as.factor))
df_total$decision <- as.factor(with(df_total,paste(D1,D2,D3,sep=",")))

df_maxmin <- data.frame(c(i_min,i_max),c(a_min,a_max))
colnames(df_maxmin) <- c("I3","A3")

my_cdf <- function(data, mapping, ...) {
  col = as.character(mapping$x)[2]
  fcol = as.character(mapping$colour)[2]
  
  dfs = process_cdf(data,df_maxmin,col,fcol)
  marks = dfs[[1]]
  df_cdf = dfs[[2]]
  
  # plot everything
  plot <- ggplot(data = df_cdf, mapping=aes(x=x_cdf,y=y_cdf,color=decision,fill=decision)) +
    geom_line() +
    geom_segment(data=marks, aes(x=interval_vector, xend=interval_vector,y=0.0, yend=interval_limits,color=NA)) +
    geom_area(position="identity",alpha=0.5) +
    scale_x_continuous(limits=c(as.numeric(df_maxmin[1,col]),
                                as.numeric(df_maxmin[2,col])))
  return(list(plot,marks))
}

labels <- levels(df_total$decision)
labels_plot <- process_labels(labels)

col <- "I3"
annotate_levels <- c(3)
outs <- my_cdf(df_total, aes_string(x=col,colour="decision",fill="decision"))
marks <- outs[[2]]
cdf1_E <- outs[[1]] +
  geom_segment(data=marks[annotate_levels,], 
               aes(x=-Inf, xend=interval_vector,
                   y=interval_limits, yend=interval_limits,
                   color=NA), linetype="longdash",size=0.5) +
  annotate(geom="text", x=0.015+as.numeric(df_maxmin[1,col]), 
           y=marks[annotate_levels,"interval_limits"]+0.15, 
           label=round(marks[annotate_levels,"interval_limits"], digits=3),
           color="black", size=5) +
  ylab("probability") +
  xlab("") +
  scale_fill_discrete(labels = labels_plot) +
  labs(fill = unname(TeX("$w$, material, $n_{struts}$"))) +
  scale_colour_discrete(guide = "none") +
  theme_void() +
  theme(legend.position="none") +
  guides(fill = guide_legend(nrow = 2)) +
  export_theme
cdf1_E

col <- "A3"
annotate_levels <- c(4)
outs <- my_cdf(df_total, aes_string(x=col,colour="decision",fill="decision"))
marks <- outs[[2]]
cdf2_E <- outs[[1]] +
  geom_segment(data=marks[annotate_levels,], 
               aes(x=-Inf, xend=interval_vector,
                   y=interval_limits, yend=interval_limits,
                   color=NA), linetype="longdash",size=0.5) +
  annotate(geom="text", x=0.015+as.numeric(df_maxmin[1,col]), 
           y=marks[annotate_levels,"interval_limits"]-0.3, 
           label=round(marks[annotate_levels,"interval_limits"], digits=3),
           color="black", size=5) +
  ylab("probability") +
  xlab("") +
  scale_fill_discrete(labels = labels_plot) +
  labs(fill = unname(TeX("$w$, material, $n_{struts}$"))) +
  scale_colour_discrete(guide = "none") +
  scale_y_continuous(breaks=c(0.0,0.5,1.0)) +
  theme_void() +
  theme(legend.position="none") +
  coord_flip() +
  export_theme
cdf2_E

plot_E <- ggplot(df_total, aes(x=I3, y=A3, color=decision)) +
  geom_point(aes(color=decision), size=2.5) +
  # geom_point(shape=1, color="black", size=3, stroke=0.1, alpha=0.2) +
  geom_line(data=ref_line_int, aes(x=x, y=y), color="black", linetype="dashed") +
  labs(color="margin node\n") +
  scale_y_continuous(limits=c(a_min,a_max), name="Change absorption capability") +
  scale_x_continuous(limits=c(i_min,i_max), name="Impact on performance") +
  theme_pubr() +
  export_theme +
  theme(legend.position=c(0.15, 0.9),
        plot.margin=unit(c(-0.1,-0.1,0.0,0.0), "cm"))
plot_E

# Remove the x tick labels from the plots
# t,r,b,l
cdf1_E <- cdf1_E + theme(axis.ticks.x=element_blank(),
                         axis.text.x=element_blank(),
                         plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm"))
cdf2_E <- cdf2_E + theme(axis.ticks.y=element_blank(),
                         axis.text.y=element_blank(),
                         plot.margin=unit(c(0.0,0.0,0.0,0.0), "cm")) +
  scale_fill_discrete(guide = "none")
plot_E <- plot_E +
  scale_color_discrete(guide = "none")

p_cdf_E <- cdf1_E + plot_spacer() + plot_E + cdf2_E +
  plot_layout(ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4), guides="collect") & theme(legend.position="top")
p_cdf_E

########################################
## 9. PDF + CDF plot of single concept
########################################
pad <- plot_spacer() + theme_void()
#
# p_cdf_pdf <- grid.arrange(cdf1, dens1, plot, dens2, cdf2, pad, ncol=3, nrow=3,
#                        layout_matrix=rbind(c(1,6,6), c(2,6,6), c(3,4,5)),
#                         widths=c(4.0, 1.0, 1.0), heights=c(1.0, 1.0, 4.0))

p_cdf_pdf <- cdf1 + pad + pad + dens1 + pad + pad + plot + dens2 + cdf2 +
  plot_layout(ncol=3, nrow=3, widths=c(4, 1, 1), heights=c(1, 1, 4), guides="collect") & theme(legend.position="top")


########################################
## 10. Export plots
########################################
dir.create(file.path(paste0(img_dir,"C1/")), showWarnings=FALSE, recursive=TRUE)

ggsave(paste0(main_wd,img_dir,"C1/scatter_mean.pdf"),
       plot=p_mean,
       dpi=320, width=8, height=8)

ggsave(paste0(main_wd,img_dir,"C1/scatter_cdf.pdf"),
       plot=p_cdf,
       dpi=320, width=8, height=8)

ggsave(paste0(main_wd,img_dir,"C1/scatter_cdf_E3.pdf"),
       plot=p_cdf_E,
       dpi=320, width=8, height=8)

ggsave(paste0(main_wd,img_dir,"C1/scatter_cdf_pdf.pdf"),
       plot=p_cdf_pdf,
       dpi=320, width=8, height=8)