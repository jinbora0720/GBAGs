rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(scico)
library(bags)
source("~/GBAGs/cov/Ctilde_ST.R")

# path
path <- "~/GBAGs/"

# plot extensions
extension <- c(".pdf", ".eps", ".png")

########
# Data #
########
# specify number of grid on each axis
ngrid3_y <- 3
ngrid3_x <- 31
n_time3 <- 31
xgrid3 <- seq(0, 1, length = ngrid3_x)
ygrid3 <- seq(0, 1, length = ngrid3_y)
tgrid3 <- seq(0, 1, length = n_time3)
xytgrid3 <- expand.grid(easting = xgrid3, northing = ygrid3, time = tgrid3) %>%
  arrange(time, easting, northing)
n3 <- nrow(xytgrid3)
coords3 <- xytgrid3

# wind directions
n_northing <- 3
n_easting <- 31
nd3 <- floor(log10(max(n_northing, n_easting, n_time3))) + 1
format3 <- paste0("%0", nd3, ".0f,","%0", nd3, ".0f,", "%0", nd3, ".0f")
coords_ptt3 <- coords3
coords_ptt3$row <- (n_northing+1) -
  as.numeric(cut_interval(coords_ptt3$northing,
                          n = n_northing,
                          labels = 1:n_northing))                               # 1 from the top
coords_ptt3$col <- as.numeric(cut_interval(coords_ptt3$easting,
                                           n = n_easting,
                                           labels = 1:n_easting))
coords_ptt3$time_d <- as.numeric(as.character(coords_ptt3$time*(n_time3-1) + 1))
coords_ptt3$partition <- sprintf(format3,
                                 coords_ptt3$row,
                                 coords_ptt3$col,
                                 coords_ptt3$time_d)

#######
# Cov #
#######
# true parameter values
a3 <- 0 # the smaller the longer the cov lasts
c <- .8 # the smaller the larger the cov
kappa <- 0
sig_sq <- 1

# wind directions
ptts3 <- sort(unique(coords_ptt3$partition))
z_true13 <- matrix("W", nrow = length(ptts3))
rownames(z_true13) <- ptts3

# Ctilde given z
Cz_W3 <- Ctilde_ST(coords_ptt = coords_ptt3, z = z_true13,
                   params = list(sig_sq = sig_sq, a = a3, c = c, kappa = kappa))

#############
# Subfigure #
#############
i_ref <- which(coords_ptt3[,"partition"] == "02,16,16")
t_diff <- n_time3 - 16
Widx_ref <- i_ref + (-t_diff:t_diff)*(ngrid3_x+1)*ngrid3_y
Widx_ref <- Widx_ref[Widx_ref <= n3 & Widx_ref > 0]
Eidx_ref <- i_ref + (-t_diff:t_diff)*(ngrid3_x-1)*ngrid3_y
Eidx_ref <- Eidx_ref[Eidx_ref <= n3 & Eidx_ref > 0]

coords_ptt3[Widx_ref,"partition"]
coords_ptt3[Eidx_ref,"partition"]

cov_df <- data.frame(dist = c(tgrid3, tgrid3),
                     dir = c(rep("W", each = length(Widx_ref)),
                             rep("E", each = length(Eidx_ref))),
                     cov = c(Cz_W3[i_ref, Widx_ref], Cz_W3[i_ref, Eidx_ref]))

gg <- cov_df %>%
  ggplot() +
  geom_line(aes(dist, cov, col = dir, linetype = dir), linewidth = 1) +
  scale_color_manual(values = c("E" = "#0072B2", "W" = "#D55E00"),
                     labels = c("opposite direction", "along wind direction")) +
  scale_linetype_manual(values = c("E" = "dashed", "W" = "solid"),
                        labels = c("opposite direction", "along wind direction")) +
  geom_vline(xintercept = coords_ptt3[i_ref,"time"]) + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11), 
        legend.position = "bottom") +
  labs(x = "Time", y = "Covariance", #title = "Winds from W (left to right)",
       title = "G-BAG nonstationary covariance: directional",
       color = "",
       linetype = "") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(0, coords_ptt3[i_ref, "time"], 1),
                     labels = c(0, 0.5, 1))

# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/cov_time_newparent", ext),
#          width = 5, height = 3.2)
# }
