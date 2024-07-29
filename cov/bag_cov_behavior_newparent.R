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
ngrid <- 30
n_time <- 4
xygrid <- expand.grid(easting = seq(-1, 1, length = ngrid),
                      northing = seq(-1, 1, length = ngrid))
knots <- data.frame(x = c(0, 0, 0, 1, -1, 1/sqrt(2),
                          1/sqrt(2), -1/sqrt(2), -1/sqrt(2)),
                    y = c(0, 1, -1, 0, 0, 1/sqrt(2),
                          -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)))
nn_data <- FNN::get.knnx(0.65*knots, xygrid, k=1)
xygrid1 <- xygrid %>% dplyr::mutate(partition = nn_data$nn.index)
# ggplot(xygrid1, aes(x = easting, y = northing, color = factor(partition))) +
#   geom_point(size = 4) +
#   geom_point(data = knots, aes(x = x, y = y), col = "red",size = 10)

xgrid <- seq(0, 1, length = ngrid)
tgrid <- seq(0, 1, length = n_time)
xygrid1 <- data.frame(expand.grid(easting = xgrid, northing = xgrid),
                      partition = xygrid1$partition)
xytgrid <- expand.grid(easting = xgrid, northing = xgrid, time = tgrid) %>%
  arrange(time, easting, northing)
n <- nrow(xytgrid)
coords <- xytgrid

#################################
# wheel of fortune (partition1) #
#################################
# wind directions
n_northing <- 3
n_easting <- 3
nd <- floor(log10(max(n_northing, n_easting, n_time))) + 1
format <- paste0("%0", nd, ".0f,","%0", nd, ".0f,", "%0", nd, ".0f")

coords_ptt <- coords %>% left_join(xygrid1, by = c("easting", "northing"))
coords_ptt$row <- ifelse(coords_ptt$partition %in% c(8,2,6), 1,
                         ifelse(coords_ptt$partition %in% c(1,4,5), 2, 3))
coords_ptt$col <- ifelse(coords_ptt$partition %in% c(5,8,9), 1,
                         ifelse(coords_ptt$partition %in% 1:3, 2, 3))
coords_ptt$time_d <- as.numeric(as.character(coords_ptt$time*(n_time-1) + 1))
coords_ptt$partition <- sprintf(format,
                                coords_ptt$row,
                                coords_ptt$col,
                                coords_ptt$time_d)

##########################
# rectangle (partition2) #
##########################
# assign partitions
easting_cut <- seq(0, 1, length = n_easting + 1)
northing_cut <- seq(0, 1, length = n_northing + 1)

# wind directions
coords_ptt2 <- coords
coords_ptt2$row <- (n_northing+1) -
  as.numeric(cut_interval(coords_ptt2$northing,
                          n = n_northing,
                          labels = 1:n_northing))                               # 1 from the top
coords_ptt2$col <- as.numeric(cut_interval(coords_ptt2$easting,
                                           n = n_easting,
                                           labels = 1:n_easting))
coords_ptt2$time_d <- as.numeric(as.character(coords_ptt2$time*(n_time-1) + 1))
coords_ptt2$partition <- sprintf(format,
                                 coords_ptt2$row,
                                 coords_ptt2$col,
                                 coords_ptt2$time_d)

#######
# Cov #
#######
# true parameter values
a <- 0 # the smaller the longer the cov lasts
c <- 0.8 # the smaller the larger the cov
kappa <- 0
sig_sq <- 1

# wind directions
ptts <- sort(unique(coords_ptt$partition))
z_true1 <- matrix("W", nrow = length(ptts))
rownames(z_true1) <- ptts
z_true2 <- matrix("N", nrow = length(ptts))
rownames(z_true2) <- ptts
z_true3 <- matrix("NW", nrow = length(ptts))
rownames(z_true3) <- ptts

ptts2 <- sort(unique(coords_ptt2$partition))
z_true12 <- matrix("W", nrow = length(ptts2))
rownames(z_true12) <- ptts2
z_true22 <- matrix("N", nrow = length(ptts2))
rownames(z_true22) <- ptts2
z_true32 <- matrix("NW", nrow = length(ptts2))
rownames(z_true32) <- ptts2

# Ctilde given z
Cz_W <- Ctilde_ST(coords_ptt = coords_ptt, z = z_true1,
                  params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_N <- Ctilde_ST(coords_ptt = coords_ptt, z = z_true2,
                  params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_NW <- Ctilde_ST(coords_ptt = coords_ptt, z = z_true3,
                   params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))

Cz_W2 <- Ctilde_ST(coords_ptt = coords_ptt2, z = z_true12,
                   params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_N2 <- Ctilde_ST(coords_ptt = coords_ptt2, z = z_true22,
                   params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_NW2 <- Ctilde_ST(coords_ptt = coords_ptt2, z = z_true32,
                    params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))

# stationary
spatdist <- dist(coords_ptt[,c("easting", "northing")], method = "euclidean")
timedist <- dist(coords_ptt[,c("time")], method = "manhattan")
invaup1 <- 1/(a*as.matrix(timedist)+1)
stCor <- invaup1*exp(-c*as.matrix(spatdist)*(invaup1^(kappa/2)))
stCov <- sig_sq*stCor

#############
# Subfigure #
#############
dataall <- data.frame(easting = rep(coords_ptt$easting, 3),
                      northing = rep(coords_ptt$northing, 3),
                      time = rep(coords_ptt$time, 3),
                      cat = rep(c("stationarity", "partition1", "partition2"),
                                each = n)) %>%
  mutate(timename = paste0("time = ",round(time,3)))
C_marg <- .5*Cz_W + .1*Cz_N + .4*Cz_NW
C_marg2 <- .5*Cz_W2 + .1*Cz_N2 + .4*Cz_NW2

i_ref <- 1365
g1 <- dataall %>%
  mutate(cov = c(stCov[i_ref,], C_marg[i_ref,], C_marg2[i_ref,])) %>%
  filter(cat == "partition2") %>%
  ggplot() +
  geom_contour_filled(aes(easting, northing, z = cov), 
                      breaks = seq(0, 1, by = 0.1)) +
  scale_fill_scico_d(palette = "lapaz", direction = -1, begin = 0.3) +
  geom_point(data = coords_ptt[i_ref,] %>%
               mutate(timename = paste0("time = ", round(time,3))),
             aes(easting, northing), size = 2, col = "red") +
  facet_grid(cat ~ timename) +
  labs(x = "", y = "", fill = "Covariance",
       title = "G-BAG nonstationary covariance: directional") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "right") +
  geom_point(data = coords_ptt[c(i_ref + 2*ngrid^2 - 262,
                                 i_ref + 2*ngrid^2 + 262),] %>%
               mutate(timename = paste0("time = ", round(time,3)),
                      cat = "partition2"),
             aes(easting, northing), col = "black", size = 2) +
  geom_text(data = coords_ptt[c(i_ref + 2*ngrid^2 - 325,
                                i_ref + 2*ngrid^2 + 324),] %>%
              mutate(timename = paste0("time = ", round(time,3)),
                     cat = "partition2"),
            aes(easting, northing), col = "black",
            label = c(expression(s[i]), expression(s[j])),
            size = 6, parse = T)
# for (ext in extension) {
#   ggsave(plot = g1, paste0(path, "plots/cov_spatial_newparent", ext),
#          width = 7, height = 2.2)
# }

################################################################################

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
