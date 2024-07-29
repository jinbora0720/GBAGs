# Nonstationarity of G-BAGs #
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(scico)
library(bags)
# devtools::install_github("marysalvana/Lagrangian")
library(Lagrangian)

# path
path <- "~/GBAGs/"

# plot extensions
extension <- c(".pdf", ".eps", ".png")

###############
# Proposition #
###############
# spatiotemporal w/ Apanasovich and Genton (2010) 
N <- 100
a <- 0.07
c <- 0.01
kappa <- 1
sig_sq <- 1
h <- seq(0, 12, length = N)
u <- seq(0, 12, length = N)
hu <- expand.grid(h, u)

basedata <- data.frame(h = hu[,1], u = hu[,2]) %>% 
  mutate(invaup1 = 1/((a*u^2+1)^kappa)) %>% 
  mutate(z = sig_sq*invaup1/((a*h^2+1)^(kappa/2))*exp(-c*h^2*invaup1-c*u^2)) %>% 
  select(-invaup1)

p <- c(0.51, 0.55, seq(0.6, 0.9, by = 0.1), 0.95, 0.99)
Ctilde12 = Ctilde13 <- matrix(0, nrow = N^2, ncol = length(p))
for (i in 1:length(p)) {
  Ctilde12[,i] <- p[i]*basedata$z
  Ctilde13[,i] <- (1-p[i])*basedata$z
}

g0 <- data.frame(h = rep(basedata$h, times = length(p)), 
                 u = rep(basedata$u, times = length(p)), 
                 pi = rep(paste("pi", "==", p), each = N^2), 
                 nonstCov12 = as.numeric(Ctilde12), 
                 nonstCov13 = as.numeric(Ctilde13)) %>% 
  pivot_longer(c("nonstCov12", "nonstCov13"), 
               names_to = "what", values_to = "Cov") %>% 
  filter(pi %in% paste("pi", "==", c(0.51, 0.55, 0.7))) %>% 
  ggplot() + 
  facet_wrap(~ pi, ncol = 3, labeller = labeller(pi = label_parsed)) + 
  geom_contour(aes(x = h, y = u, z = Cov, color = what, linetype = what), 
               breaks = 0.1) +
  scale_color_manual(values = c("nonstCov12" = "#E69F00", 
                                "nonstCov13" = "#009E73"), 
                     labels = c(expression(paste(s[1]," & ",s[2])), 
                                expression(paste(s[1]," & ",s[3])))) +
  scale_linetype_manual(values = c("nonstCov12" = 1, 
                                   "nonstCov13" = 2), 
                        labels = c(expression(paste(s[1]," & ",s[2])), 
                                   expression(paste(s[1]," & ",s[3])))) + 
  labs(x = "||h||", y = "|u|", 
       title = "G-BAG nonstationary covariance: directional", 
       color = "", linetype = "") +
  theme(legend.position = "bottom", 
        plot.margin = margin(t = 1, l = 1, r = 8, b = -2),
        legend.margin = margin(b = 0, r = 0, t = -10, l = 0), 
        legend.text = element_text(size = 11))
# for (ext in extension) {
#   ggsave(plot = g0, paste0(path, "plots/cov_proposition", ext),
#          width = 7, height = 2.5)
# }

################################################################################

##################################
# Different partitioning schemes #
##################################
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

#-- wheel of fortune (partition1) --#
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

#-- rectangle (partition2) --#
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

#-- Cov --#
# true parameter values
a <- .7 # the smaller the longer the cov lasts
c <- .8 # the smaller the larger the cov
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
Cz_W <- Ctilde(coords_ptt = coords_ptt, z = z_true1,
               params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_N <- Ctilde(coords_ptt = coords_ptt, z = z_true2,
               params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_NW <- Ctilde(coords_ptt = coords_ptt, z = z_true3,
                params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))

Cz_W2 <- Ctilde(coords_ptt = coords_ptt2, z = z_true12,
                params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_N2 <- Ctilde(coords_ptt = coords_ptt2, z = z_true22,
                params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))
Cz_NW2 <- Ctilde(coords_ptt = coords_ptt2, z = z_true32,
                 params = list(sig_sq = sig_sq, a = a, c = c, kappa = kappa))

# stationary
spatdist <- dist(coords_ptt[,c("easting", "northing")], method = "euclidean")
timedist <- dist(coords_ptt[,c("time")], method = "manhattan")
invaup1 <- 1/(a*as.matrix(timedist)+1)
stCor <- invaup1*exp(-c*as.matrix(spatdist)*(invaup1^(kappa/2)))
stCov <- sig_sq*stCor

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
  filter(cat != "stationarity") %>%
  ggplot() +
  geom_contour_filled(aes(easting, northing, z = cov)) +
  scale_fill_scico_d(palette = "lapaz", direction=-1) +
  geom_point(data = coords_ptt[i_ref,] %>%
               mutate(timename = paste0("time = ", round(time,3))),
             aes(easting, northing), size = 2, col = "red") +
  facet_grid(cat ~ timename) +
  labs(x = "", y = "", fill = "Covariance",
       title = "G-BAG nonstationary covariance: directional") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "bottom") +
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

g2 <- dataall %>%
  mutate(cov = c(stCov[i_ref,], C_marg[i_ref,], C_marg2[i_ref,])) %>%
  filter(cat == "stationarity") %>%
  ggplot() +
  geom_contour_filled(aes(easting, northing, z = cov),
                      breaks = seq(0,1,by=0.1)) +
  scale_fill_manual(values =
                      c("#A5A58E", "#7B9F9E", "#518DA6", "#3572A3",
                        "#295394", "#22327F", "#190C65")) +
  geom_point(data = coords_ptt[i_ref,] %>%
               mutate(timename = paste0("time = ", round(time,3))),
             aes(easting, northing), size = 2, col = "red") +
  facet_grid(cat ~ timename) +
  labs(x = "", y = "", fill = "Covariance",
       title = "Base covariance: symmetric") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "none")

gg <- ggpubr::ggarrange(g1, g2, nrow = 2, common.legend = T,
                        legend = "right", heights = c(2,1.3))
# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/cov_spatial", ext),
#          width = 7, height = 5.5)
# }
################################################################################

############################################
# Cov between left and right under W winds #
############################################
# specify number of grid on each axis
ngrid3 = n_northing = n_easting <- 3
n_time3 <- 30
xgrid3 <- seq(0, 1, length = ngrid3)
tgrid3 <- seq(0, 1, length = n_time3)
xytgrid3 <- expand.grid(easting = xgrid3, northing = xgrid3, time = tgrid3) %>%
  arrange(time, easting, northing)
n3 <- nrow(xytgrid3)
coords3 <- xytgrid3

# wind directions
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

#-- Cov --#
# true parameter values
a3 <- 2 # the smaller the longer the cov lasts
c <- .8 # the smaller the larger the cov
kappa <- 0
sig_sq <- 1

# wind directions
ptts3 <- sort(unique(coords_ptt3$partition))
z_true13 <- matrix("W", nrow = length(ptts3))
rownames(z_true13) <- ptts3

# Ctilde given z
Cz_W3 <- Ctilde(coords_ptt = coords_ptt3, z = z_true13,
                params = list(sig_sq = sig_sq, a = a3, c = c, kappa = kappa))

# stationary
spatdist3 <- dist(coords_ptt3[,c("easting", "northing")], method = "euclidean")
timedist3 <- dist(coords_ptt3[,c("time")], method = "manhattan")
invaup13 <- 1/(a3*as.matrix(timedist3)+1)
stCor3 <- invaup13*exp(-c*as.matrix(spatdist3)*(invaup13^(kappa/2)))
stCov3 <- sig_sq*stCor3

i_ref <- 95
Eidx_ref <- (i_ref - ngrid3) + (-n_time3:n_time3)*ngrid3^2
Eidx_ref <- Eidx_ref[Eidx_ref < n3 & Eidx_ref > 0]
Widx_ref <- (i_ref + ngrid3) + (-n_time3:n_time3)*ngrid3^2
Widx_ref <- Widx_ref[Widx_ref < n3 & Widx_ref > 0]

poly_df_base <- data.frame(x = c(0, 1, 0), y = c(0, 0.5, 1))
poly_df_list <- list()
poly_centers <- expand.grid(seq(0.1, 0.9, 0.2), seq(0.1, 0.9, 0.2))
for(i in 1:nrow(poly_centers)){
  poly_df_list[[i]] <- poly_df_base %>% mutate(group = i,
                                               x = 0.05*x + poly_centers[i,1],
                                               y = 0.05*y + poly_centers[i,2])
}
poly_df <- bind_rows(poly_df_list) %>% mutate(xy = (x==.5)&(y==.5)) %>%
  filter(xy != TRUE) %>%
  mutate(x = x-0.025)

gg1 <- coords_ptt3[Widx_ref,] %>%
  mutate(timename = paste0("time = ", round(time,3))) %>%
  filter(time == coords_ptt3[i_ref,"time"]) %>%
  ggplot() +
  geom_point(aes(easting, northing), col = "#D55E00", size = 6) +
  geom_point(data = coords_ptt3[Eidx_ref,] %>%
               mutate(timename = paste0("time = ", round(time,3))) %>%
               filter(time == coords_ptt3[i_ref,"time"]),
             aes(easting, northing), col = "#0072B2", size = 6) +
  geom_point(data = coords_ptt3[i_ref,] %>%
               mutate(timename = paste0("time = ", round(time,3))),
             aes(easting, northing), size = 6) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "", title = " ") +
  geom_polygon(data = poly_df, aes(x = x, y = y, group = group), alpha = .2)

cov_df <- data.frame(dist = rep(tgrid3,2),
                     dir = c(rep("W", each = length(Widx_ref)),
                             rep("E", each = length(Eidx_ref))),
                     cov = c(Cz_W3[i_ref, Widx_ref], Cz_W3[i_ref, Eidx_ref]))
west_max <- cov_df %>% filter(dir == "W") %>%
  arrange(desc(cov)) %>% `[`(1,"dist")
east_max <- cov_df %>% filter(dir=="E") %>%
  arrange(desc(cov)) %>% `[`(1,"dist")

gg2 <- cov_df %>%
  ggplot() +
  geom_rect(data = data.frame(xmin = coords_ptt3[i_ref,"time"],
                              xmax = west_max,
                              ymin = 0, ymax = .8),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            alpha = .2, fill = "#D55E00") +
  geom_rect(data = data.frame(xmin = east_max,
                              xmax = coords_ptt3[i_ref,"time"],
                              ymin = 0, ymax = .8),
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax), alpha = .2, fill = "#0072B2") +
  geom_line(aes(dist, cov, col = dir, linetype = dir), linewidth = 1) +
  scale_color_manual(values = c("E" = "#0072B2", "W" = "#D55E00"),
                     labels = c("left", "right")) +
  scale_linetype_manual(values = c("E" = "dashed", "W" = "solid"),
                        labels = c("left", "right")) +
  geom_vline(xintercept = coords_ptt3[i_ref,"time"]) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title=element_text(size = 12),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Covariance", #title = "Winds from W (left to right)",
       title = "G-BAG nonstationary covariance: directional",
       color = "Relative location to the reference point",
       linetype = "Relative location to the reference point") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(0, east_max, coords_ptt3[i_ref, "time"],
                              west_max, 1),
                     labels = c(0, latex2exp::TeX("t* - $\\Delta$"),
                              "t*", latex2exp::TeX("t* + $\\Delta$"), 1))

gg3 <- data.frame(dist = rep(tgrid3,2),
                  dir = c(rep("E", each = length(Eidx_ref)),
                          rep("W", each = length(Widx_ref))),
                  cov = c(stCov3[i_ref, Eidx_ref]-0.004,
                          stCov3[i_ref, Widx_ref]+0.004)) %>%
  ggplot() +
  geom_line(aes(dist, cov, col = dir, linetype = dir), linewidth = 1) +
  scale_color_manual(values = c("E" = "#0072B2", "W" = "#D55E00"),
                     labels = c("left", "right")) +
  scale_linetype_manual(values = c("E" = "dashed", "W" = "solid"),
                        labels = c("left", "right")) +
  geom_vline(xintercept = coords_ptt3[i_ref,"time"]) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Covariance", #title = "Winds from W (left to right)",
       title = "Base covariance: symmetric",
       color = "Relative location to the reference point",
       linetype = "Relative location to the reference point") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = c(0, 0.2, coords_ptt3[i_ref, "time"], 0.6, 0.8, 1),
                     labels = c(0, 0.2, "t*", 0.6, 0.8, 1))

ggg <- ggpubr::ggarrange(gg1,
                         gg2,
                         gg3 + theme(axis.title.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.text.y = element_blank()),
                         nrow = 1, widths = c(1, 1.6, 1.5),
                         common.legend = T, legend = "bottom")
# for (ext in extension) {
#   ggsave(plot = ggg, paste0(path, "plots/cov_time", ext),
#          width = 11, height = 3)
# }

################################################################################

##############
# Lagrangian #
##############
# locations
xgrid <- seq(0, 1, length = 6)
ygrid <- xgrid
tgrid <- seq(0, 1, length = 15)
coords <- expand.grid(easting = xgrid, northing = ygrid, time = tgrid) %>%
  arrange(time, easting, northing)
n <- nrow(coords)

# rotate 
theta <- pi/4
rot45_mat <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 
                    2, 2, byrow = T)
coords_45 <- data.frame(as.matrix(coords[,1:2])%*%rot45_mat, time = coords[,3]) %>%
  rename("easting" = "X1", "northing" = "X2")

# seed 
seed <- 710
set.seed(seed)

# -- Cov --#
# frozen Lagrangian nonstationary spatiotemporal covariance
## spatially-varying
cov_mat_sv <- cov_uni_lagrangian(location = coords, sigma2 = 1, scale = 1, nu = 1/2,
                                 vx = 3.1, vy = 3.1, nonstat = 2) 
# G-BAG
n_easting <- 11 
n_northing <- 1
n_time <- 15

nd <- floor(log10(max(n_easting, n_northing, n_time))) + 1
coords_ptt <- data.frame(coords_45)
coords_ptt$row <- (n_northing+1) -
  as.numeric(cut_interval(coords_ptt$northing,
                          n = n_northing,
                          labels = 1:n_northing))                           
coords_ptt$col <- as.numeric(cut_interval(coords_ptt$easting,
                                          n = n_easting,
                                          labels = 1:n_easting))
coords_ptt$time_d <- as.numeric(as.character(coords_ptt$time*(n_time-1) + 1))
coords_ptt$partition <- sprintf(bags:::ptt_format(nd),
                                coords_ptt$row,
                                coords_ptt$col,
                                coords_ptt$time_d)
ptts <- sort(unique(coords_ptt$partition))

C_sv <- Ctilde(coords_ptt = coords_ptt,
               z = matrix(rep("W", length(ptts)), ncol = 1),
               params = list(sig_sq = 1, a = 3.5, c = 0.96, kappa = 0.8)) 

time_label <- paste0("time = ", round(tgrid, 3))
names(time_label) <- tgrid

plt1 <- data.frame(rbind(coords, coords), cov = c(C_sv[1,], cov_mat_sv[1,]),
                   type = rep(c("G-BAG", "Lagrangian"), each = n)) %>%
  filter(type == "G-BAG") %>% 
  ggplot() +
  geom_contour_filled(aes(easting, northing, z = cov)) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  # geom_raster(aes(easting, northing, fill = cov)) +
  # scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  facet_wrap(~time, labeller = labeller(time = time_label), nrow = 3) +
  labs(x = "", y = "", fill = "Covariance", title = "G-BAG")

plt2 <- data.frame(rbind(coords, coords), cov = c(C_sv[1,], cov_mat_sv[1,]),
                   type = rep(c("G-BAG", "Lagrangian"), each = n)) %>%
  filter(type == "Lagrangian") %>% 
  ggplot() +
  geom_contour_filled(aes(easting, northing, z = cov)) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  # geom_raster(aes(easting, northing, fill = cov)) +
  # scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  facet_wrap(~time, labeller = labeller(time = time_label), nrow = 3) +
  labs(x = "", y = "", fill = "Covariance", title = "Lagrangian")

cor(c(cov_mat_sv), c(C_sv))
plt3 <- data.frame(x = c(cov_mat_sv), y = c(C_sv)) %>%
  ggplot() +
  geom_point(aes(x,y), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = "red", 
              linewidth = 1, linetype = "dashed") +
  labs(x = "Lagrangian", y = "G-BAG") + 
  coord_fixed()

plt <- ggpubr::ggarrange(plt3, plt1, plt2, nrow = 1, common.legend = T,
                         legend = "right", widths = c(0.8, 1, 1))
# for (ext in extension) {
#   ggsave(plot = plt, paste0(path, "plots/cov_Lagrangian", ext),
#          width = 20, height = 6)
# }
