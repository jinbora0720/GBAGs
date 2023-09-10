# Air Quality Analysis in CA #
# sensitivity to partitions; coarse partition
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf)
library(maps)
library(scales)
library(meshed)
library(INLA)
library(metR)
library(bags)
library(scico)

# path
path <- "~/GBAGs/"

# plot extensions
extension <- c(".pdf", ".eps", ".png")

########
# Data #
########
# call data
fire <- readRDS(paste0(path, "data/CA_fire_20211028.RDS"))
data0 <- readRDS(paste0(path, "data/CA_final_20211028.RDS")) %>%
  rename(time = date, logpm25 = logpm25_mean)
data0 %>% group_by(time) %>% summarise(n = n()) %>% summary()

# california map
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
ca_sf <- st_transform(states %>% filter(ID == "california"), 26911)
ca_range <- st_bbox(ca_sf)

timelist <- unique(data0$time)
timedf <- data.frame(time = timelist, time_d = 1:length(timelist))
data <- left_join(data0, timedf, by = "time") %>%
  rename(true_easting = easting, true_northing = northing) %>%
  mutate(time = (time_d - min(time_d))/(max(time_d) - min(time_d)),
         easting = true_easting/1000, northing = true_northing/1000)            # m to km
rm(data0)
tgrid <- sort(unique(data$time))
n_time <- length(tgrid)

# train data
# fully predict the space at the last time point
prob <- 0.8
set.seed(123)
data <- data %>% group_by(easting, northing) %>%
  mutate(train = ifelse(runif(1) < prob, 1, 0)) %>%
  ungroup(easting, northing)
n <- nrow(data)

tr_idx <- which(data$train == 1)
n_tr <- length(tr_idx)
x_tr <- data$dist_to_fire[tr_idx]
y_tr <- data$logpm25[tr_idx]
coords_tr <- data[tr_idx, c("easting", "northing", "time",
                             "time_d", "true_easting", "true_northing")]

tt_idx <- which(data$train == 0)
n_tt <- length(tt_idx)
x_tt <- data$dist_to_fire[tt_idx]
y_tt <- data$logpm25[tt_idx]
coords_tt <- data[tt_idx, c("easting", "northing", "time",
                             "time_d", "true_easting", "true_northing")]

# prediction data
xgrid <- seq(ca_range$xmin, ca_range$xmax, length = 80)
ygrid <- seq(ca_range$ymin, ca_range$ymax, length = 100)
xygrid <- expand.grid(easting = xgrid, northing = ygrid, time_d = 1:n_time)
xygrid_sf <- st_as_sf(xygrid, coords = c("easting", "northing"), crs = 26911)
xygrid <- xygrid %>% mutate(inCA = as.numeric(st_intersects(xygrid_sf, ca_sf$geom)))

coords_grid <- xygrid %>% filter(inCA == 1) %>%
  rename(true_easting = easting, true_northing = northing) %>%
  mutate(easting = true_easting/1000, northing = true_northing/1000,
         time = (time_d - min(time_d))/(max(time_d) - min(time_d)))
n_grid <- nrow(coords_grid)
x_grid <- rep(0, n_grid)
for (t in 1:n_time) {
  idx_tmp <- which(coords_grid$time == tgrid[t])
  data_tmp <- coords_grid[idx_tmp,]
  fire_tmp <- fire %>% filter(date == timelist[t])
  nn_tmp <- RANN::nn2(fire_tmp[,c("easting","northing")],
                      data_tmp[,c("true_easting","true_northing")], k=1)
  x_grid[idx_tmp] <- nn_tmp$nn.dists
}

# pre-processing
# CA area
# width 890.3218km = (ca_range$xmax - ca_range$xmin)/1000
# height 1075.088km = (ca_range$ymax - ca_range$ymin)/1000

n_easting <- 9 
n_northing <- 11
easting_cut <- seq(ca_range$xmin, ca_range$xmax, length = n_easting + 1)
northing_cut <- seq(ca_range$ymin, ca_range$ymax, length = n_northing + 1)
easting_midval <- 0.5*easting_cut[1:n_easting] + 0.5*easting_cut[-1]
northing_midval <- 0.5*northing_cut[1:n_northing] + 0.5*northing_cut[-1]

###########
# Methods #
###########
# directions
directions <- c("W", "NW", "N", "NE")

# mcmc
mcmc <- list(save = 1000, burn = 15000, thin = 15)

# prior
la = lc <- 0
ua = uc <- 1000

# scaling X
x_scale <- scale(c(x_tr, x_tt, x_grid), center = T, scale = T)
X_tr <- as.matrix(x_scale[1:n_tr])
X_pred <- as.matrix(x_scale[n_tr+1:(n_tt+n_grid)])

# ## bags ##
# out <- bag(y = y_tr, X = X_tr,
#            coords = as.matrix(coords_tr[,c("easting", "northing", "time")]),
#            X_pred = X_pred,
#            coords_pred = as.matrix(rbind(coords_tt[,c("easting", "northing", "time")],
#                                          coords_grid[,c("easting", "northing", "time")])),
#            n_partition = c(n_easting, n_northing, n_time),
#            breaks_partition = list(breaks_easting = NULL,
#                                    breaks_northing = NULL,
#                                    breaks_time = NULL),
#            directions = directions,
#            init = list(tau_sq = NULL,
#                        sig_sq = NULL,
#                        w = NULL,
#                        z = NULL,
#                        psi = NULL,
#                        Sn = NULL),
#            hyper = list(at = NULL, bt = NULL,
#                         as = NULL, bs = NULL,
#                         la = la, ua = ua,
#                         lc = lc, uc = uc,
#                         mu0 = NULL, invV0 = NULL),
#            mcmc = mcmc,
#            n_threads = 10,
#            seed = 123,
#            verbose = TRUE,
#            save_data = TRUE,
#            save_est = TRUE,
#            debug = list(psi_fixed = FALSE, z_fixed = FALSE))
# 
# saveRDS(out, paste0(path, "application/CA_daily_bag_v3_all.RDS"))
# 
# # effective sample size
# w_ess <- c(apply(out$w_save, 1, posterior::ess_basic),
#            apply(out$w_pred_save, 1, posterior::ess_basic))
# y_ess <- c(apply(out$y_save, 1, posterior::ess_basic),
#            apply(out$y_pred_save, 1, posterior::ess_basic))
# 
# # estimation
# what <- rowMeans(out$w_save)
# yhat <- rowMeans(out$y_save)
# yhat_qt <- apply(out$y_save, 1, function(x)
#   quantile(x, probs = c(0.025, 0.975)))
# yhat_low <- yhat_qt[1,]
# yhat_hi <- yhat_qt[2,]
# 
# # prediction
# w_bag <- rowMeans(out$w_pred_save)
# y_bag <- rowMeans(out$y_pred_save)
# y_bag_qt <- apply(out$y_pred_save, 1, function(x)
#   quantile(x, probs = c(0.025, 0.975)))
# y_bag_low <- y_bag_qt[1,]
# y_bag_hi <- y_bag_qt[2,]
# 
# saveRDS(list(tau_sq_save = out$tau_sq_save,
#              beta_save = out$beta_save,
#              sig_sq_save = out$sig_sq_save,
#              psi_save = out$psi_save,
#              z_save = out$z_save,
#              est_time_per_iter = out$est_time/out$est_iter,
#              pred_time_per_iter = out$pred_time/out$pred_iter,
#              what = what,
#              yhat = yhat,
#              yhat_low = yhat_low,
#              yhat_hi = yhat_hi,
#              w_bag = w_bag,
#              y_bag = y_bag,
#              y_bag_low = y_bag_low,
#              y_bag_hi = y_bag_hi,
#              w_ess = w_ess,
#              y_ess = y_ess,
#              coords_ptt = out$coords_ptt),
#         paste0(path, "application/CA_daily_bag_v3.RDS"))

###########
# Results #
###########
# call results
bagres <- readRDS(paste0(path, "application/CA_daily_bag_v3.RDS"))

###############
# Convergence #
###############
# traceplot of tau_sq
tracep <- data.frame(tau_sq = bagres$tau_sq_save, iter = 1:mcmc$save) %>%
  ggplot() +
  geom_line(aes(iter, tau_sq)) +
  labs(y = expression(tau^2), x = "") +
  theme(plot.margin = margin(b = -5, r = -5, t = 0, l = 0))

# running mean of beta
runmeanp <- data.frame(beta = cumsum(bagres$beta_save)/1:mcmc$save,
                       iter = 1:mcmc$save) %>%
  ggplot() +
  geom_line(aes(iter, beta)) +
  labs(y = expression(paste("Running mean of ", beta)), x = "") +
  theme(plot.margin = margin(b = 0, r = -5, t = -5, l = 0))

# ESS of w and y
essp <- data.frame(rbind(coords_tr[,c("true_easting", "true_northing", "time_d")],
                         coords_tt[,c("true_easting", "true_northing", "time_d")],
                         coords_grid[,c("true_easting", "true_northing", "time_d")]),
                   ess = bagres$y_ess) %>%
  filter(time_d == 5) %>%
  ggplot() +
  geom_point(aes(true_easting, true_northing, color = ess), size = 0.5) +
  geom_sf(data = ca_sf, fill = "NA") +
  labs(x = "", y = "", color = "ESS of y") +
  scale_color_scico(palette = "cork") +
  theme(legend.position = "bottom",
        plot.margin = margin(b = 0, r = 0, t = 0, l = -5),
        legend.margin = margin(b = 0, r = 0, t = -5, l = 0))

conv1 <- gridExtra::grid.arrange(tracep, runmeanp, nrow = 2)
conv2 <- gridExtra::grid.arrange(conv1, essp, ncol = 2)

# z Chisq test
prop <- 0.35
n_samp <- ceiling(0.25*mcmc$save)
n1 = n2 <- n_samp

z_pvalue <- rep(0, nrow(bagres$z_save))
for (ii in 1:nrow(bagres$z_save)) {
  samp1 <- bagres$z_save[ii,1:n_samp]
  samp2 <- bagres$z_save[ii,(mcmc$save - n_samp) + 1:n_samp]

  N1 <- c(sum(samp1 == directions[1]), sum(samp1 == directions[2]),
          sum(samp1 == directions[3]), sum(samp1 == directions[4]))
  p1 <- N1/n1
  N2 <- c(sum(samp2 == directions[1]), sum(samp2 == directions[2]),
          sum(samp2 == directions[3]), sum(samp2 == directions[4]))
  p2 <- N2/n2
  p <- (N1 + N2)/(n1 + n2)
  R <- which(p != 0)
  nR <- sum(p != 0)
  X_sq <- n1*sum((p1[R] - p[R])^2/p[R]) + n2*sum((p2[R] - p[R])^2/p[R])
  z_pvalue[ii] <- pchisq(X_sq, df = (nR-1)*(2-1), lower.tail = FALSE)
}
1 - mean(z_pvalue < 0.05)
