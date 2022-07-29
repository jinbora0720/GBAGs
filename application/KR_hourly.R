###################################################
# On Feb 4 - 5 (8pm - 7am), NW winds seem obvious #
# original scale                                  #
###################################################
# Air Quality Analysis in South Korea #
# Table S1, Figure S2.
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf)
library(meshed)
library(bags)
library(INLA)

# path
path <- "~/GBAGs/"

# plot extensions
extension <- c(".pdf", ".eps", ".png")

########
# Data #
########
# call data
data <- readRDS(paste0(path, "data/KR_Feb_pm25.rds"))
KR_sf <- readRDS(paste0(path, "data/KR_sf.rds"))

# 2020-02-04 ~ 2020-02-05 wind direction
timelist <- c(2020020420:2020020424, 2020020501:2020020507)                     # without much traffic

data_short <- data %>% filter(time %in% timelist) %>%
  group_by(station_id, time) %>%
  summarise(long = unique(long),
            lat = unique(lat),
            PM25 = mean(PM25, na.rm = TRUE)) %>%
  mutate(logpm25 = log(PM25))

data_short_full <- data_short %>%
  mutate(n = length(unique(time))) %>%
  filter(n == length(timelist)) %>% dplyr::select(-n) %>%
  ungroup()

# coordinates
coords <- unique(data_short_full[,c("station_id", "long", "lat")])
coords_sf <- st_transform(st_as_sf(coords, coords = c("long", "lat"),
                                   crs = 4326), crs = st_crs(KR_sf))
coords_prj <- data.frame(station_id = coords$station_id,
                         coords_sf %>% st_coordinates()) %>%
  rename(easting = X, northing = Y)
data_short_full2 <- left_join(data_short_full, coords_prj, by = "station_id")

# KR box
KR_range <- st_bbox(KR_sf)
KR_range <- data.frame(xmin = as.numeric(KR_range$xmin),
                       ymin = as.numeric(KR_range$ymin),
                       xmax = 1190000,
                       ymax = as.numeric(KR_range$ymax))

timedf <- data.frame(time = timelist, time_d = 1:length(timelist))
data0 <- left_join(data_short_full2, timedf, by = "time") %>%
  rename(true_easting = easting, true_northing = northing) %>%
  mutate(time = (time_d - min(time_d))/(max(time_d) - min(time_d)),
         easting = true_easting/1000, northing = true_northing/1000)            # m to km

tgrid <- sort(unique(data0$time))
n_time <- length(tgrid)

# train data
prob <- 0.7
set.seed(123)
data0 <- data0 %>% group_by(easting, northing) %>%
  mutate(train = ifelse(runif(1) < prob, 1, 0)) %>%
  ungroup(easting, northing)
n <- nrow(data0)

tr_idx <- which(data0$train == 1)
n_tr <- length(tr_idx)
y_tr <- data0$logpm25[tr_idx]
coords_tr <- data0[tr_idx, c("easting", "northing", "time",
                             "time_d", "true_easting", "true_northing")]

tt_idx <- which(data0$train == 0)
n_tt <- length(tt_idx)
y_tt <- data0$logpm25[tt_idx]
coords_tt <- data0[tt_idx, c("easting", "northing", "time",
                             "time_d", "true_easting", "true_northing")]

# prediction data
xgrid <- seq(KR_range$xmin, KR_range$xmax, length = 110)
ygrid <- seq(KR_range$ymin, KR_range$ymax, length = 150)
xygrid <- expand.grid(easting = xgrid, northing = ygrid, time_d = 1:n_time)
xygrid_sf <- st_as_sf(xygrid, coords = c("easting", "northing"), crs = st_crs(KR_sf))
xygrid <- xygrid %>%
  mutate(inKR = as.numeric(st_intersects(xygrid_sf, KR_sf$geometry)))

coords_grid <- xygrid %>% filter(!is.na(inKR)) %>%
  rename(true_easting = easting, true_northing = northing) %>%
  mutate(easting = true_easting/1000, northing = true_northing/1000,
         time = (time_d - min(time_d))/(max(time_d) - min(time_d)))
n_grid <- nrow(coords_grid)

# pre-processing
# south korea area
# width 443.890km = (KR_range$xmax - KR_range$xmin)/1000
# height 609.690km = (KR_range$ymax - KR_range$ymin)/1000
n_easting <- 17 # ~ 26km
n_northing <- 23 # ~ 27km
easting_cut <- seq(KR_range$xmin, KR_range$xmax, length = n_easting + 1)
northing_cut <- seq(KR_range$ymin, KR_range$ymax, length = n_northing + 1)
easting_midval <- 0.5*easting_cut[1:n_easting] + 0.5*easting_cut[-1]
northing_midval <- 0.5*northing_cut[1:n_northing] + 0.5*northing_cut[-1]

###########
# Methods #
###########
# directions
directions <- c("SW", "W", "NW", "N")

# mcmc
mcmc <- list(save = 1000, burn = 15000, thin = 15)

# prior
tmaxdist <- max(dist(coords_tr[,3], method = "manhattan"))
la <- (1/0.95-1)/tmaxdist                                                       # correlation = 0.95 at the maximum distance
ua <- (1/0.05-1)/(0.5*tmaxdist)                                                 # correlation = 0.05 at half the maximum distance
spmaxdist <- max(dist(coords_tr[,1:2]))
lc <- -log(0.5)/spmaxdist                                                       # correlation = 0.5 at the maximal distance
uc <- -log(0.01)/(0.1*spmaxdist)                                                # correlation = 0.01 at 10% of the maximum distance

# ## bags ##
# out <- bag(y = y_tr, X = NULL,
#            coords = as.matrix(coords_tr[,c("easting", "northing", "time")]),
#            X_pred = NULL,
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
# saveRDS(out, paste0(path, "application/KR_hourly_bag_all.RDS"))

# ## mgp ##
# y <- as.matrix(c(y_tr, rep(NA, n_tt + n_grid)))
# x <- matrix(1, nrow = n + n_grid, ncol = 1)
#
# coords_mgp <- rbind(coords_tr[,c("easting", "northing", "time")],
#                     coords_tt[,c("easting", "northing", "time")],
#                     coords_grid[,c("easting", "northing", "time")]) %>%
#   mutate(easting = -easting +1)
#
# set.seed(123)
# meshout <- spmeshed(family = c("gaussian"),
#                     y = y, x = x,
#                     coords = coords_mgp,
#                     k = 1,
#                     axis_partition = c(n_easting, n_northing, n_time),
#                     prior = list(phi = c(min(lc,la), max(uc,ua))),
#                     settings = list(adapting = TRUE, forced_grid = NULL,
#                                     cache = NULL, ps = TRUE, saving = TRUE),
#                     starting = list(beta = NULL, tausq = NULL, theta = NULL,
#                                     lambda = matrix(1,1,1), w = NULL,
#                                     nu = NULL, mcmcsd = .05, mcmc_startfrom = 0),
#                     debug = list(sample_beta = TRUE, sample_tausq = TRUE,
#                                  sample_theta = TRUE, sample_w = TRUE,
#                                  sample_lambda = FALSE,
#                                  verbose = TRUE, debug = FALSE),
#                     n_samples = mcmc$save, n_burn = mcmc$burn, n_thin = mcmc$thin,
#                     verbose = 100,
#                     n_threads = 10)
#
# saveRDS(meshout, paste0(path, "application/KR_hourly_mgp_all.RDS"))

# # split KR_hourly_*_all.RDS due to size
# ## bags ##
# # effective sample size
# w_ess <- c(apply(out$w_save, 1, posterior::ess_basic),
#            apply(out$w_pred_save, 1, posterior::ess_basic))
# y_ess <- c(apply(out$y_save, 1, posterior::ess_basic),
#            apply(out$y_pred_save, 1, posterior::ess_basic))
#
# # estimation
# yhat <- rowMeans(out$y_save)
# yhat_qt <- apply(out$y_save, 1, function(x)
#   quantile(x, probs = c(0.025, 0.975)))
# yhat_low <- yhat_qt[1,]
# yhat_hi <- yhat_qt[2,]
#
# # prediction
# y_bag <- rowMeans(out$y_pred_save)
# y_bag_qt <- apply(out$y_pred_save, 1, function(x)
#   quantile(x, probs = c(0.025, 0.975)))
# y_bag_low <- y_bag_qt[1,]
# y_bag_hi <- y_bag_qt[2,]
#
# saveRDS(list(tau_sq_save = out$tau_sq_save,
#              sig_sq_save = out$sig_sq_save,
#              psi_save = out$psi_save,
#              z_save = out$z_save,
#              est_time_per_iter = out$est_time/out$est_iter,
#              pred_time_per_iter = out$pred_time/out$pred_iter,
#              yhat = yhat,
#              yhat_low = yhat_low,
#              yhat_hi = yhat_hi,
#              y_bag = y_bag,
#              y_bag_low = y_bag_low,
#              y_bag_hi = y_bag_hi,
#              w_ess = w_ess,
#              y_ess = y_ess,
#              coords_ptt = out$coords_ptt),
#         paste0(path, "application/KR_hourly_bag.RDS"))
#
# ## mgp ##
# y_mgp <- meshout$yhat_mcmc %>% summary_list_mean()
# y_mgp_low <- meshout$yhat_mcmc %>% summary_list_q(q = 0.025)
# y_mgp_hi <- meshout$yhat_mcmc %>% summary_list_q(q = 0.975)
#
# saveRDS(list(tausq_mcmc = meshout$tausq_mcmc[,mcmc$thin*(1:mcmc$save)],
#              sigsq_mcmc = meshout$theta_mcmc[4,,mcmc$thin*(1:mcmc$save)],
#              psi_mcmc = meshout$theta_mcmc[-4,,mcmc$thin*(1:mcmc$save)],        # temporal decay, spatial decay, interaction
#              time_per_iter = meshout$mcmc_time/meshout$mcmc,
#              coordsdata = meshout$coordsdata,
#              y_mgp = y_mgp,
#              y_mgp_low = y_mgp_low,
#              y_mgp_hi = y_mgp_hi),
#         paste0(path, "application/KR_hourly_mgp.RDS"))

###########
# Results #
###########
# call results
bagres <- readRDS(paste0(path, "application/KR_hourly_bag.RDS"))
mgpres <- readRDS(paste0(path, "application/KR_hourly_mgp.RDS"))

# time per iteration in seconds
bagres$est_time_per_iter + bagres$pred_time_per_iter

############
# Table S2 #
############
# parameter estimation
## tausq
mean(bagres$tau_sq_save) %>% round(3)
quantile(bagres$tau_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(mgpres$tausq_mcmc) %>% round(3)
quantile(mgpres$tausq_mcmc, probs = c(0.025, 0.975)) %>% round(3)

## sigsq
mean(bagres$sig_sq_save) %>% round(3)
quantile(bagres$sig_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(mgpres$sigsq_mcmc) %>% round(3)
quantile(mgpres$sigsq_mcmc, probs = c(0.025, 0.975)) %>% round(3)

## psi
rowMeans(bagres$psi_save) %>% round(3)
apply(bagres$psi_save, 1, function(x)
  quantile(x, probs = c(0.025, 0.975))) %>% round(3)

rowMeans(mgpres$psi_mcmc) %>% round(3)
apply(mgpres$psi_mcmc, 1, function(x)
  quantile(x, probs = c(0.025, 0.975))) %>% round(3)

# prediction summaries
res <- data.frame(bagres$coords_ptt[,c("easting", "northing", "time")],
                  y = c(y_tr, y_tt, rep(NA, n_grid)),
                  y_bag = c(bagres$yhat, bagres$y_bag),
                  y_bag_low = c(bagres$yhat_low,bagres$ y_bag_low),
                  y_bag_hi = c(bagres$yhat_hi, bagres$y_bag_hi)) %>%
  left_join(data.frame(mgpres$coordsdata,
                       y_mgp = mgpres$y_mgp,
                       y_mgp_low = mgpres$y_mgp_low,
                       y_mgp_hi = mgpres$y_mgp_hi) %>%
              mutate(easting = -easting+1),
            by = c("easting", "northing", "time")) %>%
  mutate(bag_inCI = ifelse((y > y_bag_low & y < y_bag_hi), 1, 0),
         mgp_inCI = ifelse((y > y_mgp_low & y < y_mgp_hi), 1, 0),
         bag_CIwidth = y_bag_hi - y_bag_low,
         mgp_CIwidth = y_mgp_hi - y_mgp_low)

## prediction metric
res[n_tr+1:n_tt,] %>%
  summarise(bag_rmspe = sqrt(mean((y-y_bag)^2)),
            mgp_rmspe = sqrt(mean((y-y_mgp)^2)),
            bag_mape = mean(abs(y-y_bag)),
            mgp_mape = mean(abs(y-y_mgp)),
            bag_coverage = mean(bag_inCI),
            mgp_coverage = mean(mgp_inCI),
            bag_meanwidth = mean(bag_CIwidth),
            mgp_meanwidth = mean(mgp_CIwidth)) %>%
  round(3)

#############
# Figure S3 #
#############
# infer wind direction
z_postm <- apply(bagres$z_save, 1, getmode)
dir_mean <- data.frame(
  SW = apply(bagres$z_save, 1, function(x) sum(x == "SW")/mcmc$save),
  W = apply(bagres$z_save, 1, function(x) sum(x == "W")/mcmc$save),
  NW = apply(bagres$z_save, 1, function(x) sum(x == "NW")/mcmc$save),
  N = apply(bagres$z_save, 1, function(x) sum(x == "N")/mcmc$save)
)

## wind with the highest posterior probability
### wind name
wind_hpp <- apply(dir_mean, 1, function(x) directions[which.max(x)])
### posterior prob
val_hpp <- apply(dir_mean, 1, max)
### weighted average of wind directions in angle: angles should be +/- around 0 (W)
ang_mean <- pi/4*dir_mean$SW + 0*dir_mean$W - pi/4*dir_mean$NW - pi/2*dir_mean$N

## by partitions
ptts_tr <- sort(unique(bagres$coords_ptt$partition[1:n_tr]))
time_name <- c(paste0("2020-02-04 ", 20:24, ":00"),
               paste0("2020-02-05 ", 1:7, ":00"))
bag_wind_df <- data.frame(partition = ptts_tr,
                          partition1 = ptts_tr,
                          ang_mean = ang_mean,
                          wind_hpp = factor(wind_hpp,
                                            levels = directions),
                          val_hpp = val_hpp) %>%
  tidyr::separate(partition1, c("row", "col", "time_d"),
                  sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row],
         time = time_name[time_d])

# overlay wind directions and predicted plots
plt_res <-  (res %>% select(easting, northing, time, y, y_bag, y_mgp)) %>%
  left_join(rbind(coords_tr[,c("easting", "northing", "time",
                               "true_easting", "true_northing")],
                  coords_tt[,c("easting", "northing", "time",
                               "true_easting", "true_northing")],
                  coords_grid[,c("easting", "northing", "time",
                                 "true_easting", "true_northing")]),
            by = c("easting", "northing", "time"))
plt_res <- plt_res %>%
  mutate(pm_bag = ifelse(exp(y_bag)<=15,1,
                         ifelse(exp(y_bag)<=35,2,
                                ifelse(exp(y_bag)<=75,3,4))))

plt1list = plt2list = plt3list <- list()
ttselect <- c(2, 4, 8, 10)
plt_res1 <- plt_res %>% filter(time %in% tgrid[ttselect])
plt_ymin1 <- min(c(plt_res1$y, plt_res1$y_bag), na.rm = TRUE)
plt_ymax1 <- max(c(plt_res1$y, plt_res1$y_bag), na.rm = TRUE)
for (ttt in 1:length(ttselect)) {
  tt <- ttselect[ttt]
  plt1list[[ttt]] <- plt_res %>% filter(time == tgrid[tt]) %>%
    ggplot() + geom_point(aes(true_easting, true_northing, color = y)) +
    scale_color_distiller(palette = "Spectral", direction = -1,
                          limits = c(plt_ymin1, plt_ymax1),
                          na.value = "transparent") +
    geom_sf(data = KR_sf, fill = "NA") +
    labs(y="", x="",
         title = time_name[tt], col = "log(PM2.5)")

  plt2list[[ttt]] <- plt_res %>% filter(time == tgrid[tt]) %>%
    ggplot() + geom_point(aes(true_easting, true_northing, color = y_bag)) +
    scale_color_distiller(palette = "Spectral", direction = -1,
                          limits = c(plt_ymin1, plt_ymax1),
                          na.value = "transparent") +
    geom_sf(data = KR_sf, fill = "NA") +
    labs(y="", x="",
         title = time_name[tt], col = "log(PM2.5)")

  plt3list[[ttt]] <- plt_res[(n_tr+n_tt)+1:n_grid,] %>%
    filter(time == tgrid[tt]) %>%
    ggplot() + geom_tile(aes(true_easting, true_northing, fill = factor(pm_bag))) +
    scale_fill_manual(
      values = c("1" = "#009E73", "2" = "#F0E442", "3" = "#D55E00"),
      labels = c("1" = "good\n(0~15)", "2" = "moderate\n(16~35)", "3" = "unhealthy\n(36~75)")) +
    geom_sf(data = KR_sf, fill = "NA") +
    geom_spoke(data = bag_wind_df %>% filter(time_d == tt),
               aes(x = easting_center, y = northing_center, angle = ang_mean),
               radius = 20000,
               arrow = arrow(length = unit(.1, 'cm'))) +
    labs(y = "", x = "",
         title = time_name[tt], fill = "PM2.5\n(um/m3)") +
    theme(legend.text=element_text(size=13),
          legend.title=element_text(size=13))
}

g1 <- ggpubr::ggarrange(plt1list[[1]],
                        plt1list[[2]],
                        plt1list[[3]],
                        plt1list[[4]],
                        common.legend = TRUE, legend = "bottom", nrow = 4)

g2 <- ggpubr::ggarrange(plt2list[[1]],
                        plt2list[[2]],
                        plt2list[[3]],
                        plt2list[[4]],
                        common.legend = TRUE, legend = "bottom", nrow = 4)

g3 <- ggpubr::ggarrange(plt3list[[1]],
                        plt3list[[2]],
                        plt3list[[3]],
                        plt3list[[4]],
                        common.legend = TRUE, legend = "bottom", nrow = 4)

gall <- ggpubr::ggarrange(g1, g2, g3, ncol=3)
# for (ext in extension) {
#   ggsave(plot = gall, paste0(path, "plots/KR_hourly_pred_final", ext),
#          width = 14, height = 15)
# }

###############
# Convergence #
###############
# trace plot and running mean plot of tau_sq
tracep <- data.frame(tau_sq = bagres$tau_sq_save, iter = 1:mcmc$save) %>%
  ggplot() +
  geom_line(aes(iter, tau_sq)) +
  labs(y = expression(tau^2), x = "") +
  theme(plot.margin = margin(b = -2, r = -2, t = 0, l = 0))

runmeanp <- data.frame(psi = bagres$tau_sq_save/1:mcmc$save,
                       iter = 1:mcmc$save) %>%
  ggplot() +
  geom_line(aes(iter, psi)) +
  labs(y = expression(paste("Running mean of ", tau^2)), x = "") +
  theme(plot.margin = margin(b = 0, r = -2, t = -2, l = 0))

# ESS of w and y
essp <- data.frame(rbind(coords_tr[,c("true_easting", "true_northing", "time_d")],
                         coords_tt[,c("true_easting", "true_northing", "time_d")],
                         coords_grid[,c("true_easting", "true_northing", "time_d")]),
                   ess = bagres$y_ess) %>%
  filter(time_d == 5) %>%
  ggplot() +
  geom_point(aes(true_easting, true_northing, color = ess), size = 0.3) +
  geom_sf(data = KR_sf, fill = "NA") +
  labs(x = "", y = "", color = "ESS of y") +
  scale_color_scico(palette = "cork") +
  theme(legend.position = "bottom",
        plot.margin = margin(b = 0, r = -2, t = 0, l = -2),
        legend.margin = margin(b = 0, r = -2, t = -2, l = -2))

conv1 <- gridExtra::grid.arrange(tracep, runmeanp, nrow = 2)
conv2 <- gridExtra::grid.arrange(conv1, essp, ncol = 2)
# for (ext in extension) {
#   ggsave(plot = conv2, paste0(path, "plots/KR_hourly_convergence", ext),
#          width = 8, height = 4)
# }

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
