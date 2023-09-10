# Air Quality Analysis in CA #
# 1 Table, 6 Figures
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

n_easting <- 16 # (~ 55km, ~0.5 degree resolution)
n_northing <- 20
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
# saveRDS(out, paste0(path, "application/CA_daily_bag_all.RDS"))

# ## mgp ##
# set.seed(123)
# meshout <- spmeshed(family = c("gaussian"),
#                     y = c(y_tr, rep(NA, n_tt + n_grid)), x = x_scale,
#                     coords = as.matrix(
#                       rbind(coords_tr[,c("easting", "northing", "time")],
#                             coords_tt[,c("easting", "northing", "time")],
#                             coords_grid[,c("easting", "northing", "time")])),
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
#                     n_samples = mcmc$save,
#                     n_burn = mcmc$burn,
#                     n_thin = mcmc$thin,
#                     verbose = 100,
#                     n_threads = 10)
#
# saveRDS(meshout, paste0(path, "application/CA_daily_mgp_all.RDS"))

# ## INLA ##
# set.seed(123)
# ybar <- mean(y_tr)
#
# # mesh
# diff_east <- diff(range(coords_tr$easting))
# max_edge <- diff_east/15
# uniq_coords_tr <- unique(coords_tr[,1:2])
# prmesh1 <- inla.mesh.create.helper(points = uniq_coords_tr,
#                                    max.edge = c(0.1,5)*max_edge,
#                                    cutoff = max_edge/2,
#                                    offset = c(max_edge, diff_east/7))
#
# # nonstationary
# # range vary by coordinates of mesh nodes
# maxx <- max(prmesh1$loc[,1])
# minx <- min(prmesh1$loc[,1])
# maxy <- max(prmesh1$loc[,2])
# miny <- min(prmesh1$loc[,2])
# nu <- 1
# alpha <- nu + 2 / 2
# logkappa0 <- log(8 * nu) / 2                                                    # log(kappa)
# logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2                    # log(tau)
# logtau0 <- logtau0 - logkappa0
# spde_ns <- inla.spde2.matern(prmesh1,
#                              B.tau = cbind(logtau0, -1, 1,
#                                            (prmesh1$loc[,1]-minx)/(maxx-minx),
#                                            (prmesh1$loc[,2]-miny)/(maxy-miny)),
#                              B.kappa = cbind(logkappa0, 0, -1,
#                                              -(prmesh1$loc[,1]-minx)/(maxx-minx),
#                                              -(prmesh1$loc[,2]-miny)/(maxy-miny)),
#                              theta.prior.mean = rep(0, 4),
#                              theta.prior.prec = rep(1, 4))
#
# # for space-TIME model
# iset_ns <- inla.spde.make.index('i', n.spde = spde_ns$n.spde, n.group = n_time)
#
# # connect observations to mesh nodes
# A.obs <- inla.spde.make.A(mesh = prmesh1,
#                           loc = as.matrix(coords_tr[,1:2]),
#                           group = as.numeric(
#                             as.character(coords_tr$time*(n_time-1) + 1)))       # group should be an index
# stk.obs_ns <- inla.stack(
#   data = list(y = y_tr - ybar),
#   A = list(A.obs, 1),
#   effects = list(iset_ns, X = X_tr),
#   tag = 'obs')
#
# # same for prediction
# proj.pred <- inla.mesh.projector(prmesh1,
#                                  loc = as.matrix(
#                                    rbind(coords_tt[,c("easting", "northing")],
#                                          coords_grid[,c("easting", "northing")])))
# A.pred <- inla.spde.make.A(prmesh1, loc = proj.pred$loc,
#                            group = as.numeric(as.character(
#                              c(coords_tt$time, coords_grid$time)*
#                                (n_time-1) + 1)))
# stk.pred_ns <- inla.stack(data = list(y = NA),
#                           A = list(A.pred, 1),
#                           effects = list(iset_ns, X = X_pred),
#                           tag = "pred")
# stk_ns <- inla.stack(stk.obs_ns, stk.pred_ns)
#
# # model formula
# h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0.05, 0.99)))            # P(temporal auto-corr > 0.05) = 0.99 (not a problem)
# formulae_ns <- y ~ -1 + X + f(i, model = spde_ns, group = i.group,
#                               control.group =
#                                 list(model = 'ar1', hyper = h.spec))
#
# # PC prior on the autoregressive parameter
# prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))                       # P(precision of AR > 1) = 0.01 (not a problem)
# time_inla_ns <- system.time({
#   inlaout_ns <- inla(formulae_ns, data = inla.stack.data(stk_ns),
#                      control.predictor = list(compute = TRUE,
#                                               A = inla.stack.A(stk_ns)),
#                      control.family = list(hyper = list(prec = prec.prior)),
#                      num.threads = 10,
#                      verbose = TRUE)
# })
#
# index.pred <- c(inla.stack.index(stk_ns, "obs")$data,
#                 inla.stack.index(stk_ns, "pred")$data)
# y_inla <- ybar + inlaout_ns$summary.fitted.values[index.pred, "mean"]           # include Xb
# tausq_hat_inla <-
#   mean(1/inlaout_ns$marginals.hyperpar$
#          `Precision for the Gaussian observations`[,1])
# inla_sd <- sqrt(inlaout_ns$summary.fitted.values[index.pred, "sd"]^2 +          # variability from Xb fitting
#                 tausq_hat_inla)
# y_inla_low <- y_inla + qnorm(.025)*inla_sd
# y_inla_hi <- y_inla + qnorm(.975)*inla_sd
#
# saveRDS(list(time_inla_ns = time_inla_ns,
#              inlaout_ns = list(summary.fixed =
#                                  inlaout_ns$summary.fixed,
#                                marginals.hyperpar =
#                                  inlaout_ns$marginals.hyperpar,
#                                summary.fitted.values =
#                                  inlaout_ns$summary.fitted.values,
#                                summary.random =
#                                  inlaout_ns$summary.random),
#              y_inla = y_inla,
#              y_inla_low = y_inla_low,
#              y_inla_hi = y_inla_hi,
#              prmesh1 = prmesh1,
#              iset_ns = iset_ns),
#         paste0(path, "application/CA_daily_inla.RDS"))

# # split CA_daily_*_all.RDS due to size
# ## bags ##
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
#         paste0(path, "application/CA_daily_bag.RDS"))
#
# ## mgp ##
# # prediction
# y_mgp <- meshout$yhat_mcmc %>% summary_list_mean()
# y_mgp_low <- meshout$yhat_mcmc %>% summary_list_q(q = 0.025)
# y_mgp_hi <- meshout$yhat_mcmc %>% summary_list_q(q = 0.975)
#
# saveRDS(list(tausq_mcmc = meshout$tausq_mcmc[,mcmc$thin*(1:mcmc$save)],
#              beta_mcmc = meshout$beta_mcmc[,,mcmc$thin*(1:mcmc$save)],
#              sigsq_mcmc = meshout$theta_mcmc[4,,mcmc$thin*(1:mcmc$save)],
#              psi_mcmc = meshout$theta_mcmc[-4,,mcmc$thin*(1:mcmc$save)],        # temporal decay, spatial decay, interaction
#              time_per_iter = meshout$mcmc_time/meshout$mcmc,
#              coordsdata = meshout$coordsdata,
#              y_mgp = y_mgp,
#              y_mgp_low = y_mgp_low,
#              y_mgp_hi = y_mgp_hi),
#         paste0(path, "application/CA_daily_mgp.RDS"))

###########
# Results #
###########
# call results
bagres <- readRDS(paste0(path, "application/CA_daily_bag.RDS"))
mgpres <- readRDS(paste0(path, "application/CA_daily_mgp.RDS"))
inlares <- readRDS(paste0(path, "application/CA_daily_inla.RDS"))

# time per iteration in seconds
bagres$est_time_per_iter
bagres$pred_time_per_iter

#########
# Table #
#########
# parameter estimation
## beta
rowMeans(bagres$beta_save) %>% round(3)
quantile(bagres$beta_save, probs = c(0.025, 0.975)) %>% round(3)

mean(mgpres$beta_mcmc) %>% round(3)
quantile(mgpres$beta_mcmc, probs = c(0.025, 0.975)) %>% round(3)

(beta_hat_inla <- as.numeric(inlares$inlaout_ns$summary.fixed["mean"])) %>% round(3)
as.numeric(inlares$inlaout_ns$summary.fixed[c("0.025quant","0.975quant")]) %>% round(3)

## tausq
mean(bagres$tau_sq_save) %>% round(3)
quantile(bagres$tau_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(mgpres$tausq_mcmc) %>% round(3)
quantile(mgpres$tausq_mcmc, probs = c(0.025, 0.975)) %>% round(3)

mean(1/inlares$inlaout_ns$marginals.hyperpar$
       `Precision for the Gaussian observations`[,1]) %>%
  round(3)
quantile(1/inlares$inlaout_ns$marginals.hyperpar$
           `Precision for the Gaussian observations`[,1],
         probs = c(0.025, 0.975)) %>% round(3)

## sigsq
mean(bagres$sig_sq_save) %>% round(3)
quantile(bagres$sig_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(mgpres$sigsq_mcmc) %>% round(3);
quantile(mgpres$sigsq_mcmc, probs = c(0.025, 0.975)) %>% round(3)

## psi
rowMeans(bagres$psi_save) %>% round(3)
apply(bagres$psi_save, 1, function(x) quantile(x, probs = c(0.025, 0.975))) %>%
  round(3)

rowMeans(mgpres$psi_mcmc) %>% round(3)
apply(mgpres$psi_mcmc, 1, function(x) quantile(x, probs = c(0.025, 0.975))) %>%
  round(3)

# prediction summaries
res <- data.frame(bagres$coords_ptt[,c("easting", "northing", "time")],
                  y = c(y_tr, y_tt, rep(NA, n_grid)),
                  y_bag = c(bagres$yhat, bagres$y_bag),
                  y_bag_low = c(bagres$yhat_low, bagres$y_bag_low),
                  y_bag_hi = c(bagres$yhat_hi, bagres$y_bag_hi),
                  y_inla = inlares$y_inla,
                  y_inla_low = inlares$y_inla_low,
                  y_inla_hi = inlares$y_inla_hi) %>%
  left_join(data.frame(mgpres$coordsdata,
                       y_mgp = mgpres$y_mgp,
                       y_mgp_low = mgpres$y_mgp_low,
                       y_mgp_hi = mgpres$y_mgp_hi) %>%
              select(-forced_grid), by = c("easting", "northing", "time")) %>%
  mutate(bag_inCI = ifelse((y > y_bag_low & y < y_bag_hi), 1, 0),
         mgp_inCI = ifelse((y > y_mgp_low & y < y_mgp_hi), 1, 0),
         inla_inCI = ifelse((y > y_inla_low & y < y_inla_hi), 1, 0),
         bag_CIwidth = y_bag_hi - y_bag_low,
         mgp_CIwidth = y_mgp_hi - y_mgp_low,
         inla_CIwidth = y_inla_hi - y_inla_low)

## prediction metric
ressum <- res[n_tr+1:n_tt,] %>%
  summarise(bag_rmspe = sqrt(mean((y-y_bag)^2)),
            mgp_rmspe = sqrt(mean((y-y_mgp)^2)),
            inla_rmspe = sqrt(mean((y-y_inla)^2)),
            bag_mape = mean(abs(y-y_bag)),
            mgp_mape = mean(abs(y-y_mgp)),
            inla_mape = mean(abs(y-y_inla)),
            bag_coverage = mean(bag_inCI),
            mgp_coverage = mean(mgp_inCI),
            inla_coverage = mean(inla_inCI),
            bag_meanwidth = mean(bag_CIwidth),
            mgp_meanwidth = mean(mgp_CIwidth),
            inla_meanwidth = mean(inla_CIwidth))

ressum_tab <- ressum %>% as.numeric() %>% matrix(., nrow = 4, ncol = 3, byrow = T)
colnames(ressum_tab) <- c("G-BAG", "Fixed DAG", "SPDE-nonstationary")
rownames(ressum_tab) <- c("RMSPE", "MAPE", "95% CI coverage", "Mean 95% CI width")
ressum_tab %>% round(3)

######################
# Prediction by time #
######################
## prediction metric by time
res_bytime <- res[n_tr+1:n_tt,] %>%
  group_by(time) %>%
  summarise(bag_rmspe = sqrt(mean((y-y_bag)^2)),
            mgp_rmspe = sqrt(mean((y-y_mgp)^2)),
            inla_rmspe = sqrt(mean((y-y_inla)^2)),
            bag_mape = mean(abs(y-y_bag)),
            mgp_mape = mean(abs(y-y_mgp)),
            inla_mape = mean(abs(y-y_inla)),
            bag_coverage = mean(bag_inCI),
            mgp_coverage = mean(mgp_inCI),
            inla_coverage = mean(inla_inCI),
            bag_meanwidth = mean(bag_CIwidth),
            mgp_meanwidth = mean(mgp_CIwidth),
            inla_meanwidth = mean(inla_CIwidth))

# relative change from x1 to x2: (x2-x1)/x1
res_bytime %>%
  mutate(rmspe_reduction_from_mgp = (bag_rmspe - mgp_rmspe)/mgp_rmspe,
         rmspe_reduction_from_inla = (bag_rmspe - inla_rmspe)/inla_rmspe,
         mape_reduction_from_mgp = (bag_mape - mgp_mape)/mgp_mape,
         mape_reduction_from_inla = (bag_mape - inla_mape)/inla_mape) %>%
  select(rmspe_reduction_from_mgp, rmspe_reduction_from_inla,
         mape_reduction_from_mgp, mape_reduction_from_inla) %>%
  min()

res_bytime %>%
  mutate(rmspe_increase_to_mgp = (mgp_rmspe - bag_rmspe)/bag_rmspe,
         rmspe_increase_to_inla = (inla_rmspe - bag_rmspe)/bag_rmspe,
         mape_increase_to_mgp = (mgp_mape - bag_mape)/bag_mape,
         mape_increase_to_inla = (inla_mape - bag_mape)/bag_mape) %>%
  select(rmspe_increase_to_mgp, rmspe_increase_to_inla,
         mape_increase_to_mgp, mape_increase_to_inla) %>%
  max()

res_bytime2 <-
  data.frame(date = rep(timelist, 9),
             vals = c(res_bytime$bag_rmspe,
                      res_bytime$mgp_rmspe,
                      res_bytime$inla_rmspe,
                      res_bytime$bag_mape,
                      res_bytime$mgp_mape,
                      res_bytime$inla_mape,
                      res_bytime$bag_coverage,
                      res_bytime$mgp_coverage,
                      res_bytime$inla_coverage),
             model = rep(rep(c("G-BAG", "Fixed DAG", "SPDE-nonstationary"),
                             each = n_time), 3),
             what = rep(c("RMSPE", "MAPE", "95% CI coverage"), each = 3*n_time)) %>% 
  mutate(model = fct_relevel(model, c("G-BAG", "Fixed DAG", "SPDE-nonstationary")))

g <- res_bytime2 %>% ggplot() +
  geom_line(aes(date, vals, col = model, linetype = model), linewidth = 0.8) +
  facet_grid(factor(what, levels = c("MAPE", "RMSPE", "95% CI coverage"))~.,
             scales = "free") +
  labs(x = "", y = "", color = "", linetype = "") +
  scale_color_manual(values = c("G-BAG" = "#1400FFFF", "Fixed DAG" = "#C729D6FF",
                                "SPDE-nonstationary" = "#FF9C63FF")) +
  theme(legend.position = "right",
        plot.margin = margin(t = 2, l = -3, r = 2, b = -9))
# for (ext in extension) {
#   ggsave(plot = g, paste0(path, "plots/CA_daily_pred", ext),
#          width = 7, height = 4)
# }

#############
# Residuals #
#############
# examine residual plots y-w
ybar <- mean(y_tr)
e_bag <- c(y_tr-ybar, y_tt-ybar) - c(bagres$what, bagres$w_bag[1:n_tt])

uniq_coords <- unique(rbind(coords_tr[,1:2], coords_tt[,1:2]))
inla_field <- inlares$inlaout_ns$summary.random$i$mean
proj_data <- inla.mesh.projector(mesh = inlares$prmesh1,
                                 loc = as.matrix(uniq_coords))
field_proj <- list()
for (j in 1:n_time) {
  field_proj[[j]] = inla.mesh.project(proj_data,
                                      inla_field[inlares$iset_ns$i.group == j]) %>%
    as.numeric()
}
w_inla <- unlist(field_proj)
e_inla <- c(y_tr-ybar, y_tt-ybar) - w_inla

e_res <- data.frame(rbind(coords_tr[,c("easting", "northing", "time", "time_d",
                                       "true_easting", "true_northing")],
                          coords_tt[,c("easting", "northing", "time", "time_d",
                                       "true_easting", "true_northing")]),
                    e_bag = e_bag,
                    e_inla = e_inla,
                    y_inla = inlares$y_inla[1:n],
                    x = c(X_tr, X_pred[1:n_tt])) %>%
  mutate(date = timelist[time_d],
         xb_inla = x*beta_hat_inla)

ttselect <- c(21, 40, 72)
plt_eres1 <- e_res %>% filter(time %in% tgrid[ttselect])
plt_emin1 <- min(c(plt_eres1$e_bag, plt_eres1$e_inla), na.rm = TRUE)
plt_emax1 <- max(c(plt_eres1$e_bag, plt_eres1$e_inla), na.rm = TRUE)

## G-BAG
fig1 <- plt_eres1 %>%
  ggplot() + geom_point(aes(true_easting, true_northing, color = e_bag)) +
  scale_color_scico(palette = "roma", direction = -1,
                    limits = c(-6.2, 6.6),
                    na.value = "transparent") +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_wrap(~date, nrow = 1) +
  labs(y="", x="", col = "Residual   \n(G-BAG)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 2, l = -2, r = 2, b = -2),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6))

## SPDE
fig2 <- plt_eres1 %>%
  ggplot() +
  geom_point(aes(true_easting, true_northing, color = e_inla - xb_inla)) +
  scale_color_scico(palette = "roma", direction = -1,
                    limits = c(-6.2, 6.6),
                    na.value = "transparent") +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_wrap(~date, nrow = 1) +
  labs(y="", x="", col = "Residual   \n(SPDE)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 2, l = -2, r = 2, b = -2),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6))

fig3 <- plt_eres1 %>%
  ggplot() +
  geom_point(aes(true_easting, true_northing, color = y_inla)) +
  scale_color_distiller(palette = "Spectral", direction = -1,
                        na.value = "transparent") +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_wrap(~date, nrow = 1) +
  labs(y="", x="", col = "log(PM2.5)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 2, l = -2, r = 2, b = -2),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6))

fig <- gridExtra::grid.arrange(fig3, fig2, fig1, ncol = 1)
# for (ext in extension) {
#   ggsave(plot = fig, paste0(path, "plots/CA_daily_res", ext),
#          width = 7, height = 8)
# }

######################
# Predicted surfaces #
######################
# infer wind direction
z_postm <- apply(bagres$z_save, 1, getmode)
dir_mean <- data.frame(
  W = apply(bagres$z_save, 1, function(x) sum(x == "W")/mcmc$save),
  NW = apply(bagres$z_save, 1, function(x) sum(x == "NW")/mcmc$save),
  N = apply(bagres$z_save, 1, function(x) sum(x == "N")/mcmc$save),
  NE = apply(bagres$z_save, 1, function(x) sum(x == "NE")/mcmc$save)
)
table(z_postm)/length(z_postm)

## wind with the highest posterior probability
### wind name
wind_hpp <- apply(dir_mean, 1, function(x) directions[which.max(x)])
### posterior prob
val_hpp <- apply(dir_mean, 1, max)
### angle
ang_hpp <- rep(0, length(val_hpp))
ang_hpp[which(wind_hpp == "NW")] <- 7*pi/4
ang_hpp[which(wind_hpp == "N")] <- 3*pi/2
ang_hpp[which(wind_hpp == "NE")] <- 5*pi/4

## winds in partitioned regions
ptts_tr <- sort(unique(bagres$coords_ptt$partition[1:n_tr]))
bag_wind_df <- data.frame(partition = ptts_tr,
                           partition1 = ptts_tr,
                           val_hpp = val_hpp,
                           ang_hpp = ang_hpp) %>%
  tidyr::separate(partition1,
                  c("row", "col", "time_d"), sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row],
         time = timelist[time_d])

# overlay wind directions and predicted plots
## sanity check for left_join
all.equal(bagres$coords_ptt[1:n_tr, "easting"], 
          pull(coords_tr[,"easting"]))
all.equal(bagres$coords_ptt[1:n_tr, "northing"], 
          pull(coords_tr[,"northing"]))
all.equal(bagres$coords_ptt[1:n_tr, "time"], 
          pull(coords_tr[,"time"]))

all.equal(bagres$coords_ptt[n_tr + 1:n_tt, "easting"], 
          pull(coords_tt[,"easting"]))
all.equal(bagres$coords_ptt[n_tr + 1:n_tt, "northing"], 
          pull(coords_tt[,"northing"]))
all.equal(bagres$coords_ptt[n_tr + 1:n_tt, "time"], 
          pull(coords_tt[,"time"]))

all.equal(bagres$coords_ptt[(n_tr + n_tt) + 1:n_grid, "easting"], 
          coords_grid[,"easting"])
all.equal(bagres$coords_ptt[(n_tr + n_tt) + 1:n_grid, "northing"], 
          coords_grid[,"northing"])
all.equal(bagres$coords_ptt[(n_tr + n_tt) + 1:n_grid, "time"], 
          coords_grid[,"time"])
coords_grid[,"northing"] <- bagres$coords_ptt[(n_tr + n_tt) + 1:n_grid, "northing"]

plt_res <- (res %>% select(easting, northing, time,
                           y, y_bag, y_mgp, y_inla, bag_CIwidth)) %>%
  left_join(data.frame(rbind(coords_tr[,c("easting", "northing", "time",
                                          "time_d", "true_easting", "true_northing")],
                             coords_tt[,c("easting", "northing", "time",
                                          "time_d", "true_easting", "true_northing")],
                             coords_grid[,c("easting", "northing", "time",
                                            "time_d", "true_easting", "true_northing")])),
            by = c("easting", "northing", "time")) %>%
  mutate(date = timelist[time_d])

# Sep 13 - 18
ttselect <- c(44:45, 47:49)
plt_res4 <- plt_res[(n_tr + n_tt) + 1:n_grid,] %>%
  filter(time %in% tgrid[ttselect])
plt_res4_wind <- bag_wind_df %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)

## simply to make green area visible
plt_ymin4 <- max(0,min(plt_res4$y_bag))
plt_ymax4 <- max(plt_res4$y_bag)

plt_res_tmp <- plt_res %>% filter(time == tgrid[1])
tmp_nn <- RANN::nn2(data = expand.grid(easting = easting_cut/1000,
                                       northing = northing_cut/1000),
          query = plt_res_tmp[,c("easting", "northing")], k = 1)
plt_res5 <- plt_res %>%
  mutate(dist_to_bound = rep(tmp_nn$nn.dists, times = n_time)) %>%
  filter(time %in% tgrid[ttselect])
plt_pred <- plt_res5 %>%
  ggplot() + geom_point(aes(true_easting, true_northing, color = y_bag,
                            alpha = 1/dist_to_bound), size = 2) +
  scale_color_distiller(palette = "Spectral", direction = -1,
                        limits = c(plt_ymin4, plt_ymax4),
                        na.value = "transparent") +
  scale_alpha_binned(n.breaks = 3, range=c(.1,1), guide = "none") +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_wrap(~date, nrow = 1) +
  labs(y="", x="", color = "log(PM2.5)") +
  theme(legend.position = "right",
        plot.margin=margin(t = 2, l = -2, r = 2, b = -5),
        legend.margin=margin(b = 0, r = 0, t = 0, l = -6)) +
  scale_x_continuous(breaks = c(-123, -121, -119, -117, -115),
                     labels = expression(123 * degree, 121 * degree,
                                         119 * degree, 117 * degree,
                                         115 * degree * W))

plt_pred_dis <- plt_res4 %>%
  ggplot() +
  geom_contour_filled(aes(true_easting, true_northing, z = exp(y_bag)),
                      breaks = c(Inf, 65, 35, 12, 0)) +
  scale_fill_manual(values = c("#D55E00","#F0E442","#009E73", "#56B4E9"),
                    labels = c(">65","(35,65]","(12,35]","(0,12]")) +
  geom_sf(data = ca_sf, fill = "NA") +
  geom_spoke(data = plt_res4_wind,
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             color = "white",
             radius = 40000, linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res4_wind,
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  facet_wrap(~date, nrow=1) +
  labs(y = "", x = "", fill = "PM2.5\n(um/m3)", alpha = "Posterior\nprobability") +
  theme(legend.position = "right",
        plot.margin = margin(t = 2,l = -2, r = 2, b = -5),
        legend.margin = margin(b = 0,r = 0,t = 0,l = -6)) +
  scale_x_continuous(breaks = c(-123, -121, -119, -117, -115),
                     labels = expression(123 * degree, 121 * degree,
                                         119 * degree, 117 * degree,
                                         115 * degree * W))

plt_obs <- plt_res %>%
  filter(time %in% tgrid[ttselect]) %>%
  ggplot() + geom_point(aes(true_easting, true_northing, color = y)) +
  scale_color_distiller(palette = "Spectral", direction = -1,
                        limits = c(plt_ymin4, plt_ymax4),
                        na.value = "transparent") +
  geom_sf(data = ca_sf, fill = "NA") +
  theme(legend.position = "right",
        plot.margin = margin(t = 2, l = -2, r = 2, b = -5),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6)) +
  facet_wrap(~date, nrow = 1) +
  labs(y = "", x = "", col = "log(PM2.5)") +
  scale_x_continuous(breaks = c(-123, -121, -119, -117, -115),
                     labels = expression(123 * degree, 121 * degree,
                                         119 * degree, 117 * degree,
                                         115 * degree * W))

gg <- ggpubr::ggarrange(plt_obs,
                        plt_pred,
                        plt_pred_dis,
                        common.legend = FALSE, legend = "right", nrow = 3)
# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/CA_daily_pred_final", ext),
#          width = 9, height = 7)
# }

###############
# Uncertainty #
###############
plotCI <- plt_res4 %>%
  filter(date == "2020-09-13") %>%
  ggplot() +
  geom_contour_filled(aes(true_easting, true_northing, z = bag_CIwidth),
                      bins = 4) +
  scale_fill_scico_d(palette = 'nuuk') +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_wrap(~date, nrow = 1) +
  labs(y = "", x = "", fill = "95% CI width")
# for (ext in extension) {
#   ggsave(plot = plotCI, paste0(path, "plots/CA_daily_pred_CI", ext),
#          width = 5, height = 4)
# }

##################
# August Complex #
##################
# August Complex
ttselect <- c(17, 18, 20, 22)
plt_res6 <- plt_res[(n_tr + n_tt) + 1:n_grid,] %>%
  filter(time %in% tgrid[ttselect])
plt_res6_wind <- bag_wind_df %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)

ggg <- plt_res6 %>% 
  mutate(z = exp(y_bag)) %>%
  ggplot() +
  geom_contour_filled(aes(true_easting, true_northing, z = z),
                      breaks = c(0, 12, 35, 65, Inf)) +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442", "#D55E00"),
                    labels = c("(0,12]","(12,35]","(35,65]",">65")) +
  geom_sf(data = ca_sf, fill = "NA") +
  geom_spoke(data = plt_res6_wind,
            aes(x = easting_center, y = northing_center, angle = ang_hpp),
            radius = 40000, col = "white", linewidth = 1.5,
            arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res6_wind,
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  facet_wrap(~ date, nrow = 1) +
  labs(y = "", x = "", fill = "PM2.5 (um/m3)") +
  scale_x_continuous(breaks = c(-123, -121, -119)) +
  coord_sf(xlim = c(-135000, 450000), ylim = c(4100000, 4700000), expand = FALSE) +
  theme(legend.position = "bottom",
        legend.margin = margin(b = -5, r = 0, t = -25, l = 0))

for (ext in extension) {
  ggsave(plot = ggg, paste0(path, "plots/CA_daily_pred_AugustComplex", ext),
         width = 7, height = 2.5)
}

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
# for (ext in extension) {
#   ggsave(plot = conv2, paste0(path, "plots/CA_daily_convergence", ext),
#          width = 6.8, height = 4)
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

