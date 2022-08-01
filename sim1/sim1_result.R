# Fitted G-BAG is correctly specified #
# 2 Tables, 3 Figures
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(meshed)
library(INLA)
library(bags)
library(scico)

# path
path <- "~/GBAGs/"

# plot extensions
extension <- c(".pdf", ".eps", ".png")

##########
# Result #
##########
# theta1
bag1 <- readRDS(paste0(path, "sim1/sim1a_bag.RDS"))
mgp1 <- readRDS(paste0(path, "sim1/sim1a_mgp.RDS"))
inlares1 <- readRDS(paste0(path, "sim1/sim1a_inla.RDS"))

# theta2
bag2 <- readRDS(paste0(path, "sim1/sim1b_bag.RDS"))
mgp2 <- readRDS(paste0(path, "sim1/sim1b_mgp.RDS"))
inlares2 <- readRDS(paste0(path, "sim1/sim1b_inla.RDS"))

# choose either theta1 or theta2
bag <- bag1
mgp <- mgp1
inlares <- inlares1
# bag <- bag2
# mgp <- mgp2
# inlares <- inlares2

#########
# Table #
#########
# mean
mtab <-
  rbind(cbind(sapply(bag[c("beta_hat", "tausq_hat", "sigsq_hat")], mean),
              sapply(mgp[c("beta_hat", "tausq_hat", "sigsq_hat")], mean),
              c(sapply(inlares[c("beta_hat", "tausq_hat")], mean), NA),
              c(sapply(inlares[c("beta_hat_ns", "tausq_hat_ns")], mean), NA)),
        cbind(rowMeans(bag[["psi_hat"]]),
              rowMeans(mgp[["psi_hat"]]),
              rep(NA, 3),
              rep(NA, 3)),
        cbind(sapply(bag[c("rmspe", "mape", "coverage", "meanwidth")], mean),
              sapply(mgp[c("rmspe", "mape", "coverage", "meanwidth")], mean),
              sapply(inlares[c("rmspe", "mape", "coverage", "meanwidth")], mean),
              sapply(inlares[c("rmspe_ns", "mape_ns", "coverage_ns", "meanwidth_ns")], mean)))
rownames(mtab)[c(1:6)] <- c("beta", "tausq", "sigsq", "a", "c", "kappa")
colnames(mtab) <- c("G-BAG", "Q-MGP", "SPDE-stationary", "SPDE-nonstationary")
mtab %>% round(3)

# sd
sdtab <-
  rbind(cbind(sapply(bag[c("beta_hat", "tausq_hat", "sigsq_hat")], sd),
              sapply(mgp[c("beta_hat", "tausq_hat", "sigsq_hat")], sd),
              c(sapply(inlares[c("beta_hat", "tausq_hat")], sd), NA),
              c(sapply(inlares[c("beta_hat_ns", "tausq_hat_ns")], sd), NA)),
        cbind(apply(bag[["psi_hat"]], 1, sd),
              apply(mgp[["psi_hat"]], 1, sd),
              rep(NA, 3),
              rep(NA, 3)),
        cbind(sapply(bag[c("rmspe", "mape", "coverage", "meanwidth")], sd),
              sapply(mgp[c("rmspe", "mape", "coverage", "meanwidth")], sd),
              sapply(inlares[c("rmspe", "mape", "coverage", "meanwidth")], sd),
              sapply(inlares[c("rmspe_ns", "mape_ns", "coverage_ns", "meanwidth_ns")], sd)))
rownames(sdtab)[c(1:6)] <- c("beta", "tausq", "sigsq", "a", "c", "kappa")
colnames(sdtab) <- c("G-BAG", "Q-MGP", "SPDE-stationary", "SPDE-nonstationary")
sdtab %>% round(3)

# 95% ci
citab <-
  rbind(cbind(sapply(bag[c("beta_ci", "tausq_ci", "sigsq_ci")], mean),
              sapply(mgp[c("beta_ci", "tausq_ci", "sigsq_ci")], mean),
              c(sapply(inlares[c("beta_ci", "tausq_ci")], mean), NA),
              c(sapply(inlares[c("beta_ci_ns", "tausq_ci_ns")], mean), NA)),
        cbind(rowMeans(bag[["psi_ci"]]),
              rowMeans(mgp[["psi_ci"]]),
              rep(NA, 3),
              rep(NA, 3)))
rownames(citab) <- c("beta", "tausq", "sigsq", "a", "c", "kappa")
colnames(citab) <- c("G-BAG", "Q-MGP", "SPDE-stationary", "SPDE-nonstationary")
citab %>% round(3)

###############
# Directions  #
###############
## directions
directions <- c("W", "NW", "N", "NE")

## mcmc
mcmc <- list(save = 1000, burn = 10000, thin = 7)

## grid
n_time <- 8
tgrid <- seq(0, 1, length = n_time)

## posterior probability
dir_mean <- data.frame(
  W = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "W")/mcmc$save))),
  NW = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "NW")/mcmc$save))),
  N = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "N")/mcmc$save))),
  NE = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "NE")/mcmc$save)))
)

## wind with the highest posterior probability
wind_hpp <- apply(dir_mean, 1, function(x) directions[which.max(x)])
val_hpp <- apply(dir_mean, 1, max)

## irregular partitions
n_easting <- 6
n_northing <- 6
northing_sint <- seq(0, 1, length = (n_northing-2)^2+1)
northing_Lint <- seq(0, 1, length = (n_northing-2)+1)
breaks_northing <- c(northing_sint[1:2], northing_Lint[2:(n_northing-2)],
                     northing_sint[(n_northing-2)^2+0:1])
easting_sint <- seq(0, 1, length = (n_easting-2)^2+1)
easting_Lint <- seq(0, 1, length = (n_easting-2)+1)
breaks_easting <- c(easting_sint[1:2], easting_Lint[2:(n_easting-2)],
                    easting_sint[(n_easting-2)^2+0:1])

easting_width <- breaks_easting[-1] - breaks_easting[1:n_easting]
northing_height <- breaks_northing[-1] - breaks_northing[1:n_northing]
easting_midval <- 0.5*breaks_easting[1:n_easting] + 0.5*breaks_easting[-1]
northing_midval <- 0.5*breaks_northing[1:n_northing] + 0.5*breaks_northing[-1]

## train data
prob <- 0.8
n_tr <- nrow(bag$out$coords_ptt)*prob
ptts_tr <- sort(unique(bag$out$coords_ptt$partition[1:n_tr]))

## true wind directions
z_true <- rep("W", length(ptts_tr))
ptt_NW <- grep("\\d+,4,2|\\d+,5,2|\\d+,6,2|\\d+,\\d+,3|1,\\d+,4|2,\\d+,4|3,\\d+,4|4,2,8|4,3,8|4,4,8|4,5,8|4,6,8|5,6,8|6,6,8", ptts_tr)
z_true[ptt_NW] <- "NW"
ptt_N <- grep("4,\\d+,4|5,\\d+,4|6,\\d+,4|\\d+,\\d+,5|1,\\d+,6|2,\\d+,6|3,\\d+,6|3,1,8|3,2,8|3,3,8|3,4,8|3,5,8|4,1,8", ptts_tr)
z_true[ptt_N] <- "N"
ptt_NE <- grep("4,\\d+,6|5,\\d+,6|6,\\d+,6|\\d+,\\d+,7|5,1,8|5,2,8|5,3,8|5,4,8|5,5,8|6,1,8|6,2,8|6,3,8|6,4,8|6,5,8", ptts_tr)
z_true[ptt_NE] <- "NE"

bag_wind_df <- data.frame(partition = ptts_tr,
                          partition1 = ptts_tr,
                          wind = factor(z_true, levels = directions),
                          wind_hpp = factor(wind_hpp,
                                            levels = directions),
                          val_hpp = val_hpp) %>%
  tidyr::separate(partition1, c("row", "col", "time_d"),
                  sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row],
         width = easting_width[col],
         height = northing_height[n_northing + 1 - row],
         time = tgrid[time_d])

bag_wind_df$timename <- paste0("time = ", round(bag_wind_df$time,3))

g <- bag_wind_df %>%
  ggplot() +
  geom_tile(aes(x = easting_center, y = northing_center,
                width = width, height = height, fill = wind)) +
  facet_wrap(~timename, nrow = 2) +
  labs(x = "", y = "", fill = "Direction") +
  scale_fill_scico_d(palette = 'vik') +
  theme(legend.position = "none")
# for (ext in extension) {
#   ggsave(plot = g, paste0(path, "plots/sim1_directions", ext),
#          width = 9, height = 5)
# }

gg <- bag_wind_df %>% ggplot() +
  geom_tile(aes(x = easting_center, y = northing_center,
                width = width, height = height,
                fill = wind_hpp, alpha = val_hpp)) +
  labs(x = "", y = "", alpha = "Posterior probability", fill = "Direction") +
  scale_alpha(range = c(0.1, 1.5)) +
  facet_wrap(~timename, nrow = 2) +
  scale_fill_scico_d(palette = 'vik') +
  theme(legend.position = "bottom",
        legend.margin = margin(b = 0, r = 0, t = -15, l = 0))
# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/sim1a_windfill", ext),
#          width = 9, height = 5.2)
# }

###############
# Convergence #
###############
# chi-squared test for z
mean(sapply(bag1$z_pvalue, function(x) mean(x < 0.05)))
mean(sapply(bag2$z_pvalue, function(x) mean(x < 0.05)))

# effective sample size
## predicted y
summary(unlist(bag1$y_pred_ess))
summary(unlist(bag2$y_pred_ess))

## parameters
esstab <- cbind(c(sapply(bag1[c("beta_ess", "tausq_ess", "sigsqc_ess")], mean),
                  rowMeans(bag1[["psi_ess"]])[c(1,3)]),
                c(sapply(bag2[c("beta_ess", "tausq_ess", "sigsqc_ess")], mean),
                  rowMeans(bag2[["psi_ess"]])[c(1,3)]))
rownames(esstab) <- c("beta", "tausq", "sigsqc", "a", "kappa")
colnames(esstab) <- c("ESS - theta1", "ESS - theta2")
esstab %>% round(3)

# trace plots and running mean
convp <- data.frame(y = c(cumsum(bag1$out$y_pred_save[1,])/1:mcmc$save,
                          bag2$out$y_pred_save[1,]),
                    tau_sq = c(cumsum(bag1$out$tau_sq_save)/1:mcmc$save,
                               bag2$out$tau_sq_save),
                    beta = c(cumsum(bag1$out$beta_save[1,])/1:mcmc$save,
                             bag2$out$beta_save[1,]),
                    theta = rep(c("theta1", "theta2"), each = mcmc$save),
                    iter = rep(1:mcmc$save, 2)) %>%
  pivot_longer(-c(theta, iter), values_to = "draw", names_to = "param") %>%
  ggplot() +
  geom_line(aes(iter, draw)) +
  facet_grid(factor(param, labels = c(bquote(beta), bquote(tau^2), "y(t)")) ~
               factor(theta, labels = c(bquote(theta[1]), bquote(theta[2]))),
             scales = "free_y", labeller = label_parsed) +
  labs(x = "", y = "") +
  theme(plot.margin = margin(t = 0, l = -5, r = 0, b = -5))
# for (ext in extension) {
#   ggsave(plot = convp, paste0(path, "plots/sim1_convergence", ext),
#          width = 8, height = 4)
# }

################
# Running time #
################
mean(bag1$est_time_per_iter[3,] + bag1$pred_time_per_iter[3,])
mean(bag1$est_time_per_iter[3,])
mean(bag1$pred_time_per_iter[3,])

mean(bag2$est_time_per_iter[3,] + bag2$pred_time_per_iter[3,])
mean(bag2$est_time_per_iter[3,])
mean(bag2$pred_time_per_iter[3,])
