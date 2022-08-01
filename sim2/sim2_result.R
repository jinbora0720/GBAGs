# Fitted G-BAG is misspecified #
# 3 Figures
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
# theta3
bag1 <- readRDS(paste0(path, "sim2/sim2a_bag.RDS"))
mgp1 <- readRDS(paste0(path, "sim2/sim2a_mgp.RDS"))
inlares1 <- readRDS(paste0(path, "sim2/sim2a_inla.RDS"))

# theta4
bag2 <- readRDS(paste0(path, "sim2/sim2b_bag.RDS"))
mgp2 <- readRDS(paste0(path, "sim2/sim2b_mgp.RDS"))
inlares2 <- readRDS(paste0(path, "sim2/sim2b_inla.RDS"))

# choose either theta3 or theta4
bag <- bag1
mgp <- mgp1
inlares <- inlares1
# bag <- bag2
# mgp <- mgp2
# inlares <- inlares2

#####
# w #
#####
# examples of true surface
n_time0 <- 59
tgrid0 <- seq(0, 1, length = n_time0)

## with theta3
bag1$data0$timename <- paste0("time = ", round(bag1$data0$time,3))
plot1 <- bag1$data0 %>%
  filter(time > tgrid0[5] & time <= tgrid0[10]) %>%
  ggplot() +
  geom_raster(aes(easting, northing, fill = w)) +
  geom_contour(aes(easting, northing, z = w), breaks = 1.7, color = "black") +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  scale_y_continuous(breaks = c(0, .5, 1)) +
  scale_fill_viridis_c() +
  theme(plot.margin = margin(t = 2, l = -2, r = 2, b = -5),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6)) +
  facet_wrap(~timename, nrow = 1) +
  labs(x = "", y = "")

## with theta4
bag2$data0$timename <- paste0("time = ", round(bag2$data0$time,3))
plot2 <- bag2$data0 %>%
  filter(time > tgrid0[40] & time <= tgrid0[45]) %>%
  ggplot() +
  geom_raster(aes(easting, northing, fill = w)) +
  geom_contour(aes(easting, northing, z = w), breaks = 1.7, color = "black") +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  scale_y_continuous(breaks = c(0, .5, 1)) +
  scale_fill_viridis_c() +
  theme(plot.margin = margin(t = 2, l = -2, r = 2, b = -5),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6)) +
  facet_wrap(~timename, nrow = 1) +
  labs(x = "", y = "")

g <- gridExtra::grid.arrange(plot1, plot2, ncol = 1)
# for (ext in extension) {
#   ggsave(plot = g, paste0(path, "plots/sim2_truew", ext),
#          width = 10, height = 4.5)
# }

##############
# Directions #
##############
## directions
directions <- c("NW", "N", "NE", "E")

## assign partitions
n_easting <- 2
n_northing <- 6

## mcmc
mcmc <- list(save = 1000, burn = 5000, thin = 2)

## grid
n_time <- 30
tgrid <- seq(0, 1, length = n_time)

## posterior probability
dir_mean <- data.frame(
  NW = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "NW")/mcmc$save))),
  N = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "N")/mcmc$save))),
  NE = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "NE")/mcmc$save))),
  E = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "E")/mcmc$save)))
)

## weighted average of wind directions in angle
ang_mean <- 7/4*pi*dir_mean$NW + 3/2*pi*dir_mean$N +
  5/4*pi*dir_mean$NE + pi*dir_mean$E

## wind with the highest posterior probability
wind_hpp <- apply(dir_mean, 1, function(x) directions[which.max(x)])
val_hpp <- apply(dir_mean, 1, max)

## partitions
easting_cut <- seq(0, 1, length = n_easting + 1)
northing_cut <- seq(0, 1, length = n_northing + 1)
easting_midval <- 0.5*easting_cut[1:n_easting] + 0.5*easting_cut[-1]
northing_midval <- 0.5*northing_cut[1:n_northing] + 0.5*northing_cut[-1]

## train data
prob <- 0.8
n_tr <- nrow(bag$out$coords_ptt)*prob
ptts_tr <- sort(unique(bag$out$coords_ptt$partition[1:n_tr]))

bag_wind_df <- data.frame(partition = ptts_tr,
                          partition1 = ptts_tr,
                          ang_mean = ang_mean,
                          wind_hpp = factor(wind_hpp,
                                            levels = directions),
                          val_hpp = val_hpp) %>%
  tidyr::separate(partition1, c("row", "col", "time_d"),
                  sep = ",", convert = TRUE) %>%
  mutate(ang_hpp = ifelse(wind_hpp == "NW", 7/4*pi,
                          ifelse(wind_hpp == "N", 3/2*pi,
                                 ifelse(wind_hpp == "NE", 5/4*pi, pi))),
         easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row],
         time = tgrid[time_d])

# across time
bag_bypart <- data.frame(partition = ptts_tr, dir_mean) %>%
  tidyr::separate(partition, c("row", "col", "time_d"),
                  sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row]) %>%
  group_by(row, col, easting_center, northing_center) %>%
  summarise(NW = mean(NW), N = mean(N), NE = mean(NE), E = mean(E))

## plot: all arrows with prob. in each partition across times
fig1 <- bag_bypart[,directions] %>%
  reshape2::melt() %>%
  mutate(row = rep(as.numeric(unlist(bag_bypart[,1])), times = 4),
         col = rep(as.numeric(unlist(bag_bypart[,2])), times = 4),
         easting_center = rep(as.numeric(unlist(bag_bypart[,3])), times = 4),
         northing_center = rep(as.numeric(unlist(bag_bypart[,4])), times = 4),
         angle = ifelse(variable == "NW", 7/4*pi,
                        ifelse(variable == "N", 3/2*pi,
                               ifelse(variable == "NE", 5/4*pi, pi)))) %>%
  ggplot() +
  geom_spoke(aes(x = easting_center, y = northing_center+0.08, angle = angle,
                 radius = scales::rescale(value, c(.05, .15)),
                 size = value,
                 col = variable),
             arrow = arrow(length = unit(.2, 'cm'))) +
  scale_size_binned(breaks = c(0.1, 0.25, 0.65), range = c(0.3, 1.5),
                    labels = function(x) round(x,2)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  scale_color_scico_d(palette = "berlin", direction = -1) +
  labs(x = "", y = "", alpha = "Posterior\nprobability",
       size = "Posterior\nprobability", col = "Direction")

## plot: histogram of chosen wind directions
fig2 <- table(bag_wind_df$wind_hpp) %>%
  as.data.frame() %>%
  mutate(Prop = Freq/sum(Freq)) %>%
  ggplot() +
  geom_col(aes(Var1, Prop, fill = Var1)) +
  scale_fill_scico_d(palette = "berlin", direction = -1) +
  labs(x = "", y = "Percentage") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_y_continuous(labels = scales::percent)

gg <- gridExtra::grid.arrange(fig2, fig1, nrow = 1)
# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/sim2a_wind_hpp", ext),
#          width = 11, height = 4)
# }

################################
# Prediction results in figure #
################################
res <- data.frame(vals = c(bag$rmspe, mgp$rmspe,
                           inlares$rmspe, inlares$rmspe_ns,
                           bag$mape, mgp$mape,
                           inlares$mape, inlares$mape_ns,
                           bag$coverage, mgp$coverage,
                           inlares$coverage, inlares$coverage_ns),
                   model = rep(rep(c("G-BAG", "Q-MGP",
                                     "SPDE-stationary", "SPDE-nonstationary"),
                                   each = 25), 3),
                   what = rep(c("RMSPE", "MAPE", "95% CI coverage"), each = 100))

ggg <- res %>% ggplot() +
  geom_boxplot(aes(model, vals, fill = model), col = "black") +
  facet_wrap(~factor(what, levels = c("RMSPE", "MAPE", "95% CI coverage")),
             scales="free") +
  labs(x = "", y = "", fill = "") +
  scale_fill_manual(values = c("G-BAG" = "#1400FFFF", "Q-MGP" = "#C729D6FF",
                               "SPDE-nonstationary" = "#FF9C63FF",
                               "SPDE-stationary" = "#FFFF60FF")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 1, l = -5, r = 1, b = 0),
        legend.margin = margin(b = 0, r = 0, t = -10, l = 0))
# for (ext in extension) {
#   ggsave(plot = ggg, paste0(path, "plots/sim2a_pred", ext),
#          width = 9, height = 2.7)
# }

#################################
# Prediction results in numbers #
#################################
res %>% group_by(what, model) %>%
  summarize(mean = mean(vals))

#####################################
# Inferred directions being correct #
#####################################
# percentage of partitions choosing N (true direction) as the direction with
# highest posterior probability
mean(wind_hpp == "N")

# mean posterior probability of N across partitions choosing N
mean(val_hpp[wind_hpp == "N"])

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

################
# Running time #
################
mean(bag$est_time_per_iter[3,] + bag$pred_time_per_iter[3,])
mean(bag$est_time_per_iter[3,])
mean(bag$pred_time_per_iter[3,])
