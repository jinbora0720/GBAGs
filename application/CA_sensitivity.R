# Sensitivity Analysis in CA #
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

# call results
bagres <- readRDS(paste0(path, "application/CA_daily_bag.RDS"))
bagres1 <- readRDS(paste0(path, "application/CA_daily_bag_v3.RDS"))             # different partitioning
bagres2 <- readRDS(paste0(path, "application/CA_daily_bag_v2.RDS"))             # different directions

#########
# Table #
#########
# parameter estimation
## beta
rowMeans(bagres$beta_save) %>% round(3)
quantile(bagres$beta_save, probs = c(0.025, 0.975)) %>% round(3)

rowMeans(bagres1$beta_save) %>% round(3)
quantile(bagres1$beta_save, probs = c(0.025, 0.975)) %>% round(3)

rowMeans(bagres2$beta_save) %>% round(3)
quantile(bagres2$beta_save, probs = c(0.025, 0.975)) %>% round(3)

## tausq
mean(bagres$tau_sq_save) %>% round(3)
quantile(bagres$tau_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(bagres1$tau_sq_save) %>% round(3)
quantile(bagres1$tau_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(bagres2$tau_sq_save) %>% round(3)
quantile(bagres2$tau_sq_save, probs = c(0.025, 0.975)) %>% round(3)

## sigsq
mean(bagres$sig_sq_save) %>% round(3)
quantile(bagres$sig_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(bagres1$sig_sq_save) %>% round(3)
quantile(bagres1$sig_sq_save, probs = c(0.025, 0.975)) %>% round(3)

mean(bagres2$sig_sq_save) %>% round(3)
quantile(bagres2$sig_sq_save, probs = c(0.025, 0.975)) %>% round(3)

## psi
rowMeans(bagres$psi_save) %>% round(3)
apply(bagres$psi_save, 1, function(x) quantile(x, probs = c(0.025, 0.975))) %>%
  round(3)

rowMeans(bagres1$psi_save) %>% round(3)
apply(bagres1$psi_save, 1, function(x) quantile(x, probs = c(0.025, 0.975))) %>%
  round(3)

rowMeans(bagres2$psi_save) %>% round(3)
apply(bagres2$psi_save, 1, function(x) quantile(x, probs = c(0.025, 0.975))) %>%
  round(3)

# prediction summaries
res <- data.frame(bagres$coords_ptt[,c("easting", "northing", "time")],
                  y = c(y_tr, y_tt, rep(NA, n_grid)),
                  y_bag = c(bagres$yhat, bagres$y_bag),
                  y_bag_low = c(bagres$yhat_low, bagres$y_bag_low),
                  y_bag_hi = c(bagres$yhat_hi, bagres$y_bag_hi),
                  y_bag1 = c(bagres1$yhat, bagres1$y_bag),
                  y_bag1_low = c(bagres1$yhat_low, bagres1$y_bag_low),
                  y_bag1_hi = c(bagres1$yhat_hi, bagres1$y_bag_hi), 
                  y_bag2 = c(bagres2$yhat, bagres2$y_bag),
                  y_bag2_low = c(bagres2$yhat_low, bagres2$y_bag_low),
                  y_bag2_hi = c(bagres2$yhat_hi, bagres2$y_bag_hi)) %>%
  mutate(bag_inCI = ifelse((y > y_bag_low & y < y_bag_hi), 1, 0),
         bag_CIwidth = y_bag_hi - y_bag_low, 
         bag1_inCI = ifelse((y > y_bag1_low & y < y_bag1_hi), 1, 0),
         bag1_CIwidth = y_bag1_hi - y_bag1_low, 
         bag2_inCI = ifelse((y > y_bag2_low & y < y_bag2_hi), 1, 0),
         bag2_CIwidth = y_bag2_hi - y_bag2_low)

## prediction metric 
ressum <- res[n_tr+1:n_tt,] %>%
  summarise(bag_rmspe = sqrt(mean((y-y_bag)^2)),
            bag_mape = mean(abs(y-y_bag)),
            bag_coverage = mean(bag_inCI),
            bag_meanwidth = mean(bag_CIwidth), 
            bag1_rmspe = sqrt(mean((y-y_bag1)^2)),
            bag1_mape = mean(abs(y-y_bag1)),
            bag1_coverage = mean(bag1_inCI),
            bag1_meanwidth = mean(bag1_CIwidth), 
            bag2_rmspe = sqrt(mean((y-y_bag2)^2)),
            bag2_mape = mean(abs(y-y_bag2)),
            bag2_coverage = mean(bag2_inCI),
            bag2_meanwidth = mean(bag2_CIwidth))

ressum %>% round(3)

#############
# Residuals #
#############
# examine residual plots y-w
ybar <- mean(y_tr)
e_bag <- c(y_tr-ybar, y_tt-ybar) - c(bagres$what, bagres$w_bag[1:n_tt])
e_bag1 <- c(y_tr-ybar, y_tt-ybar) - c(bagres1$what, bagres1$w_bag[1:n_tt])
e_bag2 <- c(y_tr-ybar, y_tt-ybar) - c(bagres2$what, bagres2$w_bag[1:n_tt])

coords_res <- rbind(coords_tr[,c("easting", "northing", "time", "time_d",
                                 "true_easting", "true_northing")],
                    coords_tt[,c("easting", "northing", "time", "time_d",
                                 "true_easting", "true_northing")])
e_res <- data.frame(rbind(coords_res, coords_res, coords_res),
                    e_bags = c(e_bag, e_bag1, e_bag2),
                    model = rep(c("G-BAG", "G-BAG1", "G-BAG2"), each = n)) %>%
  mutate(date = timelist[time_d])

ttselect <- c(21, 40, 72)
plt_eres1 <- e_res %>% filter(time %in% tgrid[ttselect])

## G-BAG
fig <- plt_eres1 %>%
  ggplot() + geom_point(aes(true_easting, true_northing, color = e_bags)) +
  scale_color_scico(palette = "roma", direction = -1,
                    na.value = "transparent") +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_grid(model~date) +
  labs(y="", x="", col = "Residual") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 2, l = -2, r = 2, b = -2),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -6))

# for (ext in extension) {
#   ggsave(plot = fig, paste0(path, "plots/CA_sensitivity_res", ext),
#          width = 7, height = 7)
# }

###########################
# Inferred wind direction #
###########################
z_postm <- apply(bagres$z_save, 1, getmode)
table(z_postm)/length(z_postm)

z_postm1 <- apply(bagres1$z_save, 1, getmode)
table(z_postm1)/length(z_postm1)

######################
# Predicted surfaces #
######################
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

plt_res <- (res %>% select(easting, northing, time, y_bag, y_bag1, y_bag2)) %>%
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
plt_res_Sep <- plt_res[(n_tr + n_tt) + 1:n_grid,] %>%
  filter(time %in% tgrid[ttselect])
plt_res_Sep_logy_longer <- plt_res_Sep %>% 
  pivot_longer(cols = c("y_bag", "y_bag1", "y_bag2"), 
               names_to = "model0", values_to = "y_bag") %>% 
  mutate(model = case_when(model0 == "y_bag" ~ "G-BAG",
                           model0 == "y_bag1" ~ "G-BAG1",
                           .default = "G-BAG2"))

## simply to make green area visible
plt_ymin <- max(0,min(plt_res_Sep_logy_longer$y_bag))
plt_ymax <- max(plt_res_Sep_logy_longer$y_bag)

plt_pred <- plt_res_Sep_logy_longer %>%
  ggplot() + 
  geom_point(aes(true_easting, true_northing, color = y_bag), size = 2) +
  scale_color_distiller(palette = "Spectral", direction = -1,
                        limits = c(plt_ymin, plt_ymax),
                        na.value = "transparent") +
  scale_alpha_binned(n.breaks = 3, range=c(.1,1), guide = "none") +
  geom_sf(data = ca_sf, fill = "NA") +
  facet_grid(model~date) +
  labs(y="", x="", color = "log(PM2.5)") +
  theme(legend.position = "right",
        plot.margin=margin(t = 2, l = -5, r = 2, b = -5),
        legend.margin=margin(b = 0, r = 0, t = 0, l = -6)) +
  scale_x_continuous(breaks = c(-123, -121, -119, -117, -115),
                     labels = expression(123 * degree, 121 * degree,
                                         119 * degree, 117 * degree,
                                         115 * degree * W))

# for (ext in extension) {
#   ggsave(plot = plt_pred, paste0(path, "plots/CA_sensitivity_pred_final", ext),
#          width = 9.3, height = 6.5)
# }

#################
# Miscellaneous #
#################
# directions
directions <- c("W", "NW", "N", "NE")

# mcmc
mcmc <- list(save = 1000, burn = 15000, thin = 15)

# original partition 
n_easting <- 16 # (~ 55km, ~0.5 degree resolution)
n_northing <- 20
easting_cut <- seq(ca_range$xmin, ca_range$xmax, length = n_easting + 1)
northing_cut <- seq(ca_range$ymin, ca_range$ymax, length = n_northing + 1)
easting_midval <- 0.5*easting_cut[1:n_easting] + 0.5*easting_cut[-1]
northing_midval <- 0.5*northing_cut[1:n_northing] + 0.5*northing_cut[-1]

## wind with the highest posterior probability
dir_mean <- data.frame(
  W = apply(bagres$z_save, 1, function(x) sum(x == "W")/mcmc$save),
  NW = apply(bagres$z_save, 1, function(x) sum(x == "NW")/mcmc$save),
  N = apply(bagres$z_save, 1, function(x) sum(x == "N")/mcmc$save),
  NE = apply(bagres$z_save, 1, function(x) sum(x == "NE")/mcmc$save)
)

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

# G-BAG1 partition 
n_easting1 <- 9 
n_northing1 <- 11
easting_cut1 <- seq(ca_range$xmin, ca_range$xmax, length = n_easting1 + 1)
northing_cut1 <- seq(ca_range$ymin, ca_range$ymax, length = n_northing1 + 1)
easting_midval1 <- 0.5*easting_cut1[1:n_easting1] + 0.5*easting_cut1[-1]
northing_midval1 <- 0.5*northing_cut1[1:n_northing1] + 0.5*northing_cut1[-1]

## wind with the highest posterior probability
dir_mean1 <- data.frame(
  W = apply(bagres1$z_save, 1, function(x) sum(x == "W")/mcmc$save),
  NW = apply(bagres1$z_save, 1, function(x) sum(x == "NW")/mcmc$save),
  N = apply(bagres1$z_save, 1, function(x) sum(x == "N")/mcmc$save),
  NE = apply(bagres1$z_save, 1, function(x) sum(x == "NE")/mcmc$save)
)

### wind name
wind_hpp1 <- apply(dir_mean1, 1, function(x) directions[which.max(x)])
### posterior prob
val_hpp1 <- apply(dir_mean1, 1, max)
### angle
ang_hpp1 <- rep(0, length(val_hpp1))
ang_hpp1[which(wind_hpp1 == "NW")] <- 7*pi/4
ang_hpp1[which(wind_hpp1 == "N")] <- 3*pi/2
ang_hpp1[which(wind_hpp1 == "NE")] <- 5*pi/4

## winds in partitioned regions
ptts_tr1 <- sort(unique(bagres1$coords_ptt$partition[1:n_tr]))
bag_wind_df1 <- data.frame(partition = ptts_tr1,
                           partition1 = ptts_tr1,
                           val_hpp = val_hpp1,
                           ang_hpp = ang_hpp1) %>%
  tidyr::separate(partition1,
                  c("row", "col", "time_d"), sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval1[col],
         northing_center = northing_midval1[n_northing1 + 1 - row],
         time = timelist[time_d])

# G-BAG2 wind directions 
directions2 <- c("W", "SW", "S", "SE")  

## wind with the highest posterior probability
dir_mean2 <- data.frame(
  W = apply(bagres2$z_save, 1, function(x) sum(x == "W")/mcmc$save),
  SW = apply(bagres2$z_save, 1, function(x) sum(x == "SW")/mcmc$save),
  S = apply(bagres2$z_save, 1, function(x) sum(x == "S")/mcmc$save),
  SE = apply(bagres2$z_save, 1, function(x) sum(x == "SE")/mcmc$save)
)

### wind name
wind_hpp2 <- apply(dir_mean2, 1, function(x) directions2[which.max(x)])
### posterior prob
val_hpp2 <- apply(dir_mean2, 1, max)
### angle
ang_hpp2 <- rep(0, length(val_hpp2))
ang_hpp2[which(wind_hpp2 == "SW")] <- pi/4
ang_hpp2[which(wind_hpp2 == "S")] <- pi/2
ang_hpp2[which(wind_hpp2 == "SE")] <- 3*pi/4

## winds in partitioned regions
bag_wind_df2 <- data.frame(partition = ptts_tr,
                           partition1 = ptts_tr,
                           val_hpp = val_hpp2,
                           ang_hpp = ang_hpp2) %>%
  tidyr::separate(partition1,
                  c("row", "col", "time_d"), sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row],
         time = timelist[time_d])

# 1. August Complex
ttselect <- c(17, 18, 20, 22)
plt_res_AC <- plt_res[(n_tr + n_tt) + 1:n_grid,] %>%
  filter(time %in% tgrid[ttselect])
plt_res_ACwind <- bag_wind_df %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)
plt_res_ACwind1 <- bag_wind_df1 %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)

plt_res_AC_long <- plt_res_AC %>%
  select(true_easting, true_northing, date, y_bag, y_bag1) %>%
  pivot_longer(cols = c("y_bag", "y_bag1"),
               names_to = "model0", values_to = "y_bag") %>%
  mutate(z = exp(y_bag)) %>%
  select(-y_bag) %>%
  mutate(model = case_when(model0 == "y_bag" ~ "G-BAG",
                           .default = "G-BAG1"))

ggg <- plt_res_AC_long %>% ggplot() +
  geom_contour_filled(aes(true_easting, true_northing, z = z),
                      breaks = c(0, 12, 35, 65, Inf)) +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442", "#D55E00"),
                    labels = c("(0,12]","(12,35]","(35,65]",">65")) +
  geom_sf(data = ca_sf, fill = "NA") +
  geom_spoke(data = plt_res_ACwind %>% mutate(model = "G-BAG"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000, col = "white", linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_ACwind %>% mutate(model = "G-BAG"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_ACwind1 %>% mutate(model = "G-BAG1"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000, col = "white", linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_ACwind1 %>% mutate(model = "G-BAG1"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  facet_grid(model ~ date) +
  labs(y = "", x = "", fill = "PM2.5 (um/m3)") +
  scale_x_continuous(breaks = c(-123, -121, -119)) +
  coord_sf(xlim = c(-135000, 450000), ylim = c(4100000, 4700000), expand = FALSE) +
  theme(legend.position = "bottom",
        legend.margin = margin(b = 0, r = 0, t = -15, l = 0))

# 2. Discretized predicted surfaces with winds in Sep
ttselect <- c(44:45, 47:49)
plt_res_Sepwind <- bag_wind_df %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)
plt_res_Sepwind1 <- bag_wind_df1 %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)
plt_res_Sepwind2 <- bag_wind_df2 %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)

plt_pred_dis <- plt_res_Sep_logy_longer %>%
  ggplot() +
  geom_contour_filled(aes(true_easting, true_northing, z = exp(y_bag)),
                      breaks = c(Inf, 65, 35, 12, 0)) +
  scale_fill_manual(values = c("#D55E00","#F0E442","#009E73", "#56B4E9"),
                    labels = c(">65","(35,65]","(12,35]","(0,12]")) +
  geom_sf(data = ca_sf, fill = "NA") +
  geom_spoke(data = plt_res_Sepwind %>% mutate(model = "G-BAG"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             color = "white",
             radius = 40000, linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_Sepwind %>% mutate(model = "G-BAG"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_Sepwind1 %>% mutate(model = "G-BAG1"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             color = "white",
             radius = 40000, linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_Sepwind1 %>% mutate(model = "G-BAG1"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_Sepwind2 %>% mutate(model = "G-BAG2"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             color = "white",
             radius = 40000, linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_spoke(data = plt_res_Sepwind2 %>% mutate(model = "G-BAG2"),
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             radius = 40000,
             arrow = arrow(length = unit(.1, 'cm'))) +
  facet_grid(model~date) +
  labs(y = "", x = "", fill = "PM2.5\n(um/m3)", alpha = "Posterior\nprobability") +
  theme(legend.position = "right",
        plot.margin = margin(t = 2,l = -2, r = 2, b = -5),
        legend.margin = margin(b = 0,r = 0,t = 0,l = -6)) +
  scale_x_continuous(breaks = c(-123, -121, -119, -117, -115),
                     labels = expression(123 * degree, 121 * degree,
                                         119 * degree, 117 * degree,
                                         115 * degree * W))

# 3. Replicate Figure 1 with G-BAG2
ttselect <- c(34:36)
plt_res_Fig1wind2 <- bag_wind_df2 %>%
  filter(time_d %in% ttselect) %>%
  rename(date = time)

CA_smoke <- readRDS("~/Air_pollution/Data/CA_hms_smoke_20200810_20200916.rds")
tmp2 <- CA_smoke %>% 
  filter(date %in% c(as.Date("2020-09-03"), as.Date("2020-09-04"), as.Date("2020-09-05")))

ca_sf %>%
  ggplot() +
  geom_sf(data = tmp2, aes(fill = factor(Density)), alpha = 0.4) +
  scale_fill_manual(values = c("5" = "gray70", "16" = "gray40", "27" = "gray10"), 
                    name = "Smoke density") + 
  geom_spoke(data = plt_res_Fig1wind2,
             aes(x = easting_center, y = northing_center, angle = ang_hpp),
             color = "white",
             radius = 40000, # linewidth = 1.5,
             arrow = arrow(length = unit(.1, 'cm'))) +
  geom_sf(fill = "NA") +
  facet_wrap(~date, nrow=1) +
  labs(y = "Latitude", x = "Longitude") +
  theme(legend.position = "bottom", 
        legend.margin = margin(b = 0, r = 0, t = -6, l = 0), 
        legend.text = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(-124, -112, by = 4))
