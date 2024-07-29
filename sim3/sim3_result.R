rm(list = ls())

# dependencies 
library(tidyverse)
theme_set(theme_bw())
library(bags)

# path
path <- "~/GBAGs/"

# source
Rcpp::sourceCpp(paste0(path, "sim3/Stein.cpp"))

# plot extensions
extension <- c(".pdf", ".eps", ".png")

# seed
set.seed(723)

# locations
n_grid <- 25
n_time <- 10
xgrid <- seq(0, 1, length = n_grid)
ygrid <- xgrid
tgrid <- seq(0, 1, length = n_time)
coords <- expand.grid(easting = xgrid, northing = ygrid, time = tgrid) %>%
  arrange(time, easting, northing)
n <- nrow(coords)

# covariance
## alpha: decay
## eps: wind speed
## Z: direction
Z <- c(0.5, 1)
Z <- Z/sqrt(sum(Z^2))
cov_Stein <- Stein(coords = as.matrix(coords), alpha = 0.5, nu = 1/2, eps = 0.2, Z = Z)
cor_Stein <- cov2cor(cov_Stein)

idx <- which(coords[,1] == 0.5 & coords[,2] == 0 & coords[,3] == 0)
idx <- 1
data.frame(coords, cov = cor_Stein[idx,]) %>% 
  ggplot() + 
  geom_contour_filled(aes(easting, northing, z = cov)) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  # geom_raster(aes(easting, northing, fill = cov)) +
  # scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  facet_grid(~time, labeller = purrr::partial(label_both, sep = " = ")) +
  labs(x = "", y = "", fill = "Covariance")

##########
# Result #
##########
stein <- readRDS(paste0(path, "sim3/sim3_Stein.RDS"))
bag <- readRDS(paste0(path, "sim3/sim3_bag.RDS"))
mgp <- readRDS(paste0(path, "sim3/sim3_mgp.RDS"))
inlares <- readRDS(paste0(path, "sim3/sim3_inla.RDS"))

######################
# Prediction results #
######################
res <- data.frame(vals = c(stein$rmspe, bag$rmspe, mgp$rmspe,
                           inlares$rmspe, inlares$rmspe_ns,
                           stein$mape, bag$mape, mgp$mape,
                           inlares$mape, inlares$mape_ns,
                           rep(NA, 25), bag$coverage, mgp$coverage, 
                           inlares$coverage, inlares$coverage_ns,
                           rep(NA, 25), bag$meanwidth, mgp$meanwidth, 
                           inlares$meanwidth, inlares$meanwidth_ns),
                  model = rep(rep(c("Stein", "G-BAG", "Fixed DAG",
                                    "SPDE-stationary", "SPDE-nonstationary"),
                                  each = 25), 4),
                  what = rep(c("RMSPE", "MAPE", "95% CI coverage", "95% CI width"), each = 125)) %>% 
  mutate(model = fct_relevel(model, c("Stein", "G-BAG", "Fixed DAG",
                                      "SPDE-stationary", "SPDE-nonstationary")))
res %>% group_by(what, model) %>%
  summarize(mean = mean(vals), se = round(sd(vals), 3))

##############
# Directions #
##############
## directions
directions <- c("W", "SW", "S")

## assign partitions
n_easting <- 3
n_northing <- 3

## mcmc
mcmc <- list(save = 1000, burn = 2000, thin = 8)

## average posterior probability over all replicates
dir_mean <- data.frame(
  W = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "W")/mcmc$save))),
  SW = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "SW")/mcmc$save))),
  S = rowMeans(sapply(bag$z_save, function(z)
    apply(z, 1, function(x) sum(x == "S")/mcmc$save)))
)

## wind with the highest posterior probability
wind_hpp <- apply(dir_mean, 1, function(x) directions[which.max(x)])
val_hpp <- apply(dir_mean, 1, max)
ang_hpp <- rep(0, length(val_hpp))
ang_hpp[which(wind_hpp == "W")] <- 0
ang_hpp[which(wind_hpp == "SW")] <- pi/4
ang_hpp[which(wind_hpp == "S")] <- pi/2

## average wind direction
ang_hpp2 <- dir_mean$SW*pi/4 + dir_mean$S*pi/2

## partitions
easting_cut <- seq(0, 1, length = n_easting + 1)
northing_cut <- seq(0, 1, length = n_northing + 1)
easting_midval <- 0.5*easting_cut[1:n_easting] + 0.5*easting_cut[-1]
northing_midval <- 0.5*northing_cut[1:n_northing] + 0.5*northing_cut[-1]

## data w/ directions
bag_wind_df <- data.frame(partition = bag$ptts_tr,
                          partition1 = bag$ptts_tr,
                          wind_hpp = factor(wind_hpp,
                                            levels = directions),
                          val_hpp = val_hpp, 
                          ang_hpp = ang_hpp, 
                          ang_hpp2 = ang_hpp2) %>%
  tidyr::separate(partition1, c("row", "col", "time_d"),
                  sep = ",", convert = TRUE) %>%
  mutate(easting_center = easting_midval[col],
         northing_center = northing_midval[n_northing + 1 - row],
         time = tgrid[time_d])

time_label <- paste0("time = ", round(tgrid, 3))
names(time_label) <- tgrid
plt <- data.frame(rbind(bag$coords_tr, bag$coords_tt), 
           y = c(bag$y_tr, bag$y_pred)) %>% 
  ggplot() +
  geom_raster(aes(easting, northing, fill = y)) + 
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  facet_wrap(~time, ncol = 5, labeller = labeller(time = time_label)) +
  geom_spoke(data = bag_wind_df,
             aes(x = easting_center, y = northing_center, 
                 angle = ang_hpp2),
             color = "black",
             radius = 0.1, arrow = arrow(length = unit(.1, 'cm'))) +
  labs(x = "", y = "")

# for (ext in extension) {
#   ggsave(plot = plt, paste0(path, "plots/sim3_pred", ext),
#          width = 9.5, height = 4)
# }
