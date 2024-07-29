rm(list = ls())

# dependencies 
library(tidyverse)
theme_set(theme_bw())
library(bags)

# path
path <- "~/GBAGs/"

# source
Rcpp::sourceCpp(paste0(path, "sim3/Stein.cpp"))

# locations
n_grid <- 25
n_time <- 10
xgrid <- seq(0, 1, length = n_grid)
ygrid <- xgrid
xygrid <- expand.grid(easting = xgrid, northing = ygrid)
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

# kriging 
## % of train data 
prob <- 0.8

## seed
set.seed(723)
seedsave <- c(sample(10000, 24), 723)

## save results
rmspe = mape = coverage = meanwidth <- rep(0, 25)

## run
pb = txtProgressBar(style=3,width=50)
for (s in 1:25) {
  # set seed
  seed <- seedsave[s]
  set.seed(seed)
  
  # generate data 
  y <- t(chol(cov_Stein + diag(0.01, n)))%*%rnorm(n)
  data <- data.frame(coords, y)
  
  # time_label <- paste0("time = ", round(tgrid, 3))
  # names(time_label) <- tgrid
  # data %>%
  #   ggplot() +
  #   geom_raster(aes(easting, northing, fill = y)) +
  #   scale_fill_distiller(palette = "RdBu", direction = -1) +
  #   facet_wrap(~time, ncol = 5, labeller = labeller(time = time_label)) +
  #   labs(x = "", y = "")
  
  # train data vs test data
  idx_tr_tmp <- sort(sample(nrow(xygrid), nrow(xygrid)*prob))
  xygrid_tr <- list()
  for (tt in 1:n_time) {
    xygrid_tr[[tt]] <- data.frame(xygrid[idx_tr_tmp,], time = tgrid[[tt]])
  }
  xygrid_tr <- do.call(rbind, xygrid_tr) %>% arrange(easting, northing, time)
  data_tr <- left_join(xygrid_tr,
                       data %>% mutate(idx = 1:n),
                       by = c("easting", "northing", "time"))
  tr_idx <- data_tr$idx
  n_tr <- length(tr_idx)
  n_tt <- n - n_tr
  y_tr <- data_tr$y
  y_tt <- data$y[-tr_idx]
  
  coords_tr <- data_tr[, c("easting", "northing", "time")]
  coords_tt <- data[-tr_idx, c("easting", "northing", "time")]
  
  # fit 
  SSinv <- cov_Stein[-tr_idx, tr_idx]%*%solve(cov_Stein[tr_idx, tr_idx] + diag(0.01, n_tr))
  y_krig <- SSinv%*%y_tr
  rmspe[s] <- sqrt(mean((y_tt-y_krig)^2))
  mape[s] <- mean(abs(y_tt-y_krig))

  setTxtProgressBar(pb, s/25)
}
close(pb)

saveRDS(list(rmspe = rmspe,
             mape = mape),
        paste0(path, "sim3/sim3_Stein.RDS"))
