rm(list = ls())

# dependencies 
library(tidyverse)
theme_set(theme_bw())
library(meshed) # meshed_0.2.tar

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

idx <- which(coords[,1] == 1 & coords[,2] == 0 & coords[,3] == 0)
idx <- 1
data.frame(coords, cov = cor_Stein[idx,]) %>% 
  ggplot() + 
  geom_contour_filled(aes(easting, northing, z = cov)) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  # geom_raster(aes(easting, northing, fill = cov)) +
  # scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  facet_grid(~time, labeller = purrr::partial(label_both, sep = " = ")) +
  labs(x = "", y = "", fill = "Covariance")

# fit MGP
## % of train data 
prob <- 0.8

## assign partitions
n_easting <- 3
n_northing <- 3

## seed
set.seed(723)
seedsave <- c(sample(10000, 24), 723)

## save results
rmspe = mape = coverage = meanwidth = tausq_hat = sigsq_hat <- rep(0, 25)
psi_hat = tot_time <- matrix(0, nrow = 3, ncol = 25)

## mcmc
mcmc <- list(save = 1000, burn = 2000, thin = 8)

## prior
tmaxdist <- max(dist(coords[,3], method = "manhattan"))
la <- (1/0.95-1)/tmaxdist                                                       # correlation = 0.95 at the maximum distance
ua <- (1/0.05-1)/tmaxdist                                                       # correlation = 0.05 at the maximal distance
spmaxdist <- max(dist(coords[,1:2]))
lc <- -log(0.95)/spmaxdist                                                      # correlation = 0.95 at the maximum distance
uc <- -log(0.05)/spmaxdist                                                      # correlation = 0.05 at the maximum distance

## run
save_data <- FALSE
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
  
  # mgp
  y0 <- data$y
  y0[-tr_idx] <- NA
  x <- matrix(1, nrow = n, ncol = 1)
  
  time_mgp <- system.time({
    meshout <- spmeshed(family = c("gaussian"),
                        y = y0, x = x,
                        coords = data[, c("easting", "northing", "time")],
                        k = 1,
                        axis_partition = c(n_easting, n_northing, n_time),
                        prior = list(phi = c(min(lc,la), max(uc,ua))),
                        settings = list(adapting = TRUE, forced_grid = NULL,
                                        cache = NULL, ps = TRUE, saving = TRUE),
                        starting = list(beta = NULL, tausq = NULL, theta = NULL,
                                        lambda = matrix(1,1,1), w = NULL,
                                        nu = NULL, mcmcsd = .05,
                                        mcmc_startfrom = 0),
                        debug = list(sample_beta = TRUE, sample_tausq = TRUE,
                                     sample_theta = TRUE, sample_w = TRUE,
                                     sample_lambda = FALSE,
                                     verbose = FALSE, debug = FALSE),
                        n_samples = mcmc$save,
                        n_burn = mcmc$burn,
                        n_thin = mcmc$thin,
                        verbose = 0,
                        n_threads = 10)
  })
  
  # total time spent in minutes
  tot_time[,s] <- (time_mgp/60)[1:3]
  
  # parameter estimation summary
  tausq_mcmc <- meshout$tausq_mcmc[,mcmc$thin*(1:mcmc$save)]
  psi_mcmc <- meshout$theta_mcmc[-4,,mcmc$thin*(1:mcmc$save)]                   # temporal decay, spatial decay, interaction
  sigsq_mcmc <- meshout$theta_mcmc[4,,mcmc$thin*(1:mcmc$save)]
  
  tausq_hat[s] <- mean(tausq_mcmc)
  psi_hat[,s] <- rowMeans(psi_mcmc)
  sigsq_hat[s] <- mean(sigsq_mcmc)
  
  # prediction summary
  y_mgp <- meshout$yhat_mcmc %>% summary_list_mean() 
  y_mgp_low <- meshout$yhat_mcmc %>% summary_list_q(q = 0.025) 
  y_mgp_hi <- meshout$yhat_mcmc %>% summary_list_q(q = 0.975) 
  mgp_df_out <- data.frame(rbind(coords_tr[,c("easting", "northing", "time")],
                                 coords_tt[,c("easting", "northing", "time")]),
                           y = c(y_tr, y_tt)) %>%
    left_join(cbind(meshout$coordsdata, y_mgp, y_mgp_low, y_mgp_hi))
  
  mgp_res <- mgp_df_out[-c(1:n_tr),] %>%
    mutate(cover = ifelse((y > y_mgp_low & y < y_mgp_hi), 1, 0),
           range = y_mgp_hi - y_mgp_low) %>%
    summarise(rmspe = sqrt(mean((y-y_mgp)^2)),
              mape = mean(abs(y-y_mgp)),
              coverage = mean(cover),
              meanwidth = mean(range))
  
  coverage[s] <- mgp_res$coverage
  meanwidth[s] <- mgp_res$meanwidth
  rmspe[s] <- mgp_res$rmspe
  mape[s] <- mgp_res$mape
  
  setTxtProgressBar(pb, s/25)
}
close(pb)

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             tot_time = tot_time,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             meshout = meshout),                                                # seed = 3
        paste0(path, "sim3/sim3_mgp_all.RDS"))

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth),
        paste0(path, "sim3/sim3_mgp.RDS"))
