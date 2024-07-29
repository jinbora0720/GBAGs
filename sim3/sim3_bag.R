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

# fit G-BAG
## % of train data 
prob <- 0.8

## assign partitions
n_easting <- 3
n_northing <- 3

## seed
set.seed(723)
seedsave <- c(sample(10000, 24), 723)

## save results
rmspe = mape = coverage = meanwidth <- rep(0, 25)
tausq_hat = tausq_ess =
  sigsq_hat = sigsq_ess = sigsqc_ess <- rep(0, 25)
psi_hat = psi_ess <- matrix(0, nrow = 3, ncol = 25)
est_time_per_iter = pred_time_per_iter <- matrix(0, nrow = 3, ncol = 25)
rownames(est_time_per_iter) = rownames(pred_time_per_iter) <-
  c("user", "system", "elapsed")
z_pvalue = z_save = y_pred_ess <- list()

## mcmc
mcmc <- list(save = 1000, burn = 2000, thin = 8)

## directions
directions <- c("W", "S", "SW")

## prior
tmaxdist <- max(dist(coords[,3], method = "manhattan"))
la <- (1/0.95-1)/tmaxdist                                                       # correlation = 0.95 at the maximum distance
ua <- (1/0.05-1)/tmaxdist                                                       # correlation = 0.05 at the maximum distance
spmaxdist <- max(dist(coords[,1:2]))
lc <- -log(0.95)/spmaxdist                                                      # correlation = 0.95 at the maximum distance
uc <- -log(0.05)/spmaxdist                                                      # correlation = 0.05 at the maximum distance

## z Chisq test
prop <- 0.35
n_samp <- ceiling(0.25*mcmc$save)
n1 = n2 <- n_samp

## run
save_data <- FALSE
pb = txtProgressBar(style=3,width=50)
for (s in 1:25) {
  if (s == 25) {
    save_data <- TRUE
  }
  
  # set seed
  seed <- seedsave[s]
  set.seed(seed)
  
  # generate data 
  y <- t(chol(cov_Stein + diag(0.01, n)))%*%rnorm(n)
  data <- data.frame(coords, y)
  
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
  out <- bag(y = y_tr, X = NULL,
             coords = as.matrix(coords_tr),
             X_pred = NULL,
             coords_pred = as.matrix(coords_tt),
             n_partition = c(n_easting, n_northing, n_time),
             breaks_partition = list(breaks_easting = NULL,
                                     breaks_northing = NULL,
                                     breaks_time = NULL),
             directions = directions,
             init = list(tau_sq = NULL,
                         sig_sq = NULL,
                         w = NULL,
                         z = NULL,
                         psi = NULL,
                         Sn = NULL),
             hyper = list(at = NULL, bt = NULL,
                          as = NULL, bs = NULL,
                          la = la, ua = ua,
                          lc = lc, uc = uc,
                          mu0 = NULL, invV0 = NULL),
             mcmc = mcmc,
             n_threads = 10,
             seed = seed,
             verbose = FALSE,
             save_data = save_data,
             save_est = FALSE,
             debug = list(psi_fixed = FALSE, z_fixed = FALSE))
  
  # time spent per iteration
  est_time_per_iter[,s] <- as.numeric(out$est_time/out$est_iter)[1:3]
  pred_time_per_iter[,s] <- as.numeric(out$pred_time/out$pred_iter)[1:3]
  
  # parameter estimation summary
  tausq_hat[s] <- mean(out$tau_sq_save)
  psi_hat[,s] <- rowMeans(out$psi_save)
  sigsq_hat[s] <- mean(out$sig_sq_save)
  z_save[[s]] <- out$z_save
  
  # effective sample size
  tausq_ess[s] <- coda::effectiveSize(out$tau_sq_save)
  sigsq_ess[s] <- coda::effectiveSize(out$sig_sq_save)
  psi_ess[,s] <- apply(out$psi_save, 1, coda::effectiveSize)
  sigsqc_ess[s] <- coda::effectiveSize(out$sig_sq_save*out$psi_save[2,])
  y_pred_ess[[s]] <- apply(out$y_pred_save, 1, coda::effectiveSize)
  pvalue <- rep(0, nrow(out$z_save))
  for (ii in 1:nrow(out$z_save)) {
    samp1 <- out$z_save[ii,1:n_samp]
    samp2 <- out$z_save[ii,(mcmc$save - n_samp) + 1:n_samp]
    
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
    pvalue[ii] <- pchisq(X_sq, df = (nR-1)*(2-1), lower.tail = FALSE)
  }
  z_pvalue[[s]] <- pvalue
  
  # prediction summary
  y_bag <- rowMeans(out$y_pred_save)
  y_bag_qt <- apply(out$y_pred_save, 1, function(x)
    quantile(x, probs = c(0.025, 0.975)))
  y_bag_low <- y_bag_qt[1,]
  y_bag_hi <- y_bag_qt[2,]
  
  coverage[s] <- mean(ifelse(y_tt > y_bag_low & y_tt < y_bag_hi, 1, 0))
  meanwidth[s] <- mean(y_bag_hi - y_bag_low)
  rmspe[s] <- sqrt(mean((y_tt-y_bag)^2))
  mape[s] <- mean(abs(y_tt-y_bag))

  setTxtProgressBar(pb, s/25)
}
close(pb)

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             est_time_per_iter = est_time_per_iter,
             pred_time_per_iter = pred_time_per_iter,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             tausq_ess = tausq_ess,
             psi_ess = psi_ess,
             sigsq_ess = sigsq_ess,
             sigsqc_ess = sigsqc_ess,
             z_save = z_save,
             z_pvalue = z_pvalue,
             y_pred_ess = y_pred_ess,
             out = out),
        paste0(path, "sim3/sim3_bag_all.RDS"))

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth, 
             z_save = z_save,
             y_tr = y_tr,                                                       # last seed
             y_tt = y_tt, 
             coords_tr = coords_tr, 
             coords_tt = coords_tt,
             ptts_tr = sort(unique(out$coords_ptt$partition[1:n_tr])), 
             y_pred = y_bag),
        paste0(path, "sim3/sim3_bag.RDS"))

