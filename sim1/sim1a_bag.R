# Bag of DAGs Simulations #
#########
# Model #
#########
# y(l) = X(l)^T %*% beta + w(l) + eps(l), i = 1, ..., n
# l = (s, t)
# X(l)^T: (1 x p),
# eps(l) ~ N(0, tau_sq)
# w|z=h ~ N(0, sig_sq*Ctilde_h(phi)) for h = 1, ..., K where
# Ctilde_h(phi)^{-1} = (I-H_h(phi))^T %*% R_h(phi)^{-1} %*% (I-H_h(phi))
# H_h(phi) is a sparse matrix whose non-zero elements' locations are determined by the hth DAG.
# R_h(phi) is a diagonal matrix.
# H_h(phi) and R_h(phi) are conditional mean coefficient and conditional variance of normal distributions.
# Assume a nonseparable, stationary, and isotropic spatio-temporal correlation function
# Cor(w(l), w(l + (h,u))) = exp(-c||h||/(a|u|+1)^{kappa/2})/(a|u|+1)
# c: spatial decay, c > 0
# a: temporal decay, a > 0
# kappa: interaction between space and time

rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(bags)

# path
path <- "~/GBAGs/"

########
# Data #
########
# specify number of grid on each axis
ngrid <- 40
n_time <- 8
xgrid <- seq(0, 1, length = ngrid)
ygrid <- xgrid
tgrid <- seq(0, 1, length = n_time)
xytgrid <- expand.grid(easting = xgrid, northing = ygrid, time = tgrid) %>%
  arrange(time, easting, northing)
n <- nrow(xytgrid)
coords <- xytgrid

# assign partitions (irregular partitions)
n_easting <- 6
n_northing <- 6

# directions
directions <- c("W", "NW", "N", "NE")

# data
nd <- floor(log10(max(n_northing, n_easting, n_time))) + 1
format <- paste0("%0", nd, ".0f,","%0", nd, ".0f,", "%0", nd, ".0f")

coords_ptt <- coords
northing_sint <- seq(0, 1, length = (n_northing-2)^2+1)
northing_Lint <- seq(0, 1, length = (n_northing-2)+1)
breaks_northing <- c(northing_sint[1:2], northing_Lint[2:(n_northing-2)],
                     northing_sint[(n_northing-2)^2+0:1])
coords_ptt$row <- (n_northing+1) -
  as.numeric(cut(coords_ptt$northing,
                 breaks = breaks_northing,
                 labels = 1:n_northing, include.lowest = T)) # 1 from the top
easting_sint <- seq(0, 1, length = (n_easting-2)^2+1)
easting_Lint <- seq(0, 1, length = (n_easting-2)+1)
breaks_easting <- c(easting_sint[1:2], easting_Lint[2:(n_easting-2)],
                    easting_sint[(n_easting-2)^2+0:1])
coords_ptt$col <- as.numeric(cut(coords_ptt$easting,
                                 breaks = breaks_easting,
                                 labels = 1:n_easting, include.lowest = T))
coords_ptt$time_d <- as.numeric(as.character(coords_ptt$time*(n_time-1) + 1))
coords_ptt$partition <- sprintf(format,
                                coords_ptt$row, coords_ptt$col,
                                coords_ptt$time_d)

ptts <- sort(unique(coords_ptt$partition))
z_true <- rep("W", length(ptts))
ptt_NW <- grep("\\d+,4,2|\\d+,5,2|\\d+,6,2|\\d+,\\d+,3|1,\\d+,4|2,\\d+,4|3,\\d+,4|4,2,8|4,3,8|4,4,8|4,5,8|4,6,8|5,6,8|6,6,8", ptts)
z_true[ptt_NW] <- "NW"
ptt_N <- grep("4,\\d+,4|5,\\d+,4|6,\\d+,4|\\d+,\\d+,5|1,\\d+,6|2,\\d+,6|3,\\d+,6|3,1,8|3,2,8|3,3,8|3,4,8|3,5,8|4,1,8", ptts)
z_true[ptt_N] <- "N"
ptt_NE <- grep("4,\\d+,6|5,\\d+,6|6,\\d+,6|\\d+,\\d+,7|5,1,8|5,2,8|5,3,8|5,4,8|5,5,8|6,1,8|6,2,8|6,3,8|6,4,8|6,5,8", ptts)
z_true[ptt_NE] <- "NE"

# mcmc
mcmc <- list(save = 1000, burn = 10000, thin = 7)

# z Chisq test
prop <- 0.35
n_samp <- ceiling(0.25*mcmc$save)
n1 = n2 <- n_samp

# true parameters
tau_sq <- 0.01
sig_sq <- 2
a <- 5
c <- 0.5
kappa <- 0.9
beta <- 2

# prior
tmaxdist <- max(dist(coords[,3], method = "manhattan"))
la <- (1/0.2-1)/tmaxdist                                                        # correlation = 0.2 at the maximum distance
ua <- (1/0.2-1)/(0.5*tmaxdist)                                                  # correlation = 0.2 at the half maximal distance
spmaxdist <- max(dist(coords[,1:2]))
lc <- -log(0.8)/(spmaxdist)                                                     # correlation = 0.8 at the maximum distance
uc <- -log(0.8)/(0.2*spmaxdist)                                                 # correlation = 0.8 at one fifth of the maximum distance

# generate seed
set.seed(123)
seedsave <- sample(10000, 100)
seedsave <- c(seedsave[1:24], 123)

# save results
rmspe = mape = coverage = meanwidth <- rep(0, 25)
beta_hat = beta_ci = beta_ess =
  tausq_hat = tausq_ci = tausq_ess =
  sigsq_hat = sigsq_ci = sigsq_ess = sigsqc_ess <- rep(0, 25)
psi_hat = psi_ci = psi_ess <- matrix(0, nrow = 3, ncol = 25)
est_time_per_iter = pred_time_per_iter <- matrix(0, nrow = 3, ncol = 25)
rownames(est_time_per_iter) = rownames(pred_time_per_iter) <-
  c("user", "system", "elapsed")
z_pvalue = z_save = y_pred_ess <- list()

pb = txtProgressBar(style=3,width=50)
for (s in 1:25) {
  seed <- seedsave[s]
  set.seed(seed)

  data <- rbag(coords = as.matrix(coords),
               n_partition = c(n_easting, n_northing, n_time),
               breaks_partition = list(breaks_easting = breaks_easting,
                                       breaks_northing = breaks_northing,
                                       breaks_time = NULL),
               z = z_true,
               params = list(tau_sq = tau_sq,
                             sig_sq = sig_sq,
                             a = a, c = c, kappa = kappa,
                             beta = beta),
               seed = seed)

  #################
  # Training data #
  #################
  # train data
  prob <- 0.8
  grid0 <- expand.grid(easting = xgrid, northing = ygrid)

  idx_tr_tmp <- sort(sample(nrow(grid0), nrow(grid0)*prob))
  xygrid_tr <- list()
  for (tt in 1:n_time) {
    xygrid_tr[[tt]] <- data.frame(grid0[idx_tr_tmp,], time = tgrid[[tt]])
  }
  xygrid_tr <- do.call(rbind, xygrid_tr) %>% arrange(easting, northing, time)
  data_tr <- left_join(xygrid_tr,
                       data %>% mutate(idx = 1:n),
                       by = c("easting", "northing", "time"))
  tr_idx <- data_tr$idx
  n_tr <- length(tr_idx)
  n_tt <- n - n_tr
  w_tr <- data_tr$w
  w_tt <- data$w[-tr_idx]
  y_tr <- data_tr$y
  y_tt <- data$y[-tr_idx]
  X_tr <- data_tr$X
  X_tt <- data$X[-tr_idx]

  coords_tr <- data_tr[, c("easting", "northing", "time")]
  coords_tt <- data[-tr_idx, c("easting", "northing", "time")]

  ###########
  # Methods #
  ###########
  # bag
  save_data <- FALSE
  if (s == 25) {
    save_data <- TRUE
  }
  out <- bag(y = y_tr, X = as.matrix(X_tr),
             coords = as.matrix(coords_tr[,c("easting", "northing", "time")]),
             X_pred = as.matrix(X_tt),
             coords_pred = as.matrix(coords_tt[,c("easting", "northing", "time")]),
             n_partition = c(n_easting, n_northing, n_time),
             breaks_partition = list(breaks_easting = breaks_easting,
                                     breaks_northing = breaks_northing,
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
  beta_hat[s] <- mean(out$beta_save[1,])
  tausq_hat[s] <- mean(out$tau_sq_save)
  psi_hat[,s] <- rowMeans(out$psi_save)
  sigsq_hat[s] <- mean(out$sig_sq_save)
  z_save[[s]] <- out$z_save

  # 95% CI coverage
  beta_ci[s] <- (beta > quantile(out$beta_save[1,], 0.025) &
                   beta < quantile(out$beta_save[1,], 0.975))
  tausq_ci[s] <- (tau_sq > quantile(out$tau_sq_save, 0.025) &
                    tau_sq < quantile(out$tau_sq_save, 0.975))
  psi_ci[1,s] <- (a > quantile(out$psi_save[1,], 0.025) &
                    a < quantile(out$psi_save[1,], 0.975))
  psi_ci[2,s] <- (c > quantile(out$psi_save[2,], 0.025) &
                    c < quantile(out$psi_save[2,], 0.975))
  psi_ci[3,s] <- (kappa > quantile(out$psi_save[3,], 0.025) &
                    kappa < quantile(out$psi_save[3,], 0.975))
  sigsq_ci[s] <- (sig_sq > quantile(out$sig_sq_save, 0.025) &
                    sig_sq < quantile(out$sig_sq_save, 0.975))

  # effective sample size
  beta_ess[s] <- coda::effectiveSize(out$beta_save[1,])
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
             beta_hat = beta_hat,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             beta_ci = beta_ci,
             tausq_ci = tausq_ci,
             psi_ci = psi_ci,
             sigsq_ci = sigsq_ci,
             beta_ess = beta_ess,
             tausq_ess = tausq_ess,
             psi_ess = psi_ess,
             sigsq_ess = sigsq_ess,
             sigsqc_ess = sigsqc_ess,
             z_save = z_save,
             z_pvalue = z_pvalue,
             y_pred_ess = y_pred_ess,
             out = out),                                                        # seed = 123
        paste0(path, "sim1/sim1a_bag.RDS"))


