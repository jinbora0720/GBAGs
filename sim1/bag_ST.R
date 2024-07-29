bag_ST <- function(y, X = NULL, coords,
                X_pred = NULL, coords_pred = NULL,
                n_partition,
                breaks_partition = list(breaks_easting = NULL,
                                        breaks_northing = NULL,
                                        breaks_time = NULL),
                directions,
                init = list(tau_sq = NULL,
                            sig_sq = NULL,
                            w = NULL,
                            z = NULL,
                            psi = NULL,
                            Sn = NULL),
                hyper = list(at = NULL, bt = NULL,
                             as = NULL, bs = NULL,
                             la = NULL, ua = NULL,
                             lc = NULL, uc = NULL,
                             mu0 = NULL, invV0 = NULL),
                mcmc = list(save = 1000, burn = 0, thin = 1),
                n_threads = 10,
                seed = 123,
                verbose = TRUE,
                save_data = TRUE,
                save_est = FALSE,
                debug = list(psi_fixed = FALSE, z_fixed = FALSE)) {

  ##################
  # to get started #
  ##################
  # set seed
  set.seed(seed)

  # make y have zero mean
  ybar <- mean(y)
  y <- y - ybar

  # dimension of data
  if (is.numeric(y)) { n <- length(y) } else { n <- nrow(y) }
  if (!is.null(X)) p <- ncol(X)
  if (is.null(directions)) {directions <- c("SW", "W", "NW", "N")}
  K <- length(directions)
  d <- ncol(coords)
  if (d == 2) {
    coords <- cbind(coords, 0)
  }

  # number of partitions
  if (is.null(breaks_partition$breaks_easting)) {
    n_easting <- n_partition[1]
  } else {
    breaks_easting <- breaks_partition$breaks_easting
    n_easting <- length(breaks_easting) - 1
  }
  if (is.null(breaks_partition$breaks_northing)) {
    n_northing <- n_partition[2]
  } else {
    breaks_northing <- breaks_partition$breaks_northing
    n_northing <- length(breaks_northing) - 1
  }
  if (d == 2) {
    n_time <- 1
  } else {
    if (is.null(breaks_partition$breaks_time)) {
      n_time <- n_partition[3]
    } else {
      breaks_time <- breaks_partition$breaks_time
      n_time <- length(breaks_time) - 1
    }
  }

  # mcmc parameters
  thin <- mcmc$thin
  burn <- mcmc$burn
  save <- mcmc$save
  S <- burn + thin*save

  # save results
  out <- list()

  #################
  # preprocessing #
  #################
  if (verbose) {
    cat('preprocessing...\n')
  }
  nd <- floor(log10(max(n_easting, n_northing, n_time))) + 1
  coords_all_ptt <- data.frame(rbind(coords, coords_pred))
  rownames(coords_all_ptt) <- NULL
  colnames(coords_all_ptt) <- c("easting", "northing", "time")
  if (is.null(breaks_partition$breaks_northing)) {                              # row first col later
    coords_all_ptt$row <- (n_northing+1) -
      as.numeric(cut_interval(coords_all_ptt$northing,
                              n = n_northing,
                              labels = 1:n_northing))                           # 1 from the top

  } else {
    coords_all_ptt$row <- (n_northing+1) -
      as.numeric(cut(coords_all_ptt$northing,
                     breaks = breaks_northing,
                     labels = 1:n_northing,
                     include.lowest = T))                                       # 1 from the top
  }
  if (is.null(breaks_partition$breaks_easting)) {
    coords_all_ptt$col <- as.numeric(cut_interval(coords_all_ptt$easting,
                                                  n = n_easting,
                                                  labels = 1:n_easting))
  } else {
    coords_all_ptt$col <- as.numeric(cut(coords_all_ptt$easting,
                                         breaks = breaks_easting,
                                         labels = 1:n_easting,
                                         include.lowest = T))
  }
  if (is.null(breaks_partition$breaks_time)) {
    coords_all_ptt$time_d <-
      as.numeric(as.character(coords_all_ptt$time*(n_time-1) + 1))
  } else {
    coords_all_ptt$time_d <- as.numeric(cut(coords_all_ptt$time,
                                            breaks = breaks_time,
                                            labels = 1:n_time,
                                            include.lowest = T))
  }

  coords_all_ptt$partition <- sprintf(bags:::ptt_format(nd),
                                      coords_all_ptt$row,
                                      coords_all_ptt$col,
                                      coords_all_ptt$time_d)

  coords_ptt <- coords_all_ptt[1:n,]
  coords_pred_ptt <- coords_all_ptt[-c(1:n),]

  # index set
  ptts <- sort(unique(coords_ptt$partition))
  idx <- list()
  for (m in ptts) {
    idx[[m]] <- which(coords_ptt$partition == m)
  }

  # for prototypical partitions
  ptts_proto <- grep(paste0("\\d+,\\d+,", sprintf(paste0("%0", nd, ".0f"), 1),
                            "|\\d+,\\d+,", sprintf(paste0("%0", nd, ".0f"), 2)),
                     ptts, value = TRUE)

  # parent partition and idx
  pptts_wind = pptts_list <- list()
  for (m in ptts) {
    # level 1: partition
    tmp <- data.frame(x = m) %>%
      tidyr::separate(x,
                      c("row", "col", "time"),
                      sep = ",", convert = TRUE)

    pptts_wind[1:K] <- list(NULL)                                               # necessary for Rcpp
    names(pptts_wind) <- directions

    for (h in directions) {
      # level 2: wind direction
      # ppttstmp <- c(parentS(tmp, coords_ptt, h, nd),
      #               parentT(tmp, coords_ptt, nd))
      ppttstmp <- c(bags:::parentS(tmp, coords_ptt, h, nd),
                    parentST(tmp, coords_ptt, h, nd))

      if (sum(is.na(ppttstmp)) < 2) {                                           # $ppartition and $idx exist only when there is at least one parents
        pptts_wind_inf <- list()
        pptts_wind_inf[['ppartition']] <- ppttstmp
        pptts_wind_inf[['pidx']] <- c(idx[[ppttstmp[1]]],
                                      idx[[ppttstmp[2]]])                       # important to keep ppartition order
        pptts_wind[[h]] <- pptts_wind_inf
      }
    }
    pptts_list[[m]] <- pptts_wind
  }

  if (save_data) {
    out$coords_ptt <- coords_all_ptt
  }

  ##################
  # initial values #
  ##################
  if (is.null(init$tau_sq)) { tau_sq <- 0.1 }  else { tau_sq <- init$tau_sq }
  if (is.null(init$sig_sq)) { sig_sq <- 1 }  else { sig_sq <- init$sig_sq }
  if (is.null(init$w)) { w <- rnorm(n) }  else { w <- init$w }

  if (debug$z_fixed) {
    if (is.null(init$z)) {
      stop("z should be given for a non-mixture model.")
    } else z <- init$z
  } else {
    if (is.null(init$z)) {
      z <- rep(directions[1], length(ptts))
    }  else { z <- init$z }
    pi_i <- matrix(1/K, ncol = K, nrow = length(z))
  }

  if (debug$psi_fixed) {
    if (is.null(init$psi)) {
      stop("psi should be given for a model with fixed psi.")
    } else psi <- init$psi
  } else {
    if (is.null(init$psi)) { psi <- rep(0.5, 3) }  else { psi <- init$psi }
    if (is.null(init$Sn)) { Sn <- diag(1, 3) }  else { Sn <- init$Sn }
  }

  ###################
  # hyperparameters #
  ###################
  if (is.null(hyper$at)) { at <- 2 } else { at <- hyper$at }
  if (is.null(hyper$bt)) { bt <- 0.1 } else { bt <- hyper$bt }
  if (is.null(hyper$as)) { as <- 2 } else { as <- hyper$as }
  if (is.null(hyper$bs)) { bs <- 1 } else { bs <- hyper$bs }
  if (!is.null(X)) {
    if (is.null(hyper$mu0)) { mu0 <- rep(0, p) } else { mu0 <- hyper$mu0 }
    if (is.null(hyper$invV0)) {
      invV0 <- diag(0.01, p)
    } else { invV0 <- hyper$invV0 }
  }

  if (!debug$psi_fixed) {
    if (is.null(hyper$la)) {
      la <- 0.001; ua <- 1000
    } else {
      la <- hyper$la; ua <- hyper$ua
    }
    if (is.null(hyper$lc)) {
      lc <- 0.001; uc <- 1000
    } else {
      lc <- hyper$lc; uc <- hyper$uc
    }

    # update initial values using bound info
    if (is.null(init$psi)) {
      psi[1] <- runif(1, la, ua)
      psi[2] <- runif(1, lc, uc)
    }

    # use invS to store unif bound info
    invS <- matrix(c(la, ua, lc, uc), byrow = TRUE, 2, 2)

    theta <- bags:::logloglogit(psi, unifprior = TRUE,
                         la = la, ua = ua,
                         lc = lc, uc = uc)
  }

  ################
  # store output #
  ################
  out$tau_sq_save <- rep(0, save)
  out$sig_sq_save <- rep(0, save)
  if (!is.null(X)) out$beta_save <- matrix(0, p, save)
  if (!debug$z_fixed) {
    out$z_save <- matrix(0, length(ptts), save)
  }
  if (!debug$psi_fixed) out$psi_save <- matrix(0, 3, save)
  w_save <- matrix(0, n, save)
  if (save_est) {
    out$y_save <- matrix(0, n, save)
  }

  RnH_list <- bags:::createRnH(coords = coords, idx = idx,
                        pptts_list = pptts_list,
                        ptts_proto = ptts_proto,
                        a = psi[1], c = psi[2], kappa = psi[3],
                        produce_R = FALSE, directions)
  res <- y

  ##############
  # estimation #
  ##############
  if (verbose) {
    cat('burnin...\n')
    pb = txtProgressBar(style=3,width=50)
    pb2 = txtProgressBar(style=3,width=50)
  }

  est_time <- 0
  for (s in 1:S) {
    time_est <- system.time({
      # z update
      if (!debug$z_fixed) {
        z <- bags:::zupdate(pi_i = pi_i, w = w, sig_sq = sig_sq,
                     RnH_list = RnH_list, ptts = ptts, idx = idx,
                     pptts_list = pptts_list, ptts_proto = ptts_proto,
                     nd = nd, directions = directions)
      }

      # w and sig_sq update
      wsig <- bags:::wsigupdate(res = res, w = w, sig_sq = sig_sq,
                         tau_sq = tau_sq, z = z,
                         RnH_list = RnH_list, ptts = ptts, idx = idx,
                         pptts_list = pptts_list, ptts_proto = ptts_proto,
                         nd = nd, as = as, bs = bs)
      w <- as.numeric(wsig$w)
      sig_sq <- wsig$sig_sq

      # psi update
      if (!debug$psi_fixed) {
        psiall <- bags:::psiupdate(w = w, z = z, sig_sq = sig_sq,
                            theta = theta, Sn = Sn, invS = invS,
                            RnH_list = RnH_list, coords = coords,
                            idx = idx, pptts_list = pptts_list,
                            ptts_proto = ptts_proto, ptts = ptts,
                            nd = nd, iter = s,
                            adaptive = TRUE, unifprior = TRUE,
                            directions = directions)
        psi <- psiall$psi
        theta <- psiall$theta
        Sn <- psiall$Sn
        RnH_list <- psiall$RnH_list
      }

      # beta update
      if (!is.null(X)) {
        beta <- bags:::beta_update(y = y, X = X, w = w,
                            tau_sq = tau_sq, mu0 = mu0, invV0 = invV0)
        res <- y - X%*%beta
      }

      # tau_sq update
      tau_sq <- bags:::tau_sq_update(res = res, w = w, at = at, bt = bt)
    })

    est_time <- est_time + time_est

    # save results
    if ((s > burn) & (s-burn) %% thin == 0) {
      out$tau_sq_save[(s-burn)/thin] <- tau_sq
      out$sig_sq_save[(s-burn)/thin] <- sig_sq
      if (!debug$psi_fixed) out$psi_save[,(s-burn)/thin] <- psi
      if (!debug$z_fixed) {
        out$z_save[,(s-burn)/thin] <- z
      }
      if (!is.null(X)) out$beta_save[,(s-burn)/thin] <- beta
      w_save[,(s-burn)/thin] <- w
      if (save_est) {
        if (is.null(X)) {
          out$y_save[,(s-burn)/thin] <- ybar + w + sqrt(tau_sq)*rnorm(n)
        } else {
          out$y_save[,(s-burn)/thin] <- ybar + X%*%beta + w + sqrt(tau_sq)*rnorm(n)
        }
      }
    }

    if (verbose) {
      setTxtProgressBar(pb, s/burn)
      if (s > burn) {
        if (s == burn + 1) {
          close(pb)
          cat('saving...\n')
          setTxtProgressBar(pb2,(s-burn)/(S-burn))
        } else setTxtProgressBar(pb2,(s-burn)/(S-burn))
      }
    }
  }

  if (save_est) {
    out$w_save <- w_save
  }
  out$est_iter <- S
  out$est_time <- est_time

  ##############
  # prediction #
  ##############
  if (verbose) {
    cat('\nprediction...\n')
  }

  # prediction
  ptts_pred <- sort(unique(coords_pred_ptt$partition))
  idx_pred <- list()
  for (l in ptts_pred) {
    idx_pred[[l]] <- which(coords_pred_ptt$partition == l)
  }

  # parent partitions and idx
  z_postm <- apply(out$z_save, 1, getmode)
  pptts_pred_list <- list()
  for (l in ptts_pred) {
    tmp <- data.frame(x = l) %>%
      tidyr::separate(x, c("row", "col", "time"),
                      sep = ",", convert = TRUE)
    if (l %in% ptts) {
      h <- z_postm[which(ptts == l)]
      pptts <- c(l,
                 (pptts_list[[l]])[[h]]$ppartition,                          # parent chosen by wind
                 sprintf(bags:::ptt_format(nd), tmp$row, tmp$col, tmp$time+1))         # future
    } else {
      pptts <- parent8NN(tmp, coords_ptt, nd)
    }
    pptts_idx <- NULL
    for (i in 1:length(pptts)) {
      pptts_idx <- c(pptts_idx, idx[[pptts[i]]])
    }
    pptts_pred_list[[l]] <- list(ppartition = pptts, pidx = pptts_idx)
  }

  # sample w_pred and y_pred from posterior predictive dist
  pred_time <- system.time({
    pred <- bags:::bag_predict(coords_tr = coords,
                        coords_pred = coords_pred,
                        ptts_pred = ptts_pred,
                        idx_pred = idx_pred,
                        pptts_pred_list = pptts_pred_list,
                        w_save = w_save,
                        psi_save = out$psi_save,
                        sigsq_save = out$sig_sq_save,
                        tausq_save = out$tau_sq_save,
                        verbose = verbose,
                        num_threads = n_threads)
  })

  out$w_pred_save <- pred$w_pred_save
  if (is.null(X_pred)) {
    out$y_pred_save <- pred$y_pred_save + ybar
  } else {
    out$y_pred_save <- pred$y_pred_save + ybar + X_pred%*%out$beta_save
  }

  out$pred_iter <- mcmc$save
  out$pred_time <- pred_time

  if (verbose) { close(pb2) }
  return(out)
}
