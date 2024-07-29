parentST <- function(coords, coords_ptt,
                    wind_dir = c("W", "NW", "S", "SW", "N", "NE", "E", "SE"),
                    nd) {
  # coords should have one row, col, and time
  # coords_ptt should have row, col, and time
  
  if (wind_dir == "W") {
    prow <- coords$row
    candd <- coords_ptt[coords_ptt$row == prow, c("row","col")]
    if (sum(candd$col < coords$col) == 0) {
      pcol <- NA
    } else {
      pcol <- max(candd$col[candd$col < coords$col])
    }
  }
  
  if (wind_dir == "E") {
    prow <- coords$row
    candd <- coords_ptt[coords_ptt$row == prow, c("row","col")]
    if (sum(candd$col > coords$col) == 0) {
      pcol <- NA
    } else {
      pcol <- min(candd$col[candd$col > coords$col])
    }
  }
  
  if (wind_dir == "NW") {
    diff <- coords$row - coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) -
                          as.numeric(coords_ptt$col) == diff, c("row", "col")]
    if (sum(candd$col < coords$col)*sum(candd$row < coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- max(candd$row[candd$row < coords$row])
      pcol <- max(candd$col[candd$col < coords$col])
    }
  }
  
  if (wind_dir == "SE") {
    diff <- coords$row - coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) -
                          as.numeric(coords_ptt$col) == diff, c("row", "col")]
    if (sum(candd$col > coords$col)*sum(candd$row > coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- min(candd$row[candd$row > coords$row])
      pcol <- min(candd$col[candd$col > coords$col])
    }
  }
  
  if (wind_dir == "S") {
    pcol <- coords$col
    candd <- coords_ptt[coords_ptt$col == pcol, c("row","col")]
    if (sum(candd$row > coords$row) == 0) {
      prow <- NA
    } else {
      prow <- min(candd$row[candd$row > coords$row])
    }
  }
  
  if (wind_dir == "N") {
    pcol <- coords$col
    candd <- coords_ptt[coords_ptt$col == pcol, c("row","col")]
    if (sum(candd$row < coords$row) == 0) {
      prow <- NA
    } else {
      prow <- max(candd$row[candd$row < coords$row])
    }
  }
  
  if (wind_dir == "SW") {
    total <- coords$row + coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) +
                          as.numeric(coords_ptt$col) == total, c("row", "col")]
    if (sum(candd$col < coords$col)*sum(candd$row > coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- min(candd$row[candd$row > coords$row])
      pcol <- max(candd$col[candd$col < coords$col])
    }
  }
  
  if (wind_dir == "NE") {
    total <- coords$row + coords$col
    candd <- coords_ptt[as.numeric(coords_ptt$row) +
                          as.numeric(coords_ptt$col) == total, c("row", "col")]
    if (sum(candd$col > coords$col)*sum(candd$row < coords$row) == 0) {
      pcol <- NA; prow <- NA
    } else {
      prow <- max(candd$row[candd$row < coords$row])
      pcol <- min(candd$col[candd$col > coords$col])
    }
  }
  
  if (is.null(coords_ptt$time_d)) {
    candd <- coords_ptt[which(coords_ptt$row == prow &
                                coords_ptt$col == pcol), "time"]
  } else {
    candd <- coords_ptt[which(coords_ptt$row == prow &
                                coords_ptt$col == pcol), "time_d"]
  }
  if (sum(candd < coords$time) == 0) {
    ptime <- NA
  } else{
    ptime <- max(candd[candd < coords$time])
  }
  
  return(ifelse(is.na(prow) + is.na(pcol) + is.na(ptime) > 0, NA,
                sprintf(bags:::ptt_format(nd), prow, pcol, ptime)))
}

Ctilde_ST <- function(coords_ptt,
                   z,
                   params = list(sig_sq = 2,
                                 a = 5,
                                 c = 0.5,
                                 kappa = 0.9)) {
  ##################
  # to get started #
  ##################
  # data size
  n <- nrow(coords_ptt)
  H = R <- matrix(0, n, n)
  
  # select directions
  directions <- unique(as.vector(z))
  K <- length(directions)
  
  # index set
  nd <- nchar(strsplit(coords_ptt$partition[1], ",")[[1]][1])
  ptts <- sort(unique(coords_ptt$partition))
  idx <- list()
  for (m in ptts) {
    idx[[m]] <- which(coords_ptt$partition == m)
  }
  
  # assign rownames
  if (is.null(rownames(z))) {
    rownames(z) <- ptts
  }
  
  # parent partition and idx
  pptts_wind = pptts_list <- list()
  for (m in ptts) {
    # level 1: partition
    tmp <- data.frame(x = m) %>%
      tidyr::separate(x, c("row", "col", "time"),
                      sep = ",", convert = TRUE)
    
    pptts_wind[1:K] <- list(NULL)
    names(pptts_wind) <- directions
    
    for (h in directions) {
      # level 2: wind direction
      ppttstmp <- c(
        bags:::parentS(tmp, coords_ptt, h, nd),
        # bags:::parentT(tmp, coords_ptt, nd)
        parentST(tmp, coords_ptt, h, nd)
        )
      # $ppartition and $idx exist only when there is at least one parents
      if (sum(is.na(ppttstmp)) < 2) {
        pptts_wind_inf <- list()
        pptts_wind_inf[['ppartition']] <- ppttstmp
        pptts_wind_inf[['pidx']] <- c(idx[[ppttstmp[1]]],
                                      idx[[ppttstmp[2]]])                       # change: important to keep ppartition order
        pptts_wind[[h]] <- pptts_wind_inf
      }
    }
    pptts_list[[m]] <- pptts_wind
  }
  
  a <- params$a
  c <- params$c
  kappa <- params$kappa
  sig_sq <- params$sig_sq
  
  ##########
  # Ctilde #
  ##########
  for (l in ptts) {
    idx_l <- idx[[l]]
    n_l <- length(idx_l)
    h_l <- z[l,]
    if (length(h_l) > 1) {
      idx_pl <- NULL
      for (s in 1:length(h_l)) {
        idx_pl <- c(idx_pl, pptts_list[[l]][[h_l[s]]]$pidx)
      }
    } else {
      idx_pl <- pptts_list[[l]][[h_l]]$pidx
    }
    coords_tmp <- rbind(coords_ptt[idx_l, c("easting", "northing", "time")],
                        coords_ptt[idx_pl, c("easting", "northing", "time")])
    spatdist <- dist(coords_tmp[,c("easting", "northing")], method = "euclidean")
    timedist <- dist(coords_tmp[,c("time")], method = "manhattan")
    invaup1 <- 1/(a*as.matrix(timedist)+1)
    Cor <- invaup1*exp(-c*as.matrix(spatdist)*(invaup1^(kappa/2)))
    if (is.null(idx_pl)) {
      R[idx_l, idx_l] <- sig_sq*Cor[1:n_l, 1:n_l]
    } else {
      CinvC <- Cor[1:n_l, -c(1:n_l)] %*%
        solve(Cor[-c(1:n_l), -c(1:n_l)] + diag(0.0005, length(idx_pl)))
      R[idx_l, idx_l] <- sig_sq*(Cor[1:n_l, 1:n_l] - CinvC%*%Cor[-c(1:n_l), 1:n_l])
      H[idx_l, idx_pl] <- CinvC
    }
  }
  
  invImH <- solve(diag(1, n) - H)
  return(invImH %*% R %*% t(invImH))
}
