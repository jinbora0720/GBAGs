# Fitted G-BAG is correctly specified #
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
# theta2
bag <- readRDS(paste0(path, "sim1/sim1b_bag.RDS"))
bag_new <- readRDS(paste0(path, "sim1/sim1b_bag_newparent.RDS"))

#########
# Table #
#########
# mean
mtab <-
  rbind(cbind(sapply(bag[c("beta_hat", "tausq_hat", "sigsq_hat")], mean),
              sapply(bag_new[c("beta_hat", "tausq_hat", "sigsq_hat")], mean)),
        cbind(rowMeans(bag[["psi_hat"]]),
              rowMeans(bag_new[["psi_hat"]])),
        cbind(sapply(bag[c("rmspe", "mape", "coverage", "meanwidth")], mean),
              sapply(bag_new[c("rmspe", "mape", "coverage", "meanwidth")], mean)))
rownames(mtab)[c(1:6)] <- c("beta", "tausq", "sigsq", "a", "c", "kappa")
colnames(mtab) <- c("G-BAG", "G-BAG (new)")
mtab %>% round(3)

# sd
sdtab <-
  rbind(cbind(sapply(bag[c("beta_hat", "tausq_hat", "sigsq_hat")], sd),
              sapply(bag_new[c("beta_hat", "tausq_hat", "sigsq_hat")], sd)),
        cbind(apply(bag[["psi_hat"]], 1, sd),
              apply(bag_new[["psi_hat"]], 1, sd)),
        cbind(sapply(bag[c("rmspe", "mape", "coverage", "meanwidth")], sd),
              sapply(bag_new[c("rmspe", "mape", "coverage", "meanwidth")], sd)))
rownames(sdtab)[c(1:6)] <- c("beta", "tausq", "sigsq", "a", "c", "kappa")
colnames(sdtab) <- c("G-BAG", "G-BAG (new)")
sdtab %>% round(3)
