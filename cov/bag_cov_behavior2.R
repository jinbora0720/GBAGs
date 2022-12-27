# Proposition 5

library(tidyverse)
theme_set(theme_bw())

# spatial w/ exponential covariance 
N <- 100
sig_sq <- 1
phi <- 0.5
h <- seq(0, 10, length = N)
C <- sig_sq*exp(-phi*h)

p <- seq(0, 1, by = 0.1)
Ctilde12 = Ctilde13 <- matrix(0, nrow = N, ncol = length(p))
for (i in 1:length(p)) {
  Ctilde12[,i] <- p[i]*C
  Ctilde13[,i] <- (1-p[i])*C
}

data.frame(dist = rep(h, times = length(p)), 
           p = rep(paste0("p = ", p), each = N), 
           nonstCov12 = as.numeric(Ctilde12), 
           nonstCov13 = as.numeric(Ctilde13)) %>% 
  pivot_longer(-c("dist", "p"), names_to = "what", values_to = "Cov") %>% 
  filter(p %in% paste0("p = ", c(0, 0.2, 0.7, 0.9))) %>% 
  ggplot() + 
  facet_grid(~p) + 
  geom_line(aes(dist, Cov, color = what, linetype = what)) +
  scale_color_manual(values = c("nonstCov12" = "#E69F00", 
                                "nonstCov13" = "#009E73"), 
                     labels = c(expression(paste(s[1],"&",s[2])), 
                                expression(paste(s[1],"&",s[3])))) +
  scale_linetype_manual(values = c("nonstCov12" = 1, 
                                   "nonstCov13" = 2), 
                        labels = c(expression(paste(s[1],"&",s[2])), 
                                   expression(paste(s[1],"&",s[3])))) + 
  labs(x = "Distance", y = "G-BAG induced covariance", 
       color = "", linetype = "") 

# spatiotemporal w/ Apanasovich and Genton (2010) 
a <- 0.07
c <- 0.01
kappa <- 1
sig_sq <- 1
h <- seq(0, 12, length = N)
u <- seq(0, 12, length = N)
hu <- expand.grid(h, u)

basedata <- data.frame(h = hu[,1], u = hu[,2]) %>% 
  mutate(invaup1 = 1/((a*u^2+1)^kappa)) %>% 
  mutate(z = sig_sq*invaup1/((a*h^2+1)^(kappa/2))*exp(-c*h^2*invaup1-c*u^2)) %>% 
  select(-invaup1)

basedata %>% 
  mutate(nonstCov12 = 0.9*z,
         nonstCov13 = (1-0.9)*z) %>% 
  pivot_longer(c("nonstCov12", "nonstCov13"), 
               names_to = "what", values_to = "Cov") %>% 
  ggplot() + 
  geom_contour(aes(x = h, y = u, z = Cov, color = what, linetype = what)) +
  scale_color_manual(values = c("nonstCov12" = "#E69F00", 
                                "nonstCov13" = "#009E73"), 
                     labels = c(expression(paste(s[1],"&",s[2])), 
                                expression(paste(s[1],"&",s[3])))) +
  scale_linetype_manual(values = c("nonstCov12" = 1, 
                                   "nonstCov13" = 2), 
                        labels = c(expression(paste(s[1],"&",s[2])), 
                                   expression(paste(s[1],"&",s[3])))) + 
  labs(x = "h", y = "u", title = "G-BAG induced covariance", 
       color = "", linetype = "") 

p <- c(0.51, 0.55, seq(0.6, 0.9, by = 0.1), 0.95, 0.99)
Ctilde12 = Ctilde13 <- matrix(0, nrow = N^2, ncol = length(p))
for (i in 1:length(p)) {
  Ctilde12[,i] <- p[i]*basedata$z
  Ctilde13[,i] <- (1-p[i])*basedata$z
}

data.frame(h = rep(basedata$h, times = length(p)), 
           u = rep(basedata$u, times = length(p)), 
           p = rep(paste0("p = ", p), each = N^2), 
           nonstCov12 = as.numeric(Ctilde12), 
           nonstCov13 = as.numeric(Ctilde13)) %>% 
  pivot_longer(c("nonstCov12", "nonstCov13"), 
               names_to = "what", values_to = "Cov") %>% 
  filter(p %in% paste0("p = ", c(0.51, 0.55, 0.7))) %>% 
  ggplot() + 
  facet_wrap(~p, ncol = 3) + 
  geom_contour(aes(x = h, y = u, z = Cov, color = what, linetype = what), 
               breaks = 0.1) +
  scale_color_manual(values = c("nonstCov12" = "#E69F00", 
                                "nonstCov13" = "#009E73"), 
                     labels = c(expression(paste(s[1],"&",s[2])), 
                                expression(paste(s[1],"&",s[3])))) +
  scale_linetype_manual(values = c("nonstCov12" = 1, 
                                   "nonstCov13" = 2), 
                        labels = c(expression(paste(s[1],"&",s[2])), 
                                   expression(paste(s[1],"&",s[3])))) + 
  labs(x = "||h||", y = "|u|", 
       title = "G-BAG induced covariance: nonstationarity", 
       color = "", linetype = "") +
  theme(legend.position = "bottom")
# for (ext in extension) {
#   ggsave(plot = ggg, paste0(path, "plots/cov_time", ext),
#          width = 11, height = 3)
# }