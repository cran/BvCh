## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5)
set.seed(2024)

## ----load---------------------------------------------------------------------
library(BvCh)

## ----bvch---------------------------------------------------------------------
a1 <- 1.5; a2 <- 2.5; a3 <- 0.5; b <- 1.8
z1 <- c(0.4, 0.9, 0.7, 1.3)
z2 <- c(0.6, 0.4, 0.7, 1.1)
cbind(z1, z2,
      pdf = dbvch(z1, z2, a1, a2, a3, b),
      surv = sbvch(z1, z2, a1, a2, a3, b),
      haz = hbvch(z1, z2, a1, a2, a3, b))

## ----check--------------------------------------------------------------------
S <- sbvch(z1, z2, a1, a2, a3, b)
Sref <- schen(z1, b, a1) * schen(z2, b, a2) * schen(pmax(z1, z2), b, a3)
c(max_abs_diff_surv = max(abs(S - Sref)),
  max_abs_diff_haz  = max(abs(hbvch(z1, z2, a1, a2, a3, b) -
                                dbvch(z1, z2, a1, a2, a3, b) / S)))

## ----mixture------------------------------------------------------------------
K <- a1 + a2 + a3
ac <- dbvch_ac(z1, z2, a1, a2, a3, b)
sg <- dbvch_sing(z1, z2, a1, a2, a3, b)
mix <- (a1 + a2) / K * ac + a3 / K * sg
data.frame(ac = ac, sing = sg, mixture = mix, dbvch = dbvch(z1, z2, a1, a2, a3, b))

## ----mass---------------------------------------------------------------------
m1 <- integrate(function(z) vapply(z, function(v)
  integrate(function(v2) dbvch(v, v2, a1, a2, a3, b), v, Inf,
            rel.tol = 1e-10, subdivisions = 2000L)$value, numeric(1)),
  0, Inf, rel.tol = 1e-8, subdivisions = 2000L)$value
m2 <- integrate(function(z) vapply(z, function(v) if (v > 0)
  integrate(function(v2) dbvch(v, v2, a1, a2, a3, b), 0, v,
            rel.tol = 1e-10, subdivisions = 2000L)$value else 0, numeric(1)),
  0, Inf, rel.tol = 1e-8, subdivisions = 2000L)$value
u <- seq(1e-8, 6, length.out = 200001L)
m0 <- a3 / K * sum(reliaR::dchen(u, b, K)) * (u[2] - u[1])
rbind(f1 = m1, f2 = m2, diagonal = m0, total = m1 + m2 + m0,
      expected = c(a1/K, a2/K, a3/K, 1))

## ----sim----------------------------------------------------------------------
sim <- rbvch(100000, a1, a2, a3, b)
q <- seq(0.4, 1.6, length.out = 6)
data.frame(q = q,
           theo = 1 - schen(q, b, a1) * schen(q, b, a3),
           emp  = sapply(q, function(v) mean(sim[, 1] <= v)))

## ----bec----------------------------------------------------------------------
x1 <- c(0.4, 0.9, 0.7, 1.3); x2 <- c(0.6, 0.4, 0.7, 1.1)
cbind(x1, x2,
      pdf  = dbec(x1, x2, a1, a2, a3, b, lambda = 2),
      surv = sbec(x1, x2, a1, a2, a3, b, lambda = 2),
      haz  = hbec(x1, x2, a1, a2, a3, b, lambda = 2))

## ----becgrad------------------------------------------------------------------
off <- x1 != x2
h12 <- hbec1(x1, x2, a1, a2, a3, b, 2) * hbec2(x1, x2, a1, a2, a3, b, 2)
data.frame(hbec = hbec(x1, x2, a1, a2, a3, b, 2),
           hbec1_times_hbec2 = h12,
           off_diagonal = off)

## ----ifgm---------------------------------------------------------------------
t1 <- 1.2; b1 <- 1.3; t2 <- 1.5; b2 <- 1.1; g <- 0.6; om <- 0.4
cbind(x = 0.5, y = 0.9,
      reliability = relifgmchen(0.5, 0.9, t1, b1, t2, b2, g, om),
      pdf  = difgmchen(0.5, 0.9, t1, b1, t2, b2, g, om),
      haz  = hifgmchen(0.5, 0.9, t1, b1, t2, b2, g, om))

## ----copula-------------------------------------------------------------------
data.frame(C = c(cifgm(0.3, 0.7, g, om), cifgm(0.5, 0.5, g, om)),
           tau = tau_ifgm(g, om), rho = rho_ifgm(g, om))

## ----tauscan------------------------------------------------------------------
data.frame(gamma = c(0, 0.3, 0.6, 0.9, 1),
           tau = sapply(c(0, 0.3, 0.6, 0.9, 1), tau_ifgm, omega = 0.4),
           rho = sapply(c(0, 0.3, 0.6, 0.9, 1), rho_ifgm, omega = 0.4))

## ----fit----------------------------------------------------------------------
d <- rbvch(1000, a1, a2, a3, b)
fit <- fitbvch(d)
fit

## ----se-----------------------------------------------------------------------
coef(fit)
fit$se

## ----ks-----------------------------------------------------------------------
ksfit_bvch(fit, d)

## ----compare------------------------------------------------------------------
f1 <- fitbvch(rbvch(500, a1, a2, a3, b))
f2 <- fitbec(rbec(500, a1, a2, a3, b, 2))
f3 <- fitifgmchen(rifgmchen(500, t1, b1, t2, b2, g, om))
compare_models(f1, f2, f3)

## ----surface, fig.height = 6--------------------------------------------------
plot_surface_bvch(alpha1 = a1, alpha2 = a2, alpha3 = a3,
                  beta.vec = c(0.8, 1.8, 5, 10))

## ----contour------------------------------------------------------------------
plot_contour_bvch(alpha1 = a1, alpha2 = a2, alpha3 = a3, beta = b)

## ----relsurface---------------------------------------------------------------
plot_reliability_ifgmchen(theta1 = t1, beta1 = b1, theta2 = t2, beta2 = b2,
                          gamma = g, omega = om)

## ----ties---------------------------------------------------------------------
c(empirical = mean(sim[, 1] == sim[, 2]), expected = a3 / K)

