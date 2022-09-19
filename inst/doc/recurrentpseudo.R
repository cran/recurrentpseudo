## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning = F, error = T, message = F--------------------------------------
# Development version
#require(devtools)
#devtools::install_github("JulieKFurberg/recurrentpseudo", force = TRUE)

require(recurrentpseudo)

## -----------------------------------------------------------------------------
# Main functions
#?pseudo.onedim
#?pseudo.twodim
#?pseudo.threedim

#?pseudo.geefit


## ----cache = F----------------------------------------------------------------
# Example: Bladder cancer data from survival package
require(survival)

# Make a three level status variable
bladder1$status3 <- ifelse(bladder1$status %in% c(2, 3), 2, bladder1$status)

# Add one extra day for the two patients with start=stop=0
# subset(bladder1, stop <= start)
bladder1[bladder1$id == 1, "stop"] <- 1
bladder1[bladder1$id == 49, "stop"] <- 1

# Restrict the data to placebo and thiotepa
bladdersub <- subset(bladder1, treatment %in% c("placebo", "thiotepa"))

# Make treatment variable two-level factor
bladdersub$Z <- as.factor(ifelse(bladdersub$treatment == "placebo", 0, 1))
levels(bladdersub$Z) <- c("placebo", "thiotepa")
head(bladdersub)


## ----cache = F----------------------------------------------------------------
# Pseudo observations at t = 30
pseudo_bladder_1d <- pseudo.onedim(tstart = bladdersub$start,
                                   tstop = bladdersub$stop,
                                   status = bladdersub$status3,
                                   id = bladdersub$id,
                                   covar_names = "Z",
                                   tk = c(30),
                                   data = bladdersub)
head(pseudo_bladder_1d$outdata)

# GEE fit
fit_bladder_1d <- pseudo.geefit(pseudodata = pseudo_bladder_1d,
                                covar_names = c("Z"))
fit_bladder_1d

# Treatment differences
xi_diff_1d <- as.matrix(c(fit_bladder_1d$xi[2]), ncol = 1)

mslabels <- c("treat, mu")
rownames(xi_diff_1d) <- mslabels
colnames(xi_diff_1d) <- ""
xi_diff_1d

# Variance matrix for differences
sigma_diff_1d <- matrix(c(fit_bladder_1d$sigma[2,2]),
                          ncol = 1, nrow = 1,
                          byrow = T)

rownames(sigma_diff_1d) <- colnames(sigma_diff_1d) <- mslabels
sigma_diff_1d

## ----cache = F----------------------------------------------------------------
# Pseudo observations at t = 30
pseudo_bladder_2d <- pseudo.twodim(tstart = bladdersub$start,
                                   tstop = bladdersub$stop,
                                   status = bladdersub$status3,
                                   id = bladdersub$id,
                                   covar_names = "Z",
                                   tk = c(30),
                                   data = bladdersub)
head(pseudo_bladder_2d$outdata)

# GEE fit
fit_bladder_2d <- pseudo.geefit(pseudodata = pseudo_bladder_2d,
                                covar_names = c("Z"))
fit_bladder_2d

# Treatment differences
xi_diff_2d <- as.matrix(c(fit_bladder_2d$xi[2],
                          fit_bladder_2d$xi[4]), ncol = 1)

mslabels <- c("treat, mu", "treat, surv")
rownames(xi_diff_2d) <- mslabels
colnames(xi_diff_2d) <- ""
xi_diff_2d


# Variance matrix for differences
sigma_diff_2d <- matrix(c(fit_bladder_2d$sigma[2,2],
                          fit_bladder_2d$sigma[2,4],
                          fit_bladder_2d$sigma[2,4],
                          fit_bladder_2d$sigma[4,4]),
                          ncol = 2, nrow = 2,
                          byrow = T)

rownames(sigma_diff_2d) <- colnames(sigma_diff_2d) <- mslabels
sigma_diff_2d

## ----eval = T-----------------------------------------------------------------
# Add deathtype variable to bladder data
# Deathtype = 1 (bladder disease death), deathtype = 2 (other death reason)
bladdersub$deathtype <- with(bladdersub, ifelse(status == 2, 1, ifelse(status == 3, 2, 0)))
table(bladdersub$deathtype, bladdersub$status)

# Pseudo-observations
pseudo_bladder_3d <- pseudo.threedim(tstart = bladdersub$start,
                                     tstop = bladdersub$stop,
                                     status = bladdersub$status3,
                                     id = bladdersub$id,
                                     deathtype = bladdersub$deathtype,
                                     covar_names = "Z",
                                     tk = c(30), 
                                     data = bladdersub)
head(pseudo_bladder_3d$outdata_long)

# GEE fit
fit_bladder_3d <- pseudo.geefit(pseudodata = pseudo_bladder_3d,
                                covar_names = c("Z"))
fit_bladder_3d

# Treatment differences
xi_diff_3d <- as.matrix(c(fit_bladder_3d$xi[2],
                          fit_bladder_3d$xi[4],
                          fit_bladder_3d$xi[6]), ncol = 1)

mslabels <- c("treat, mu", "treat, cif1", "treat, cif2")
rownames(xi_diff_3d) <- mslabels
colnames(xi_diff_3d) <- ""
xi_diff_3d


# Variance matrix for differences
sigma_diff_3d <- matrix(c(fit_bladder_3d$sigma[2,2],
                          fit_bladder_3d$sigma[2,4],
                          fit_bladder_3d$sigma[2,6],

                          fit_bladder_3d$sigma[2,4],
                          fit_bladder_3d$sigma[4,4],
                          fit_bladder_3d$sigma[4,6],

                          fit_bladder_3d$sigma[2,6],
                          fit_bladder_3d$sigma[4,6],
                          fit_bladder_3d$sigma[6,6]

                          ),
                        ncol = 3, nrow = 3,
                        byrow = T)

rownames(sigma_diff_3d) <- colnames(sigma_diff_3d) <- mslabels
sigma_diff_3d


## ----eval = T-----------------------------------------------------------------
# Compare - should match for mu elements 
xi_diff_1d
xi_diff_2d
xi_diff_3d

sigma_diff_1d
sigma_diff_2d
sigma_diff_3d

## ----eval = T-----------------------------------------------------------------
require(dplyr)

## One-dim
# A binary variable, Z1_
# A continuous variable, Z2_
# A categorical variable, Z3_
set.seed(0308)

bladdersub <- as.data.frame(
  bladdersub %>% group_by(id) %>% 
  mutate(Z1_ = Z,
         Z2_ = rnorm(1, mean = 3, sd = 1),
         Z3_ = sample(x = c("A", "B", "C"), 
                      size = 1, replace = TRUE, 
                      prob = c(1/4, 1/2, 1/4))
         ))
# head(bladdersub, 20)

# Make pseudo obs at more timepoints (more data)
# Pseudo observations at t = 20, 30, 40
pseudo_bladder_1d_3t <- pseudo.onedim(tstart = bladdersub$start,
                                      tstop = bladdersub$stop,
                                      status = bladdersub$status3,
                                      id = bladdersub$id,
                                      covar_names = c("Z1_", "Z2_", "Z3_"),
                                      tk = c(20, 30, 40),
                                      data = bladdersub)

fit1 <- pseudo.geefit(pseudodata = pseudo_bladder_1d_3t, 
                      covar_names = c("Z1_", "Z2_", "Z3_"))

fit1$xi
fit1$sigma

fit1$xi[4]

## ----eval = T-----------------------------------------------------------------
## Two-dim
# Pseudo observations at t = 20, 30, 40
pseudo_bladder_2d_3t <- pseudo.twodim(tstart = bladdersub$start,
                                      tstop = bladdersub$stop,
                                      status = bladdersub$status3,
                                      id = bladdersub$id,
                                      covar_names = c("Z1_", "Z2_", "Z3_"),
                                      tk = c(20, 30, 40),
                                      data = bladdersub)

fit2 <- pseudo.geefit(pseudodata = pseudo_bladder_2d_3t, 
                      covar_names = c("Z1_", "Z2_", "Z3_"))

# fit2$xi
# fit2$sigma


## Three-dim
pseudo_bladder_3d_3t <- pseudo.threedim(tstart = bladdersub$start,
                                        tstop = bladdersub$stop,
                                        status = bladdersub$status3,
                                        id = bladdersub$id,
                                        covar_names = c("Z1_", "Z2_", "Z3_"),
                                        deathtype = bladdersub$deathtype,
                                        tk = c(20, 30, 40),
                                        data = bladdersub)

fit3 <- pseudo.geefit(pseudodata = pseudo_bladder_3d_3t, 
                      covar_names = c("Z1_", "Z2_", "Z3_"))

# fit3$xi
# fit3$sigma

## Compare for mu
fit1$xi[4]
fit2$xi[4]
fit3$xi[4]

