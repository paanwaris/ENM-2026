# ---------------------------------------------------------------------------
# Demo: reshape a niche's rotation & correlation WITHOUT rebuilding it
# ---------------------------------------------------------------------------
# NOTE: update_ellipsoid_covariance() is a function that already ships with the
# nicheR package (R/ellipsoid_utilities.R). This file does NOT redefine it — it
# only shows how to use it. Do not `source()` a competing definition.
#
#   update_ellipsoid_covariance(object, covariance, tol = 1e-6, verbose = TRUE)
#
#   object     : a nicheR_ellipsoid (from build_ellipsoid())
#   covariance : either ONE number applied to every off-diagonal element,
#                or a NAMED vector giving specific variable pairs,
#                e.g. c("bio_1-bio_12" = -100)
#   tol        : tolerance used when checking positive-definiteness limits
#
# The function swaps in the new covariance, then re-runs ellipsoid_calculator()
# so the centroid, semi-axes, eigen-structure, volume and inverse are all
# recomputed consistently. A valid covariance must keep the matrix positive
# definite; check the safe range first with `object$cov_limits`.
# ---------------------------------------------------------------------------

library(nicheR)

## 1. Build a starting (axis-aligned) niche -------------------------------
ell <- build_ellipsoid(
  range = data.frame(bio_1 = c(3, 28), bio_12 = c(200, 3500)),
  cl    = 0.95
)

## 2. Look at the SAFE covariance range for the variable pair --------------
ell$cov_limits            # the min/max covariance that keeps the niche valid

## 3. Introduce a POSITIVE correlation (tilts warm <-> wet together) -------
# A guaranteed-valid magnitude is anything with |cov| < sqrt(var1 * var2):
v         <- diag(ell$cov_matrix)                # the two marginal variances
safe_cov  <- 0.5 * sqrt(v[[1]] * v[[2]])         # correlation ~ 0.5 (well inside limits)

ell_pos <- update_ellipsoid_covariance(ell, covariance =  safe_cov)
ell_neg <- update_ellipsoid_covariance(ell, covariance = -safe_cov)

## 4. Named form — set a specific pair only --------------------------------
ell_named <- update_ellipsoid_covariance(
  ell,
  covariance = c("bio_1-bio_12" = -safe_cov)
)

## 5. What changed? --------------------------------------------------------
ell$volume                 # original
ell_pos$volume             # correlation shrinks the volume
ell_pos$cov_limits_remaining   # remaining safe room for still-zero pairs

## 6. See all three in E-space --------------------------------------------
plot_ellipsoid(ell,     col_ell = "#7f8c8d", lwd = 2.5,
               main = "Reshaping a niche with update_ellipsoid_covariance()")
add_ellipsoid(ell_pos,  col_ell = "#5E78C4", lwd = 2.5)   # tilted up
add_ellipsoid(ell_neg,  col_ell = "#D9559E", lwd = 2.5)   # tilted down
