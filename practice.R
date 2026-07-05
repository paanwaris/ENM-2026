# ==========================================================================
# ENM Practice Script  —  nicheR (niches) + bean (bias thinning)
# ==========================================================================
# Run it top to bottom. Four parts:
#
#   PART 1     Virtual niches with nicheR  — generalist vs specialist (2-D)
#   PART 2     Burn (fire) occurrences — niche WITHOUT thinning
#   PART 3     Sambar deer (bean) — thin, then compare ORIGINAL vs THINNED
#   YOUR TURN  Template to import your own species from GBIF
#
# Two variables are used throughout to keep things simple:
#   bio_1  = mean annual temperature
#   bio_12 = annual precipitation
# ==========================================================================

# install.packages(c("terra", "nicheR", "bean"))   # run once, if needed
library(terra)
library(bean)
library(nicheR)

# Environmental layers for Thailand (these ship with the bean package):
env <- rast(system.file("extdata", "thai_env.tif", package = "bean"))
env                      # a quick look: layer names, extent, resolution


# ==========================================================================
# PART 1  —  Virtual niches with nicheR  (2-D: bio_1, bio_12)
# ==========================================================================
# A GENERALIST tolerates a WIDE range of conditions; a SPECIALIST only a
# NARROW range. We describe each niche by its min/max tolerance ranges.

generalist <- build_ellipsoid(
  range = data.frame(bio_1  = c(21, 28),      # wide temperature
                     bio_12 = c(800, 1500)),  # wide precipitation
  cl = 0.95
)

specialist <- build_ellipsoid(
  range = data.frame(bio_1  = c(24, 27),       # narrow: only warm
                     bio_12 = c(1500, 2500)),  # narrow: only moderate rain
  cl = 0.95
)

# Plot the two niche ellipses in E-space (bio_1 vs bio_12), over the Thai
# background cloud. plot_ellipsoid() opens the plot; add_ellipsoid() overlays.
back_data <- as.data.frame(env, xy = TRUE)      # background = every pixel's climate
plot_ellipsoid(generalist,
               background = back_data[, c("bio_1", "bio_12")],
               col_ell = "#5E78C4", col_bg = "grey85",   # blue ellipse, grey cloud
               lwd = 2.5, pch = 20, cex_bg = 0.4,
               xlab = "bio_1 (temperature)", ylab = "bio_12 (precipitation)",
               main = "Generalist (blue) vs Specialist (pink)")
add_ellipsoid(specialist, col_ell = "#D9559E", lwd = 2.5)   # overlay the specialist

# Predict suitability for every pixel (0 = unsuitable, 1 = ideal):
gen_pred  <- predict(generalist, newdata = env[[c("bio_1", "bio_12")]],
                     suitability_truncated = TRUE)
spec_pred <- predict(specialist, newdata = env[[c("bio_1", "bio_12")]],
                     suitability_truncated = TRUE)

# Map them side by side — the generalist lights up a much larger area:
par(mfrow = c(1, 2))
plot(gen_pred[["suitability_trunc"]],  main = "Generalist — suitability")
plot(spec_pred[["suitability_trunc"]], main = "Specialist — suitability")
par(mfrow = c(1, 1))


# ==========================================================================
# PART 2  —  Burn (fire) occurrences, niche WITHOUT thinning
# ==========================================================================
# Fire occurrence points (columns: year, x, y). This is a big file, so the
# read may take a moment.
burn <- read.csv("burn/point_fire.csv")

# See the raw fire points on the map (subsample for a quick, readable plot):
plot(env[["bio_1"]], main = "Fire (burn) occurrences over Thai temperature")
show <- if (nrow(burn) > 5000) burn[sample(nrow(burn), 5000), ] else burn
points(show$x, show$y, pch = 20, cex = 0.3, col = "#c0392b")

# Attach the environment to each fire point. transform = "none" keeps raw units
# so we can predict directly onto the raster (no re-scaling needed):
burn_prep <- prepare_bean(
  data        = burn,
  env_rasters = env,
  longitude   = "x",
  latitude    = "y",
  transform   = "none"
)

# Fit a niche to the fire points WITHOUT thinning, then predict + plot:
burn_fit  <- fit_ellipsoid(burn_prep, env_vars = c("bio_1","bio_4", "bio_12"),
                           method = "covmat", level = 0.95)
plot(burn_fit)
burn_suit <- predict(burn_fit, newdata = env[[c("bio_1", "bio_4", "bio_12")]],
                     suitability_truncated = TRUE)

plot(burn_suit[["suitability_trunc"]],
     main = "Fire (burn) — suitability (WITHOUT thinning)")


# ==========================================================================
# PART 3  —  Sambar deer (bean): thin, then compare ORIGINAL vs THINNED
# ==========================================================================
# Biased occurrence records that ship with bean:
sambar <- read.csv(system.file("extdata", "Rusa_unicolor.csv", package = "bean"))

# Attach the environment (raw units):
sambar_prep <- prepare_bean(
  data        = sambar,
  env_rasters = env,
  longitude   = "x",
  latitude    = "y",
  transform   = "none"
)

# 1) Pick an objective grid size, then thin in environmental space:
res     <- find_env_resolution(sambar_prep, env_vars = c("bio_1", "bio_12"),
                               method = "sheather-jones")
thinned <- thin_env_nd(sambar_prep, env_vars = c("bio_1", "bio_12"),
                       grid_resolution = res$suggested_resolution, seed = 123)

print(thinned)   # how many points were kept?

# 2) Fit a niche to the ORIGINAL and to the THINNED points, then predict both:
fit_original  <- fit_ellipsoid(sambar_prep, env_vars = c("bio_1", "bio_12"),
                               method = "covmat", level = 0.95)
fit_thinned   <- fit_ellipsoid(thinned$thinned_data, env_vars = c("bio_1", "bio_12"),
                               method = "covmat", level = 0.95)

suit_original <- predict(fit_original, newdata = env[[c("bio_1", "bio_12")]],
                         suitability_truncated = TRUE)
suit_thinned  <- predict(fit_thinned, newdata = env[[c("bio_1", "bio_12")]],
                         suitability_truncated = TRUE)

# 3) Compare the two suitability maps: original (biased) vs thinned (corrected):
par(mfrow = c(1, 2))
plot(suit_original[["suitability_trunc"]], main = "Sambar — without thinning")
plot(suit_thinned[["suitability_trunc"]],  main = "Sambar — with thinning")
par(mfrow = c(1, 1))


# ==========================================================================
# YOUR TURN  —  import your own species from GBIF
# ==========================================================================
# Create a folder to store the files you download from GBIF:
dir.create("GBIF", showWarnings = FALSE)

# ---- Uncomment and edit the block below (needs internet + the rgbif package) ----
# install.packages("rgbif")
# library(rgbif)
#
# my_species <- "Rusa unicolor"          # <-- put YOUR species name here
#
# # Download occurrences that have coordinates, and save them into GBIF/:
# hits   <- occ_search(scientificName = my_species,
#                      hasCoordinate  = TRUE,
#                      limit          = 2000)
# my_occ <- hits$data[, c("decimalLongitude", "decimalLatitude")]
# names(my_occ) <- c("x", "y")           # bean expects columns named x / y
# my_occ <- na.omit(my_occ)
# write.csv(my_occ, "GBIF/my_species.csv", row.names = FALSE)
#
# # Run the SAME pipeline as Part 3 on your own data:
# my_prepared <- prepare_bean(my_occ, env_rasters = env,
#                             longitude = "x", latitude = "y", transform = "none")
# my_res      <- find_env_resolution(my_prepared, env_vars = c("bio_1", "bio_12"))
# my_thinned  <- thin_env_nd(my_prepared, env_vars = c("bio_1", "bio_12"),
#                            grid_resolution = my_res$suggested_resolution, seed = 1)
#
# my_fit_o  <- fit_ellipsoid(my_prepared,             env_vars = c("bio_1", "bio_12"))
# my_fit_t  <- fit_ellipsoid(my_thinned$thinned_data, env_vars = c("bio_1", "bio_12"))
# my_suit_o <- predict(my_fit_o, newdata = env[[c("bio_1", "bio_12")]], suitability_truncated = TRUE)
# my_suit_t <- predict(my_fit_t, newdata = env[[c("bio_1", "bio_12")]], suitability_truncated = TRUE)
#
# par(mfrow = c(1, 2))
# plot(my_suit_o[["suitability_trunc"]], main = "Your species — without thinning")
# plot(my_suit_t[["suitability_trunc"]], main = "Your species — with thinning")
# par(mfrow = c(1, 1))
#
# NOTE: if your species lives outside Thailand, swap `env` for climate layers
# that cover its range (e.g. WorldClim), keeping the bio_1 / bio_12 layer names.
# ==========================================================================
