# set number of CPU cores for use
Sys.setenv("MC_CORES" = 4)

# load packages
packages <- c("dplyr", "parallel", "rstan")
load <- lapply(packages, "library", character.only = TRUE)

# check number of in-use cores
options("mc.cores")

# read tide gauges, tsunami amplitude data, and covariate data
gauges <- read.csv("tide-gauges.csv")
tsunami.data <- read.csv("tsunami-data.csv")
covariate.data <- read.csv("covariate-data.csv")

# prepare amplitude data
y <- tsunami.data[["Amplitude"]]
N <- nrow(tsunami.data)

# prepare number of events per tide gauge
persite <- (tsunami.data %>% count(Location))
persite <- (persite %>% arrange(desc(Location)))[["n"]]

# prepare threshold of record completeness for each tide gauge
u <- gauges[["Amin"]]

# prepare number of tide gauges
M <- nrow(gauges)

# standardize covariates (i.e., transform to z-scores)
covariate.data[["Latitude"]] <-
  (covariate.data[["Latitude"]] - mean(covariate.data[["Latitude"]])) /
  sd(covariate.data[["Latitude"]])

covariate.data[["Longitude"]] <-
  (covariate.data[["Longitude"]] - mean(covariate.data[["Longitude"]])) /
  sd(covariate.data[["Longitude"]])

covariate.data[["Shelf.Width"]] <-
  (covariate.data[["Shelf.Width"]] - mean(covariate.data[["Shelf.Width"]])) /
  sd(covariate.data[["Shelf.Width"]])

# prepare matrix of covariates / vector of distances
x <- as.matrix(covariate.data[, c(2, 3, 4, 5)])
dist <- covariate.data[, 6]

# compile stan model
lgp.model <- stan_model("TsunamiBHM.stan")

# list Stan input data
feed.data <-
  list(
    y = y,
    N = N,
    persite = persite,
    M = M,
    x = x,
    dist = dist,
    u =  u
  )

# initialize (hyper)parameters at random for all chains
a.init <- list(
  zb = rep(0.25, M),
  zac = rep(0.35, M),
  w0b = 1.5,
  w1b = 1.5,
  w2b = 1.5,
  w3b = 1.5,
  ab = 0.5,
  rhob = 100,
  w0ac = 1.5,
  w1ac = 1.5,
  w2ac = 1.5,
  w3ac = 1.5,
  aac = 0.5,
  rhoac = 100
)

b.init <- list(
  zb = rep(0.35, M),
  zac = rep(0.45, M),
  w0b = 0.5,
  w1b = 0.5,
  w2b = 0.5,
  w3b = 0.5,
  ab = 0.25,
  rhob = 200,
  w0ac = 0.5,
  w1ac = 0.5,
  w2ac = 0.5,
  w3ac = 0.5,
  aac = 0.25,
  rhoac = 200
)

c.init <- list(
  zb = rep(0.15, M),
  zac = rep(0.55, M),
  w0b = 0.25,
  w1b = 0.25,
  w2b = 0.25,
  w3b = 0.25,
  ab = 0.7,
  rhob = 300,
  w0ac = 0.25,
  w1ac = 0.25,
  w2ac = 0.25,
  w3ac = 0.25,
  aac = 0.7,
  rhoac = 300
)

d.init <- list(
  zb = rep(0.55, M),
  zac = rep(0.65, M),
  w0b = 4,
  w1b = 4,
  w2b = 4,
  w3b = 4,
  ab = 0.1,
  rhob = 350,
  w0ac = 4,
  w1ac = 4,
  w2ac = 4,
  w3ac = 4,
  aac = 0.1,
  rhoac = 350
)

# put initialization into a list
init <- list(a.init, b.init, c.init, d.init)

# define number of chains
nc <- 4

# sample from the posterior distribution
fit.model <-
  sampling(
    object = lgp.model,
    data = feed.data,
    chains = nc,
    cores = nc,
    iter = 6000,
    warmup = 3000,
    seed = 1234,
    control = list(adapt_delta = 0.98, max_treedepth = 13),
    sample_file = "STAN/Chain.txt",
    init = init,
    verbose = TRUE,
    show_messages = TRUE,
    include = FALSE,
    pars = c("zb",
             "zac",
             "gfb",
             "gfac",
             "Kb",
             "Lb",
             "Kac",
             "Lac")
  )

# get results and store them in a dataframe
results <- summary(fit.model)$summary
results <- as.data.frame(results)
write.csv(results, "STAN/Results.csv")

# save stan model
fit.model@stanmodel@dso <- new("cxxdso")
saveRDS(fit.model, file = "STAN/rstanmodel.rds")

# save post warm-up samples from all chains altogether
combined <-
  extract(
    fit.model,
    inc_warmup = FALSE,
    include = TRUE,
    permuted = TRUE
  )
frame <- do.call(cbind.data.frame, combined)

write.csv(frame, "STAN/Combined-chains.csv", row.names = FALSE)
