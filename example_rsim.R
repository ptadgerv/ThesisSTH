exp_basehaz <- function(t, lambda = 0.5) lambda * 1 * t^0
exp_weibull <- function(t, lambda = 0.5, gamma = 1.5) lambda * gamma * t^(gamma - 1)
curve(exp_basehaz, from = 0, to = 5, lty = 1, ylim = c(0, 2), ylab = expression(h[0](t)), xlab = "Follow-up time t")
curve(exp_weibull, from = 0, to = 5, lty = 2, add = TRUE)
legend(x = "topleft", lty = 1:2, legend = c("Exponential baseline hazard", "Weibull baseline hazard"), bty = "n")

simulate_data <- function(dataset, n, baseline, params = list(), coveff = -0.50) {
  # Simulate treatment indicator variable
  x <- rbinom(n = n, size = 1, prob = 0.5)
  # Draw from a U(0,1) random variable
  u <- runif(n)
  # Simulate survival times depending on the baseline hazard
  if (baseline == "Exponential") {
    t <- -log(u) / (params$lambda * exp(x * coveff))
  } else {
    t <- (-log(u) / (params$lambda * exp(x * coveff)))^(1 / params$gamma)
  }
  # Winsorising tiny values for t (smaller than one day on a yearly-scale, e.g. 1 / 365.242), and adding a tiny amount of white noise not to have too many concurrent values
  t <- ifelse(t < 1 / 365.242, 1 / 365.242, t)
  t[t == 1 / 365.242] <- t[t == 1 / 365.242] + rnorm(length(t[t == 1 / 365.242]), mean = 0, sd = 1e-4)
  # ...and make sure that the resulting value is positive
  t <- abs(t)

  # Make event indicator variable applying administrative censoring at t = 5
  d <- as.numeric(t < 5)
  t <- pmin(t, 5)
  # Return a data.frame
  data.frame(dataset = dataset, x = x, t = t, d = d, n = n, baseline = baseline, stringsAsFactors = FALSE)
}

set.seed(755353002)
reps <- 1:100
data <- list()
data[["n = 50, baseline = Exp"]] <- lapply(
  X = reps,
  FUN = simulate_data,
  n = 50,
  baseline = "Exponential",
  params = list(lambda = 0.5)
)
data[["n = 250, baseline = Exp"]] <- lapply(
  X = reps,
  FUN = simulate_data,
  n = 250,
  baseline = "Exponential",
  params = list(lambda = 0.5)
)
data[["n = 50, baseline = Wei"]] <- lapply(
  X = reps,
  FUN = simulate_data,
  n = 50,
  baseline = "Weibull",
  params = list(lambda = 0.5, gamma = 1.5)
)
data[["n = 250, baseline = Wei"]] <- lapply(
  X = reps,
  FUN = simulate_data,
  n = 250,
  baseline = "Weibull",
  params = list(lambda = 0.5, gamma = 1.5)
)


library(survival)
library(rstpm2)
library(eha)

fit_models <- function(data, model) {
  # Fit model
  if (model == "Cox") {
    fit <- survival::coxph(Surv(t, d) ~ x, data = data)
  } else if (model == "RP(2)") {
    fit <- rstpm2::stpm2(Surv(t, d) ~ x, data = data, df = 2)
  } else {
    fit <- eha::phreg(Surv(t, d) ~ x, data = data, dist = "weibull", shape = 1)
  }
  # Return relevant coefficients
  data.frame(
    dataset = unique(data$dataset),
    n = unique(data$n),
    baseline = unique(data$baseline),
    theta = coef(fit)["x"],
    se = sqrt(ifelse(model == "Exp", fit$var["x", "x"], vcov(fit)["x", "x"])),
    model = model,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


results <- list()
results[["n = 50, baseline = Exp, model = Cox"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 50, baseline = Exp"]],
    FUN = fit_models,
    model = "Cox"
  )
)
results[["n = 250, baseline = Exp, model = Cox"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 250, baseline = Exp"]],
    FUN = fit_models,
    model = "Cox"
  )
)
results[["n = 50, baseline = Wei, model = Cox"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 50, baseline = Wei"]],
    FUN = fit_models,
    model = "Cox"
  )
)
results[["n = 250, baseline = Wei, model = Cox"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 250, baseline = Wei"]],
    FUN = fit_models,
    model = "Cox"
  )
)

results[["n = 50, baseline = Exp, model = Exp"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 50, baseline = Exp"]],
    FUN = fit_models,
    model = "Exp"
  )
)
results[["n = 250, baseline = Exp, model = Exp"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 250, baseline = Exp"]],
    FUN = fit_models,
    model = "Exp"
  )
)
results[["n = 50, baseline = Wei, model = Exp"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 50, baseline = Wei"]],
    FUN = fit_models,
    model = "Exp"
  )
)
results[["n = 250, baseline = Wei, model = Exp"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 250, baseline = Wei"]],
    FUN = fit_models,
    model = "Exp"
  )
)

results[["n = 50, baseline = Exp, model = RP(2)"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 50, baseline = Exp"]],
    FUN = fit_models,
    model = "RP(2)"
  )
)
results[["n = 250, baseline = Exp, model = RP(2)"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 250, baseline = Exp"]],
    FUN = fit_models,
    model = "RP(2)"
  )
)
results[["n = 50, baseline = Wei, model = RP(2)"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 50, baseline = Wei"]],
    FUN = fit_models,
    model = "RP(2)"
  )
)
results[["n = 250, baseline = Wei, model = RP(2)"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data[["n = 250, baseline = Wei"]],
    FUN = fit_models,
    model = "RP(2)"
  )
)

class(results[["n = 50, baseline = Exp, model = Cox"]])
verlo<-results[["n = 50, baseline = Exp, model = Cox"]]


relhaz <- do.call(
  rbind.data.frame,
  results
)
view(relhaz)
row.names(relhaz) <- NULL


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
library(usethis)
usethis::use_data(relhaz, overwrite = TRUE)
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################


library(rsimsum)
s <- rsimsum::simsum(data = relhaz, estvarname = "theta", se = "se", true = -0.50, methodvar = "model", ref = "Cox", by = c("n", "baseline"))
s
summary(s)
