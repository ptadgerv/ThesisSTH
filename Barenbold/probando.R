# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))


library(cmdstanr)
library(rstan)


## Not run:

library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

# Set path to CmdStan
# (Note: if you installed CmdStan via install_cmdstan() with default settings
# then setting the path is unnecessary but the default below should still work.
# Otherwise use the `path` argument to specify the location of your
# CmdStan installation.)
set_cmdstan_path(path =  "C:\\Users\\32498\\Downloads\\cmdstan-2.31.0\\cmdstan-2.31.0")

# Create a CmdStanModel object from a Stan program,
# here using the example model that comes with CmdStan
setwd("C:/Users/32498/Documents/Master EPI/SEM6/Thesis/Barenbold")
file <- file.path( "Stan code of statistical model fit to observed prevalence simulation.stan")
mod <- cmdstan_model(file)
mod$print()

# Data as a named list (like RStan)
stan_data <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

# Run MCMC using the 'sample' method
fit_mcmc <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 2,
  parallel_chains = 2
)

# Use 'posterior' package for summaries
fit_mcmc$summary()

# Get posterior draws
draws <- fit_mcmc$draws()
print(draws)

# Convert to data frame using posterior::as_draws_df
as_draws_df(draws)

# Plot posterior using bayesplot (ggplot2)
mcmc_hist(fit_mcmc$draws("theta"))

# Call CmdStan's diagnose and stansummary utilities
fit_mcmc$cmdstan_diagnose()
fit_mcmc$cmdstan_summary()

# For models fit using MCMC, if you like working with RStan's stanfit objects
# then you can create one with rstan::read_stan_csv()

# stanfit <- rstan::read_stan_csv(fit_mcmc$output_files())


# Run 'optimize' method to get a point estimate (default is Stan's LBFGS algorithm)
# and also demonstrate specifying data as a path to a file instead of a list
my_data_file <- file.path(cmdstan_path(), "examples/bernoulli/bernoulli.data.json")
fit_optim <- mod$optimize(data = my_data_file, seed = 123)

fit_optim$summary()


# Run 'variational' method to approximate the posterior (default is meanfield ADVI)
fit_vb <- mod$variational(data = stan_data, seed = 123)

fit_vb$summary()

# Plot approximate posterior using bayesplot
mcmc_hist(fit_vb$draws("theta"))


# Specifying initial values as a function
fit_mcmc_w_init_fun <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 2,
  refresh = 0,
  init = function() list(theta = runif(1))
)
fit_mcmc_w_init_fun_2 <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 2,
  refresh = 0,
  init = function(chain_id) {
    # silly but demonstrates optional use of chain_id
    list(theta = 1 / (chain_id + 1))
  }
)
fit_mcmc_w_init_fun_2$init()

# Specifying initial values as a list of lists
fit_mcmc_w_init_list <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 2,
  refresh = 0,
  init = list(
    list(theta = 0.75), # chain 1
    list(theta = 0.25)  # chain 2
  )
)
fit_optim_w_init_list <- mod$optimize(
  data = stan_data,
  seed = 123,
  init = list(
    list(theta = 0.75)
  )
)
fit_optim_w_init_list$init()

## End(Not run)