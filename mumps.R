## Stats 531, W21
##
## Author:  Jessica Leviton, Hongfan Chen
## Updated: April, 19, 2021
# 79: -------------------------------------------------------------------------
# Packages: 
library(tidyverse)
library(pomp)
library(doParallel)
library(doRNG)
# directories: ----------------------------------------------------------------
path = './'
## ----------------------------------------------------------------------------
## -------------------------------- DATA Import -------------------------------
## ----------------------------------------------------------------------------
## Project Tycho: Contagious Diseases, mumps data
mumps_file = sprintf('%s/mumps.csv', path)
mumps = read_delim(mumps_file, delim = ',' ) %>%
  select(week, state_name, cases)

mumps_data = mumps %>%
  filter(state_name == "MICHIGAN") %>%
  filter(week >= 197137 & week <= 197334) %>%
  select(cases) %>%
  mutate(week = 1:100) %>%
  relocate(week, cases)

mumps_data_t0 <- 0

mumps_data %>%
  ggplot(aes(x = week, y = cases))+
  geom_line(col = "tomato3")
## ----------------------------------------------------------------------------
## -------------------------------- SEIR Model: -------------------------------
## ----------------------------------------------------------------------------
seir_step <- Csnippet("
double Beta;
Beta = exp(b1 + b2 * cos(M_2PI/52*t - Phi));
double dN_SE = rbinom(S, 1-exp(-Beta*I/N*dt));
double dN_EI = rbinom(E, 1-exp(-mu_EI*dt));
double dN_IR = rbinom(I, 1-exp(-mu_IR*dt));
S -= dN_SE;
E += dN_SE - dN_EI;
I += dN_EI - dN_IR;
H += dN_IR;
")

seir_init <- Csnippet("
S = nearbyint(eta*N);
E = 20;
I = 10;
H = 0;
")

dmeas <- Csnippet("
lik = dnbinom(cases, H, rho, give_log);
")

rmeas <- Csnippet("
cases = rnbinom(H, rho);
")

mumpSEIR = mumps_data %>%
  pomp(
    times = "week",
    t0 = 0,
    rprocess = euler(seir_step, delta.t=1/7),
    rinit = seir_init,
    rmeasure = rmeas,
    dmeasure = dmeas,
    accumvars = "H",
    partrans = parameter_trans(
      logit = c("rho", "eta")
    ),
    statenames = c("S", "E", "I", "H"),
    paramnames = c("b1", "b2", "Phi", "mu_EI",
                   "mu_IR", "eta", "rho", "N"),
    cdir = "./",
    cfile = "mumpSEIR"
  )
## First attempts: ------------------------------------------------------------
mumps_fixed_params = c(N = 8881826, mu_EI = 0.412, mu_IR = 0.714)

params = c(b1 = 1, b2 = 1, Phi = 0.1,
           rho = 0.8, eta = 0.0216,
           mumps_fixed_params)

y = mumpSEIR %>%
  simulate(params = params,
           nsim = 10,
           format = "data.frame",
           include.data = TRUE)

y %>%
  ggplot(aes(x = week, y = cases,
             group = .id, color = .id=="data")
  ) +
  geom_line() +
  labs(x = "Weeks",
       y = "Reporting Cases",
       color = "Original Data")
## ----------------------------------------------------------------------------
## ------------------------- Computational Level ------------------------------
## ----------------------------------------------------------------------------
run_level = 2
mumps_Np = switch(run_level, 100, 1e3, 2e3)
mumps_Nmif = switch(run_level, 10, 100, 150)
mumps_Nreps_eval = switch(run_level, 2, 10, 20)
mumps_Nreps_local = switch(run_level, 10, 30, 40)
mumps_Nreps_global = switch(run_level, 10, 60, 100)
mumps_Nsim = switch(run_level, 50, 70, 100)

## ----------------------------------------------------------------------------
## -------------------- Local Maximization of Likelihood ----------------------
## ----------------------------------------------------------------------------
cl = makeCluster(8)
registerDoParallel(cl)
registerDoRNG(2021531)
mifs_local = foreach(i = 1:mumps_Nreps_local,
                     .packages = c("pomp", "tidyverse"),
                     .combine = c) %dopar% { 
  mumpSEIR %>%
    mif2(
      params = params,
      Np = mumps_Np, 
      Nmif = mumps_Nmif,
      cooling.fraction.50 = 0.5,
      rw.sd = rw.sd(b1 = 0.02, b2 = 0.02, Phi = 0.02,
                    rho = 0.02, eta = ivp(0.02)
                    )
    )
}
stopCluster(cl)

cl = makeCluster(8)
registerDoParallel(cl)
registerDoRNG(2021531)
lik_local = foreach(i = 1:mumps_Nreps_local,
                     .packages = c("pomp", "tidyverse"),
                     .combine=rbind) %dopar% {
  logmeanexp(
    replicate(mumps_Nreps_eval,
              logLik(pfilter(mumpSEIR,
                             params =  coef(mifs_local[[i]]),
                             Np = mumps_Np)
                     )
              ),
    se = TRUE)
}
stopCluster(cl)
## Local: Parameter and likelihood --------------------------------------------
r_local = t(sapply(mifs_local, coef)) %>%
  as_tibble() %>%
  bind_cols(tibble(logLik = lik_local[,1],
                   logLik_se = lik_local[,2])
  ) %>%
  arrange(-logLik) %>%
  head(10)
## Local: Parameter movement --------------------------------------------------
mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x = iteration,
             y = value,
             group = L1,
             color = factor(L1)
             )
         )+
  geom_line()+
  guides(color = FALSE)+
  facet_wrap(~variable,
             scales = "free_y")
## Fit model using params selected from local search and plot it out: ---------
params_name = c("b1", "b2", "Phi", "rho", "eta")
local_params = r_local[1,]
best_local = unlist(c(local_params[params_name], mumps_fixed_params))

set.seed(643345567)
mod_local = mumpSEIR %>%
  simulate(params = best_local,
           nsim = 1,
           format = "data.frame",
           include.data = TRUE)


mod_local %>%
  ggplot(aes(x = week,
             y = cases,
             group = .id,
             color = factor(.id)
  )
  ) +
  geom_line() +
  scale_color_brewer(type = "qual",
                     palette = 6
  )+
  guides(color = FALSE)

## ----------------------------------------------------------------------------
## -------------------- Global Maximization of Likelihood ---------------------
## ----------------------------------------------------------------------------
mumps_box = rbind(
  b1 = c(0,5), b2 = c(0,5), Phi = c(0, 2*pi),
  eta = c(0,0.10), rho = c(0, 0.9)
)

cl = makeCluster(8)
registerDoParallel(cl)
registerDoRNG(2021531)
mifs_global = foreach(i = 1:mumps_Nreps_global,
                       .packages = 'pomp', 
                       .combine = c) %dopar%{
  mif2(mifs_local[[1]],
       params = c(apply(mumps_box,
                        1,
                        function(x) runif(1, x[1], x[2])),
                  mumps_fixed_params
                  )
       )
}
stopCluster(cl)

cl = makeCluster(8)
registerDoParallel(cl)
registerDoRNG(2021531)
lik_global = foreach(i = 1:mumps_Nreps_global,
                     .packages = 'pomp',
                     .combine = rbind) %dopar% {
  logmeanexp(
    replicate(mumps_Nreps_eval, 
              logLik(pfilter(mumpSEIR,
                             params = coef(mifs_global[[i]]),
                             Np = mumps_Np)
                     )
              ),
    se = TRUE
    )
}
stopCluster(cl)
## Global: Parameter and likelihood -------------------------------------------
r_global = t(sapply(mifs_global, coef)) %>%
  as_tibble() %>%
  bind_cols(tibble(logLik = lik_global[,1],
                   logLik_se = lik_global[,2])
  ) %>%
  arrange(-logLik) %>%
  head(10)

pairs( ~ logLik + b1 + b2 + Phi+ rho + eta ,
       data = r_global, pch = 16)
## Global: Parameter movement -------------------------------------------------
mifs_global %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x = iteration,
             y = value,
             group = L1,
             color = factor(L1)
  )
  )+
  geom_line()+
  guides(color = FALSE)+
  facet_wrap(~variable,
             scales = "free_y")

## Fit model using params selected from global search and plot it out: --------
global_params = r_global[1,]
best_global = unlist(c(global_params[params_name], mumps_fixed_params))

set.seed(238765234)
mod_global = mumpSEIR %>%
  simulate(params = best_global,
           nsim = 1,
           format = "data.frame",
           include.data = TRUE)

mod_global %>%
  ggplot(aes(x = week,
             y = cases,
             group = .id,
             color = .id=="data")
  ) +
  geom_line() +
  labs(x = "Weeks",
       y = "Reporting Cases",
       color = "Original Data")
