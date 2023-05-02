climate <- readRDS('Output/standat.rds')

N <- length(climate$period)
NP <- 9
P<- as.numeric(as.factor(climate$era.region))
K <- 1
x <- data.frame(climate$stand_Bakun)
y<-climate$stand_PDO
data <- list(N = N, #total number of opservations
             NP=NP, #total number of time periods
             P = P, #time period pointer vector
             K = K,#number of covariates, starting with one but can add for final model structure
             x=x, #Upwelling
             y=y #response variable
)

warmups <- 1000
total_iterations <- 3000
max_treedepth <-  10
n_chains <-  3
n_cores <- 4
adapt_delta <- 0.95

bhfit <- stan(
  file = here::here("AtmosClimateProj/Src/BayesianLinearHierarchicalModels.stan"),
  data = data,
  chains = n_chains,
  warmup = warmups,
  iter = total_iterations,
  cores = n_cores,
  refresh = 250,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta))

posteriors<-data.frame(summary(bhfit, prob=c(0.025, 0.25,0.75, 0.975, 0.1, 0.9))$summary)
saveRDS(posteriors, file = "Output/posteriors.rds")

