library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
#Organize SCC biology data
dat<- read.csv("data/biologydata_south.central_2023.csv")


n1 <- names(dat)[grepl('calcofi.',names(dat))]
n2 <- names(dat)[grepl('rreas.',names(dat))]
ids <- c(n1,n2,"ZALOPHUS.PUPCT","ZALOPHUS.PUPWT")
n3 <- names(dat)[grepl('SBRD.',names(dat))]
n4 <- names(dat)[grepl('ZALOPHUS.',names(dat))]

for(i in 1:ncol(dat)){
  if(names(dat)[i] %in% ids){
    dat[,i] <- log(dat[,i])
  }
}


for(i in 1:ncol(dat)){
  if(names(dat)[i] %in% ids){
  dat[,i] <- dat[,i]-mean(na.omit(dat[,i])/sd(na.omit(dat[,i])))
  }
}
#dat<-dat%>%select(c(year,n4,n3,n1))
#### Running SCC  model####
dat.scc<-dat%>%select(c(year,n4,n1))
remelt = melt(dat.scc,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(TRUE, FALSE),
                       var_index = c("survey"), num_trends = 1:3,
                       elpd_loo = NA, se_elpd_loo=NA)

varIndx = c(rep(1,3),rep(2,14))

for(i in 1:nrow(model_df)) {
    fit.mod = fit_dfa(y = Y,
                      num_trends = model_df$num_trends[i],
                      iter=n_iter,
                      varIndx = varIndx,
                      chains=n_chains, estimate_nu=model_df$est_nu[i],
                      estimate_trend_ma = model_df$estimate_trend_ma[i],
                      estimate_trend_ar = model_df$estimate_trend_ar[i],
                      estimate_process_sigma = model_df$estimate_process_sigma[i],
                      seed=123)
    # extract log_lik for new data
    pars = rstan::extract(fit.mod$model)
  loo = loo(fit.mod)
  model_df$elpd_loo[i] = loo$estimates[1,1]
  model_df$se_elpd_loo[i] = loo$estimates[1,2]
  
  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))
  
}  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))
  

# finally run the best model from the lfo-cv above -- largest is best
#model_df = dplyr::arrange(model_df,-elpd_loo)

#saveRDS(fit.mod, file = paste0("/biology_",model_df$num_trends[1],
#                               "_",model_df$var_index[1],"_",str,theta_str,phi_str,sigma_str,".rds"))
varIndx = rep(1,nspp)

fit.mod.ccc = fit_dfa(y = Y,
                  num_trends = 1,
                  iter=n_iter,
                  varIndx = varIndx,
                  chains=n_chains, estimate_nu=model_df$est_nu[1],
                  estimate_trend_ma = model_df$estimate_trend_ma[1],
                  estimate_trend_ar = model_df$estimate_trend_ar[1],
                  estimate_process_sigma = model_df$estimate_process_sigma[1],
                  seed=123)

pars = rstan::extract(fit.mod.scc$model)
r.scc <- rotate_trends(fit.mod.scc)
p.scc <- plot_trends(r.scc)
p.scc
l.scc <- plot_loadings(r.scc,names=names)
l.scc
is_converged(fit.mod.scc)
summary(fit.mod.scc)

#### Running CCC  model####
dat.ccc<-dat%>%select(c(year,n3,n2))%>%
  filter(year>=1970)
remelt = melt(dat.ccc,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

n_chains = 3
n_iter = 8000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

varIndx = c(rep(1,8),rep(2,14))
varIndx = rep(1,nspp)

fit.mod.ccc = fit_dfa(y = Y,
                      num_trends = 1,
                      iter=n_iter,
                      varIndx = varIndx,
                      chains=n_chains, estimate_nu=model_df$est_nu[1],
                      estimate_trend_ma = model_df$estimate_trend_ma[1],
                      estimate_trend_ar = model_df$estimate_trend_ar[1],
                      estimate_process_sigma = model_df$estimate_process_sigma[1],
                      seed=123)

pars = rstan::extract(fit.mod.ccc$model)
r.ccc <- rotate_trends(fit.mod.ccc)
p.ccc <- plot_trends(r.ccc)
p.ccc
l.ccc <- plot_loadings(r.ccc,names=names)
l.ccc
is_converged(fit.mod.ccc)
summary(fit.mod.ccc)

