########################################################################### 
######################### PLOTTING FOR SER POSTER ######################### 
########################################################################### 


#-------------------------------------------------------------------------------------------------------------------------------------
#PREPARATION
#-------------------------------------------------------------------------------------------------------------------------------------
# INSTALL AND LOAD PACKAGES
#INSTALL AND LOAD PACKAGES
packages <- c("foreach","doParallel","boot","rmutil","mvtnorm","gam","sandwich","ggplot2", "SuperLearner","xgboost",
              "devtools","glmnet","tmle","data.table","rpart","ranger","nnet","arm","earth","e1071","tidyverse",
              "sl3", "tlverse", "tmle3", "beepr")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}
for (package in packages) {
  library(package, character.only=T)
}


######### Read in the TMLE object we created
pree_fruit <- readRDS(file="/Users/abigailcartus/Box Sync/Temporary TMLE folder/tmle3_fruit_pree.rds")
dat <- read.csv(file="/Users/abigailcartus/Box Sync/Temporary TMLE folder/numom_small_2020.08.21.csv", header=TRUE, sep=",")
# Y = outcome (pree_acog), A = exposure (f_totdens80), W = covariates
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------
# PREPARE FOR BOOTSTRAPPING
#-------------------------------------------------------------------------------------------------------------------------------------
# Get stuff out of the tmle3 object:

# Propensity scores
g.cf_task1 <- pree_fruit$tmle_task$generate_counterfactual_task("cf1",data.frame(A=1))
pihat_0 <- pree_fruit$likelihood$get_likelihood(g.cf_task1,"A")

# Outcome predictions
Q.cf_task0 <- pree_fruit$tmle_params[[1]]$cf_likelihood_control$cf_tasks[[1]]
Q.cf_task1 <- pree_fruit$tmle_params[[1]]$cf_likelihood_treatment$cf_tasks[[1]]
mu0 <- as.numeric(pree_fruit$likelihood$get_likelihood(Q.cf_task0,"Y", "validation"))
mu1 <- as.numeric(pree_fruit$likelihood$get_likelihood(Q.cf_task1,"Y", "validation"))
muhat_0 <-  mu1*(dat$f_totdens80) + mu0*(1-dat$f_totdens80) 

# Take the AIPW EFF 
X_ <- dat$f_totdens80
Y_ <- dat$pree_acog
aipw_EFF <- as.numeric((((2*X_-1)*(Y_ - muhat_0))/((2*X_-1)*pihat_0 + (1-X_)) + mu1 - mu0))

# Add AIPW to data so we're resampling it too
# I know I'm reading the data in a second time -- may change later
analysis <- read.csv(file="/Users/abigailcartus/Box Sync/Temporary TMLE folder/numom_small_2020.08.21.csv", 
                     sep=",", header=TRUE)
analysis <- cbind(analysis, aipw_EFF)

# We have to get rid of missing in BMI or we get an error. I am doing mean imputation:
mean(analysis$bmi, na.rm=TRUE)
analysis <- analysis %>% replace_na(list(bmi = 26.09507))


#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# BOOTSTRAPPING ------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
bootNum = 5
res <- NULL
for(jj in 1:bootNum){
  set.seed(jj)
  boot_dat <- analysis %>% sample_n(.,nrow(analysis), replace=T)

  # have to add as.numeric to get rid of error 
  M_ = as.numeric(boot_dat$bmi)
  y = boot_dat$aipw_EFF
  #M0_ = cbind(0, M_)
  x <- M_
  D <- data.frame(y,x)
  
  #set.seed(123)
  mm_numom <- seq(min(M_), max(M_), by = 1)
  print(jj)
  
  folds=10
  index<-split(1:nrow(D),1:folds)
  splt <- lapply(1:folds,function(ind) D[index[[ind]],])
  SL.library <- c("SL.glm","SL.step", "SL.earth", "SL.gam", "SL.mean", "SL.bayesglm")
  
  # hand coded algorithms
  m1 <- lapply(1:folds,function(ii) gam(y~s(x,5),family="gaussian",data=rbindlist(splt[-ii])))
  m2 <- lapply(1:folds,function(ii) gam(y~s(x,4),family="gaussian",data=rbindlist(splt[-ii])))
  m3 <- lapply(1:folds,function(ii) gam(y~s(x,3),family="gaussian",data=rbindlist(splt[-ii])))
  m4 <- lapply(1:folds,function(ii) mean(rbindlist(splt[-ii])$y))
  m5 <- lapply(1:folds,function(ii) bayesglm(y~x,data=rbindlist(splt[-ii]),family=gaussian))
  
  
  p1 <- lapply(1:folds,function(ii) predict(m1[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p2 <- lapply(1:folds,function(ii) predict(m2[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p3 <- lapply(1:folds,function(ii) predict(m3[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p4 <- lapply(1:folds,function(ii) m4[[ii]])
  p5 <- lapply(1:folds,function(ii) predict(m5[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  
  for(i in 1:folds){
    splt[[i]] <- cbind(splt[[i]],p1[[i]],p2[[i]],p3[[i]],p4[[i]],p5[[i]])
  }
  
  X <- data.frame(do.call(rbind,splt))[,-2]
  names(X) <- c("y","gam1","gam2","gam3","mean","bayesglm")
  head(X)
  
  SL.r <- nnls(cbind(X[,2],X[,3],X[,4],X[,5],X[,6]),X[,1])$x
  alpha <- as.matrix(SL.r/sum(SL.r))
  round(alpha,3)
  
  m1_new <- gam(y~s(x,5),family="gaussian",data=D)
  m2_new <- gam(y~s(x,4),family="gaussian",data=D)
  m3_new <- gam(y~s(x,3),family="gaussian",data=D)
  m4_new <- mean(D$y)
  m5_new <- bayesglm(y~x,data=D,family=gaussian)
  
  p1_new <- predict(m1_new,newdata=data.frame(x=mm_numom),type="response")
  p2_new <- predict(m2_new,newdata=data.frame(x=mm_numom),type="response")
  p3_new <- predict(m3_new,newdata=data.frame(x=mm_numom),type="response")
  p4_new <- m4_new
  p5_new <- predict(m5_new,newdata=data.frame(x=mm_numom),type="response")
  
  SL.predict <- cbind(p1_new,p2_new,p3_new,p4_new,p5_new)%*%as.matrix(round(alpha,3))
  
  #identical(SL.predict,as.matrix(p1_new))
  
  res <- rbind(res,cbind(jj,SL.predict,mm_numom))
}

head(res)
tail(res)


d <- tibble(boot = res[,1],psi = res[,2], m = res[,3])
head(d)

d <- d %>% mutate(m_round = round(m, 0))
#d_test <- d %>% filter(boot==1 | boot==3)
d <- na.omit(d)
#psi is missing (NA) for some of the bootstrap resamples 
#is it okay to exclude them? 

d_mean <- aggregate(d$psi, list(d$m_round), mean)
colnames(d_mean) <- c("bmi", "rd_mean")


d_sd <- aggregate(d$psi, list(d$m_round), sd)
colnames(d_sd) <- c("bmi","sd_mean") 
d_sd <- d_sd %>% select("sd_mean")

plot_dat <- cbind(d_mean, d_sd)
plot_dat <- plot_dat %>% mutate(sd_plus = rd_mean + sd_mean,
                                sd_minus = rd_mean - sd_mean)

#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# PLOTTING ------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
ggplot(plot_dat) + 
  geom_smooth(aes(y = rd_mean, x = bmi),color="black", alpha = 1, size = 0.9, se=F) +
  geom_ribbon(aes(ymin = sd_minus, ymax =sd_plus, x = bmi), linetype = 2, size = 1.2, color="gray20", alpha = 0.1) +
  labs(title = "Outcome: Preeclampsia",
       subtitle = "Exposure: Total fruit above 80th percentile",
       x = "BMI (kg/m^2)",
       y = "Risk difference") +
  geom_hline(yintercept = 0,color="red",size=.9) +
  theme_bw()+
  scale_x_continuous(limits = c(15,45), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.09, 0.09), expand = c(0, 0)) +
  theme(plot.title=element_text(size=75),
        plot.subtitle=element_text(size=40),
        axis.title.x=element_text(size=60),
        axis.title.y=element_text(size=60),
        axis.text.x=element_text(size=60),
        axis.text.y=element_text(size=60),
        plot.margin = unit(c(1,1,1,1), "cm"))
#p + geom_ribbon(aes(ymin = plot_dat$sd_minus, ymax = plot_dat$sd_plus, x = plot_dat$heix_tot), colour="gray23", linetype = 2, alpha = 0.1)

ggsave("/Users/abigailcartus/Box Sync/Temporary TMLE folder/Results/Plots/FruitPREE_BMI.png", width = 20, height = 20, dpi = 350)
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
