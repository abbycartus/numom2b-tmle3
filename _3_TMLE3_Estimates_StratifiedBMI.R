########################################################################################################################################
##################################### STRATIFIED TMLE                         ##########################################################
##################################### AUTHOR: ARC                             ##########################################################
##################################### DATE: 07.09.2020                        ##########################################################
#########################################################################################################################################
#--------------------------------------------------------------------------------------------------------------------
#PREPARATION --------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#INSTALL TLVERSE
#install.packages("devtools")
#devtools::install_github("tlverse/tlverse")

#INSTALL AND LOAD PACKAGES
packages <- c("foreach","doParallel","boot","rmutil","mvtnorm","gam","sandwich","ggplot2", "SuperLearner","xgboost",
              "devtools","glmnet","tmle","data.table","rpart","ranger","nnet","arm","earth","e1071","tidyverse",
              "npcausal", "sl3", "tlverse", "tmle3", "beepr", "fastDummies")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}
for (package in packages) {
  library(package, character.only=T)
}


# LOAD IN THE DATA HERE
setwd("filepath")
analysis <- read.csv("data", header=TRUE, sep=",")
analysis <- select(analysis, -X_merge)

# Dummy-code age (we don't have to worry about BMI here)
analysis <- dummy_cols(analysis, select_columns="agecat3")


# CREATE SUPERLEARNER LIBRARY
# Have to update this for sl3 or else we get an error in tmle3 
ranger_1 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=500, mtry=2, replace=T)
ranger_2 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=500, mtry=3, replace=T)
ranger_3 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=500, mtry=4, replace=T)
ranger_4 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=500, mtry=2, replace=F)
ranger_5 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=500, mtry=3, replace=F)
ranger_6 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=500, mtry=4, replace=F)
ranger_7 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=2500, mtry=2, replace=T)
ranger_8 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=2500, mtry=3, replace=T)
ranger_9 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=2500, mtry=4, replace=T)
ranger_10 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=2500, mtry=2, replace=F)
ranger_11 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=2500, mtry=3, replace=F)
ranger_12 <- make_learner(Lrnr_ranger, min.node.size=10, num.trees=2500, mtry=4, replace=F)

xgboost_1 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=4, shrinkage=0.01)
xgboost_2 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=4, shrinkage=0.001)
xgboost_3 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=4, shrinkage=0.0001)
xgboost_4 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=5, shrinkage=0.01)
xgboost_5 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=5, shrinkage=0.001)
xgboost_6 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=5, shrinkage=0.0001)
xgboost_7 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=6, shrinkage=0.01)
xgboost_8 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=6, shrinkage=0.001)
xgboost_9 <- make_learner(Lrnr_xgboost, ntrees=200, max_depth=6, shrinkage=0.0001)
xgboost_10 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=4, shrinkage=0.01)
xgboost_11 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=4, shrinkage=0.001)
xgboost_12 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=4, shrinkage=0.0001)
xgboost_13 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=5, shrinkage=0.01)
xgboost_14 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=5, shrinkage=0.001)
xgboost_15 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=5, shrinkage=0.0001)
xgboost_16 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=6, shrinkage=0.01)
xgboost_17 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=6, shrinkage=0.001)
xgboost_18 <- make_learner(Lrnr_xgboost, ntrees=500, max_depth=6, shrinkage=0.0001)
xgboost_19 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=4, shrinkage=0.01)
xgboost_20 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=4, shrinkage=0.001)
xgboost_21 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=4, shrinkage=0.0001)
xgboost_22 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=5, shrinkage=0.01)
xgboost_23 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=5, shrinkage=0.001)
xgboost_24 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=5, shrinkage=0.0001)
xgboost_25 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=6, shrinkage=0.01)
xgboost_26 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=6, shrinkage=0.001)
xgboost_27 <- make_learner(Lrnr_xgboost, ntrees=1000, max_depth=6, shrinkage=0.0001)
glmnet_learner <- make_learner(Lrnr_glmnet, alpha = seq(0,1,.2))
mean_learner <- make_learner(Lrnr_mean)
glm_learner <- make_learner(Lrnr_glm)
#knn_learner <- make_learner(Lrnr_knn)
# It doesn't look like there's a KNN learner functionality in sl3. 

# CREATE A LIST WITH ALL LEARNERS
lrn_list <- list(ranger_1, ranger_2, ranger_3, ranger_4, ranger_5, ranger_6, ranger_7, ranger_8, ranger_9,
                 ranger_10, ranger_11, ranger_12, xgboost_1, xgboost_2, xgboost_3, xgboost_4, xgboost_5, 
                 xgboost_6, xgboost_7, xgboost_8, xgboost_9, xgboost_10, xgboost_11, xgboost_12, xgboost_13, 
                 xgboost_14, xgboost_15,xgboost_16, xgboost_17, xgboost_18, xgboost_19, xgboost_20, xgboost_21, 
                 xgboost_22, xgboost_23, xgboost_24, xgboost_25, xgboost_26, xgboost_27, 
                 glmnet_learner, mean_learner, glm_learner)
#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------
# TMLE3 (WITH SAMPLE SPLITTING) -------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
# DEFINE SL_Y AND SL_A 
# We only need one, because they're the same
sl_lib <- Lrnr_sl$new(learners = lrn_list, metalearner = make_learner(Lrnr_nnls))

# PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
# Here, we need to take bmicat out of the node list because we're restricting to different levels of bmi
ate_spec <- tmle_ATE(treatment_level = 1, control_level = 0)
learner_list <- list(A = sl_lib, Y = sl_lib)
nodes_fruit <- list(W = c("smokerpre","married","insurpub","gdm","pct_emptyc", "college","black",
                          "v_totdens", #This is the only thing that changes
                          "d_totdens","p_totdens","g_whldens","g_nwhldens","fa_dens","seaplant_dens","sodium_dens"), 
                    A = "f_totdens80", 
                    Y="pree_acog")

nodes_veg <- list(W = c("smokerpre","married","insurpub","gdm","pct_emptyc", "college","black",
                        "f_totdens", #This is the only thing that changes 
                        "d_totdens","p_totdens","g_whldens","g_nwhldens","fa_dens","seaplant_dens","sodium_dens"), 
                  A = "v_totdens80", 
                  Y="pree_acog")

#We don't have to process missing because there isn't any
any(is.na(analysis))

# BMI: UW+NW, OW, OB
bmidat1 <- analysis %>% filter(bmicat == 1 | bmicat== 2)
bmidat2 <- analysis %>% filter(bmicat == 3)
bmidat3 <- analysis %>% filter(bmicat == 4 | bmicat == 5 | bmicat == 6)

setwd("filepath_to_results")

#-------------------------------------------------------------------------------
# FRUIT
#-------------------------------------------------------------------------------
set.seed(123)
tmle_fruit_bmiUWNW <- tmle3(ate_spec, bmidat1, nodes_fruit, learner_list)
saveRDS(tmle_fruit_bmiUWNW, "tmle_fruit_bmiUWNW.rds")
beep("mario")

set.seed(123)
tmle_fruit_bmiOW <- tmle3(ate_spec, bmidat2, nodes_fruit, learner_list)
saveRDS(tmle_fruit_bmiOW, "tmle_fruit_bmiOW.rds")
beep("mario")

set.seed(123)
tmle_fruit_bmiOB<- tmle3(ate_spec, bmidat3, nodes_fruit, learner_list)
saveRDS(tmle_fruit_bmiOB, "tmle_fruit_bmiOB.rds")
beep("mario")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# VEGETABLES
#-------------------------------------------------------------------------------
set.seed(123)
tmle_veg_bmiUWNW <- tmle3(ate_spec, bmidat1, nodes_veg, learner_list)
saveRDS(tmle_veg_bmiUWNW, "tmle_veg_bmiUWNW.rds")
beep("mario")

set.seed(123)
tmle_veg_bmiOW <- tmle3(ate_spec, bmidat2, nodes_veg, learner_list)
saveRDS(tmle_veg_bmiOW, "tmle_veg_bmiOW.rds")
beep("mario")

set.seed(123)
tmle_veg_bmiOB <- tmle3(ate_spec, bmidat3, nodes_veg, learner_list)
saveRDS(tmle_veg_bmiOB, "tmle_veg_bmiOB.rds")
beep("mario")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------