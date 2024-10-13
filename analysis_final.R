#######STEP 0: Package loading#######
rm(list = ls(all.names = TRUE))
list.of.packages <- c( "grf", "metafor", "splitstackshape", "dplyr", "tidyverse", "foreach", "cowplot",  "bcf", "tidybayes", 
                       "doParallel", "survival", "readstata13", "ggplot2", "rsample", "DiagrammeR", "e1071", "BART", "tidytreatment", 
                       "pROC", "caret", "ModelMetrics", "MatchIt", "Hmisc", "scales", "lmtest", "sandwich", "cobalt", "tableone")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only = TRUE)
library(lmtest)
library(sandwich)

#######STEP 1: DATA Load & CLeaning#######
d<-read.csv("ALL_causal_forest_20230404.csv")
summary(d)
d$treat=(1-d$ric_mac5) #treat=1, intensified; treat=0, conventional
d$outcome=d$survival_1y #1, death; 0, survive
d$all_phenotype_b<-0
d$all_phenotype_b[d$all_phenotype==0]<-1
d$all_phenotype_t<-0
d$all_phenotype_t[d$all_phenotype==1]<-1
d$sct_type2_1<-0
d$sct_type2_1[d$sct_type2==1]<-1
d$sct_type2_2<-0
d$sct_type2_2[d$sct_type2==2]<-1
d$sct_type2_3<-0
d$sct_type2_3[d$sct_type2==3]<-1 ##unrelated PB: new strategy 
d$sct_type2_4<-0
d$sct_type2_4[d$sct_type2==4]<-1
d$sex_mismatch2_1<-0
d$sex_mismatch2_1[d$sex_mismatch2==1]<-1
d$sex_mismatch2_2<-0
d$sex_mismatch2_2[d$sex_mismatch2==2]<-1
d$abo_mismatch[d$abo_mismatch>=1 & d$abo_mismatch<=3]<-1
d$ks_tx<-as.numeric(d$ks_tx)

d.gf0<- d[, which (colnames(d) %in% 
                    c("treat", "outcome", 
                      "tyear_2011", 
                      "dx_to_sct",  
                      "pt_age", 
                      "pt_sex", 
                      "ps_tx", 
                      "pt_cmvab", 
                      "bac_fung_inf_tx", 
                      "all_phenotype_b",
                      "all_phenotype_t", 
                      "chromosome_ph", 
                      "rdri", 
                      "sct_type2_1",  
                      "sct_type2_2", 
                      "sct_type2_3", 
                      "sct_type2_4", 
                      "donor_sex", 
                      "sex_mismatch2_1",
                      "sex_mismatch2_2",
                      "abo_mismatch", 
                      "hla_mismatch",
                      "gvhd_pro2" 
                    ))]

#######STEP 2: Random forest imputation for missing values#######
library(missRanger)
d.gf.imp <- missRanger(d.gf0, num.trees = 100, pmm.k=3, verbose = 0, seed = 123)
datapre<-d.gf.imp

#######STEP 3: Propensity Score matching#######
set.seed(123)
fit.psm <- matchit(treat ~ .-outcome, method = "nearest", caliper =0.1, 
                   distance = "linear.logit",  data=datapre)
data.psm <- match.data(fit.psm)

summary(fit.psm)
bal.tab(fit.psm, thresholds = c(m = 0.1), un = TRUE)
bal.plot(fit.psm, var.name = "pt_age")
bal.plot(fit.psm, var.name = "pt_age", mirror = TRUE, type = "histogram")
love.plot(fit.psm, stats = c("mean.diffs"), 
          thresholds = c(m = 0.1), abs = TRUE, binary = "std", var.order = "unadjusted")

#######STEP 4: Bayesian Causal Forest#######
data<- data.psm[, which (colnames(data.psm) %in% 
                           c("treat", "outcome", 
                             "tyear_2011", 
                             "dx_to_sct", 
                             "pt_age", 
                             "pt_sex", 
                             "ps_tx", 
                             "pt_cmvab", 
                             "bac_fung_inf_tx", 
                             "all_phenotype_b",
                             "all_phenotype_t", 
                             "chromosome_ph", 
                             "rdri", 
                             "sct_type2_1", 
                             "sct_type2_2", 
                             "sct_type2_3", 
                             "sct_type2_4", 
                             "donor_sex", 
                             "sex_mismatch2_1", 
                             "sex_mismatch2_2",
                             "abo_mismatch", 
                             "hla_mismatch",
                             "gvhd_pro2" 
                           ))]
data<-as.data.frame(data)
summary(lm(outcome~treat+., data=data))  

  dat1<- data %>%
    select(treat, outcome, dx_to_sct, rdri, 
           pt_age, ps_tx, tyear_2011, 
           pt_sex, 
           pt_cmvab, 
           bac_fung_inf_tx, 
           all_phenotype_b,
           all_phenotype_t, 
           chromosome_ph, 
           sct_type2_1, 
           sct_type2_2, 
           sct_type2_3, 
           sct_type2_4, 
           donor_sex, 
           sex_mismatch2_1, 
           sex_mismatch2_2,
           abo_mismatch, 
           hla_mismatch,
           gvhd_pro2)
  Y <-  (dat1$outcome)
  W <-  (dat1$treat)
  X0 <-  subset(dat1, select=-Y)
  X1 <-  subset(X0, select=-T)
  X2 <- as.matrix(X1)
  model <- glm(treat ~ tyear_2011 + dx_to_sct + pt_age + pt_sex + ps_tx + pt_cmvab + 
                 bac_fung_inf_tx + all_phenotype_b + all_phenotype_t + chromosome_ph + 
                 rdri + sct_type2_1 + sct_type2_2 + sct_type2_3 + sct_type2_4 + donor_sex + 
                 sex_mismatch2_1 + sex_mismatch2_2 + abo_mismatch + hla_mismatch + gvhd_pro2, 
               family = binomial, data=dat1)
  pi <- predict(model, type = "response")
  
  model_bcf <- bcf(y             = Y,
                   z                = W,
                   x_control        = X2,
                   x_moderate       = X2,
                   pihat            = pi,
                   nburn            = 2500, 
                   nsim             = 2500,  
                   save_tree_directory = 'log',
                   log_file = file.path("log", sprintf("bcf_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S"))),
                   random_seed = 12345,
                   n_threads = 12)
  
summary(model_bcf)
tau.hat<-colMeans(model_bcf$tau)*(-1)

#######STEP 5: Calibration plot#######
  num.rankings<-4
  ranking <- rep(NA, nrow(dat1))
  tau.hat.quantiles <- quantile(tau.hat, probs = seq(0, 1, by=1/num.rankings))
  ranking <- cut(tau.hat, tau.hat.quantiles, include.lowest=TRUE,labels=seq(num.rankings))
  treatment <- "W"
  outcome <- "Y"
  fmla <- paste0(outcome, " ~ 0 + ranking + ranking:", treatment, "+tyear_2011+dx_to_sct+pt_age+pt_sex+ps_tx+pt_cmvab+bac_fung_inf_tx+all_phenotype_b+all_phenotype_t+chromosome_ph+
rdri+sct_type2_1+sct_type2_2+sct_type2_3+sct_type2_4+donor_sex+sex_mismatch2_1+sex_mismatch2_2+abo_mismatch+hla_mismatch+gvhd_pro2")
  ols.ate <- lm(fmla, data=transform(data, ranking=factor(ranking)))
  ols.ate <- coeftest(ols.ate, vcov=vcovHC(ols.ate, type='HC2'))
  interact <- which(grepl(":", rownames(ols.ate)))
  ols.ate <- data.frame("ols", paste0("Q", seq(num.rankings)), ols.ate[interact, 1:2])
  rownames(ols.ate) <- NULL # just for display
  colnames(ols.ate) <- c("method", "ranking", "estimate", "std.err")
  ols.ate
  res<-ols.ate
  res$estimate<-res$estimate
  ggplot(res) +
    #aes(x = ranking, y = estimate, group=method, color=method) + 
    aes(x = ranking, y = estimate*(-1)) + 
    geom_point(position=position_dodge(0.5), size=3) +
    geom_errorbar(aes(ymin=estimate*(-1)+1.96*std.err*(-1), ymax=estimate*(-1)-1.96*std.err*(-1)), width=.5, position=position_dodge(0.5)) +
    ylab("") + xlab("") +
    ggtitle("Average CATE within each ranking (as defined by predicted CATE)") +
    theme_minimal() +
    theme(legend.position="bottom", legend.title = element_blank())
  
#######STEP 6: Demographic#######
  newd.table<- subset(data, select=-outcome)
  vapply(newd.table, function(x) mean(!is.na(x)), numeric(1))
  
  i<-c(1:ncol(newd.table))
  newd.table[ , i] <- apply(newd.table[ , i], 2,         
                            function(x) as.numeric(as.character(x)))
  
  newd.table<-na.omit(newd.table)
  
  myVars <- c("dx_to_sct",   "pt_age")
  newd.table$dx_to_sct<-as.numeric( newd.table$dx_to_sct)
  newd.table$pt_age<-as.numeric( newd.table$pt_age)
  table(newd.table$treat)
  catVars <- c( "tyear_2011",
                "pt_sex", 
                "ps_tx",
                "rdri",
    "pt_cmvab", 
    "bac_fung_inf_tx", 
    "all_phenotype_b",
    "all_phenotype_t", 
    "chromosome_ph", 
    "sct_type2_1", 
    "sct_type2_2", 
    "sct_type2_3", 
    "sct_type2_4", 
    "donor_sex", 
    "sex_mismatch2_1", 
    "sex_mismatch2_2",
    "abo_mismatch", 
    "hla_mismatch",
    "gvhd_pro2"
  )
  newd.table[catVars] <- lapply(newd.table[catVars], factor)
  Table1 <- CreateTableOne(data = newd.table, strata="treat") 
  print(Table1) %>%
    write.csv("Table1_R1.csv")
  
  newd.table$benefit<-0
  newd.table$benefit[ranking==3| ranking==4]<-1 #notice this code should be changed based on outcome (increase or decrease)
  table(newd.table$benefit)
  Table2 <- CreateTableOne(data = newd.table, strata="benefit") 
  print(Table2) %>%
    write.csv("Table2_R1.csv")

#######STEP 7: Boostrap#######
  dat_boot<-dat1
  dat_boot$cate<-tau.hat
  write.csv(dat_boot, "psm_cate_R1.csv")
  
  models <- list()
  hba<-data.frame(1)
  estimate<-data.frame(1)
  difference<-data.frame(1)
  
  boot_hba =function(d){
    num.rankings=2
    ranking <- rep(NA, nrow(d))
    cate<-d$cate
    cate.quantiles <- quantile(cate, probs = seq(0, 1, by=1/num.rankings))
    ranking <- cut(cate, cate.quantiles, include.lowest=TRUE,labels=seq(num.rankings))
    
    model1<-(glm(outcome~treat+., family=gaussian(), data=d))
    
    model2<-(glm(outcome~treat+., family=gaussian(), data=d[d$rdri>=2, ]))
    
    model3<-(glm(outcome~treat+., family=gaussian(), data=d[d$pt_age<median(d$pt_age), ]))
    
    model4<-(glm(outcome~treat+., family=gaussian(), data=d[ranking==2, ]))
    
    for (i in seq(1, 4, by=1)) {
      model <- get(paste0("model", i))
      coef_target <- summary(model)$coefficients
      coef_all <- summary(model1)$coefficients
      hba[[i]] <- paste0("model", i)
      estimate[[paste0("model", i)]] <- coef_target["treat", "Estimate"]*(-1)
      difference[[paste0("model", i)]] <- coef_target["treat", "Estimate"]*(-1)-coef_all["treat", "Estimate"]*(-1)
    }
    testimate<-data.frame(t(estimate))
    testimate<-testimate[-1, ]
    tdifference<-data.frame(t(difference))
    tdifference<-tdifference[-1, ]
    thba<-data.frame(t(hba))
    diff.hba<-as.data.frame(cbind(thba, testimate,tdifference))
    diff.hba$t.hba. <- factor(diff.hba$t.hba., levels = paste0("model", i), ordered = TRUE)
    flattened_diff.hba <- as.vector(unlist(diff.hba))
    results<-flattened_diff.hba[5:12]
    return(results)
  }
  
  boot = function(d)
  { 
    #set.seed(Sys.time())
    return(d[sample(1:dim(d)[1], replace=TRUE), ])
  }
  
  set.seed(777)
  RD = replicate(1000, boot_hba(boot(dat_boot))) 
  results_est<-boot_hba(dat_boot)
  results_ci<-apply(RD, 1, quantile, c(0.025,  0.975))
  results<-rbind(results_est, results_ci)
  print(results) %>%
    write.csv(paste0("results_R1.csv"))