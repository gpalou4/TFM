rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(ggplot2)
cat("\n\n ggplot2 LOADED \n\n", file = stdout())
# library(MOFA2)
# cat("\n\n MOFA LOADED \n\n", file = stdout())
library(survival)
cat("\n\n survival LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: COX SURVIVAL SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

cat("\n\n... 1) LOADING DATA ...\n\n", file = stdout())

samples <- "914_2055"
groups <- "nogroups"

cat("\n\n... 1.1) LOADING SAMPLES METADATA ...\n\n", file = stdout())

MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv", sep = ""),
                            header=TRUE, stringsAsFactors = F, sep=",")
# Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
colnames(MOFA.covariates)[1] <- "sample"
# Change sample names from integers to characters
MOFA.covariates$sample <- as.character(MOFA.covariates$sample)
MOFA.covariates[1:3,]
dim(MOFA.covariates)
# 2055 24

samples <- "914_2055"
groups <- "ewas_nogroups"

# cat("\n\n... 1.2) LOADING MOFA TRAINED OBJECT ...\n\n", file = stdout())
# 
# MOFA.model <- load_model("MOFA/2.training_model/MOFA_trained_ewas_nogroups_914_2055.hdf5")
# MOFA.model
# 
# cat("\n\n... 2) ADDING SAMPLES METADATA (COVARIABLES) TO THE TRAINED MODEL OBJECT ...\n\n", file = stdout())
# 
# samples_metadata(MOFA.model) <- MOFA.covariates
# 
# cat("\n\n... 3) INSPECTING THE MOFA OBJECT ...\n\n", file = stdout())
# 
# cat("\nSlot names of the MOFA model\n", file = stdout())
# slotNames(MOFA.model)
# cat("\nMetadata stored in the MOFA model\n", file = stdout())
# head(MOFA.model@samples_metadata)

# Add another covariable: CHD incidence
chd <- read.table("/projects/regicor/data/FHS/phenotype/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", header=T, stringsAsFactors=F, sep="\t")
head(chd)
dim(chd)

chd.cov <- chd[chd$shareid%in%MOFA.covariates$sample,c("shareid","chd")]
colnames(chd.cov)[1] <- "sample"
MOFA.covariates.chd <- merge(MOFA.covariates, chd.cov, by = "sample")
MOFA.covariates.chd.filt <- MOFA.covariates.chd[MOFA.covariates.chd$cvd==1,]

# Order
MOFA.covariates.chd <- MOFA.covariates.chd[match(MOFA.covariates$sample,MOFA.covariates.chd$sample),]
MOFA.covariates <- MOFA.covariates.chd

cat("\n\n... 4) COX SURVIVAL ANALYSIS ...\n\n", file = stdout())

cat("\n\n... 4.1) Preparing the data: biomarker, cvd and risk factors ...\n\n", file = stdout())

##################################################################################################
####### Biomarker (Factor), cardiovascular outcome (CVD) and risk factors (covariables) #######
##################################################################################################

# factors <- get_factors(MOFA.model,
#                        groups = "all",
#                        factors = "all"
# )
# lapply(factors,dim)
# lapply(factors,head)
# 
# cat("\n\nFactor 14\n\n", file = stdout())
# factor14 <- factors$group1[,14]
# head(factor14)

factors.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/factors_values.RData", sep = ""))
factors.values <- get(factors.name)
factor <- 23
factors.values$group1[1:5,factor]
dim(factors.values$group1)

# factor <- 9
# factor.cov.dates <- read.csv(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/factor_cov/factor",factor,"_cov_dates.csv", sep = ""))
# head(factor.cov.dates)
# dim(factor.cov.dates)

## Stratification by sex

# Male
# factor.cov.dates <- factor.cov.dates[factor.cov.dates$sex==1,]
# sex <- "Male"
# Female
#factor.cov.dates <- factor.cov.dates[factor.cov.dates$sex==2,]
#sex <- "Female"

# features expression

# features.expression.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/features_expression_f12.RData", sep = ""))
# features.expression <- get(features.expression.name)
# features.expression[1:5,1:5]
# dim(features.expression)

# Add the interested factor values to the covariates dataframe

coxRegression <- function(biomarker) {
  
  # table(MOFA.covariates$smoke)
  # which(MOFA.covariates$smoke=="f>5")
  # MOFA.covariates[which(MOFA.covariates$smoke=="f>5"),"smoke"] = "f1-5"
  # MOFA.covariates[MOFA.covariates$smoke == "f1-5","smoke"] = "f>1"
  # table(MOFA.covariates$smoke)
  
  MOFA.covariates$factor <- biomarker
  head(MOFA.covariates)
  dim(MOFA.covariates)
  
  # Risk factors + cvd
  cat("\n\nFilter necessary variables: CVD, risk factors, factor and shareid\n\n", file = stdout())
  colnames(MOFA.covariates)[1] <- "shareid"
  filter.cov <- c("shareid","factor", "cvd", "chd", "age", "sex", "tot_chol", "hdl_chol", "sbp", "dbp", "glucose", "smoke", "CD8T","CD4T","NK","Bcell","Mono","Gran","sva1_METH")
  factor.cov <- MOFA.covariates[,filter.cov]
  head(factor.cov)
  dim(factor.cov)
  
  cat("\n\nObtain initial follow-up date\n\n", file = stdout())
  initial.followup <- read.table("/projects/regicor/data/FHS/phenotype/phs000007.v25.pht003099.v3.p9.c1.vr_dates_2011_m_0689s.HMB-IRB-MDS.txt", header=T, stringsAsFactors=F, sep="\t")
  initial.followup <- initial.followup[,c("shareid","date8")]
  head(initial.followup)
  dim(initial.followup)
  
  cat("\n\nObtain final follow-up date\n\n", file = stdout())
  end.followup <- read.table("/projects/regicor/data/FHS/phenotype/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", header=T, stringsAsFactors=F, sep="\t")
  end.followup <- end.followup[,c("shareid","cvddate","chddate")]
  head(end.followup)
  dim(end.followup)
  
  cat("\n\nCalculate CVD follow-up date\n\n", file = stdout())
  dates.followup <- merge(initial.followup, end.followup, by = "shareid")
  factor.cov.dates <- merge(dates.followup, factor.cov, by = "shareid")
  factor.cov.dates$cvdfollowup = factor.cov.dates$cvddate - factor.cov.dates$date8
  factor.cov.dates$chdfollowup = factor.cov.dates$chddate - factor.cov.dates$date8
  
  cat("\n\nFinal variables dataframe to calculate cox regression\n\n", file = stdout())
  head(factor.cov.dates)
  dim(factor.cov.dates)
  
  write.csv(factor.cov.dates, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/factor_cov/factor",factor,"_cov_dates.csv", sep = ""),
            row.names = FALSE)
  
  cat("\nSummary of follow-up time for CVD individuals\n", file = stdout())
  summary(factor.cov.dates[factor.cov.dates$cvd == 1, "cvdfollowup"])
  cat("\nSummary of follow-up time for non-CVD individuals\n", file = stdout())
  summary(factor.cov.dates[factor.cov.dates$cvd == 0, "cvdfollowup"])
  
  cat("\nSummary of follow-up time for CVD individuals\n", file = stdout())
  summary(factor.cov.dates[factor.cov.dates$cvd == 1, "chdfollowup"])
  cat("\nSummary of follow-up time for non-CVD individuals\n", file = stdout())
  summary(factor.cov.dates[factor.cov.dates$cvd == 0, "chdfollowup"])
  
  ####################################################################################
  ############################## Follow-up median ####################################
  ####################################################################################
  
  irq.cvd <- (IQR(factor.cov.dates$followup))/365.25
  median.cvd <- (median(factor.cov.dates$followup, na.rm=T))/365.25
  iqr.cvd.low <- median.cvd-irq.cvd
  iqr.cvd.up <- median.cvd+irq.cvd
  
  
  ##################################################################################################
  ### Association between biomarker (factor 14) and incidence of CVD (we need the time to event) ###
  ##################################### COX REGRESSION #############################################
  ##################################################################################################
  
  ## prueba
  
  glm.res <- glm(cvd ~ factor, data = factor.cov.dates, family = "binomial")
  summary(glm.res)
  # El t-test debería dar lo mismo que el hecho en statistical_tests!
  t.test(factor ~ cvd, data = factor.cov.dates)
  
  mod.surv.cvd.1 <- coxph(Surv(followup, cvd) ~
                            factor,
                          data=factor.cov.dates)
  cat("\n\nResults for model 1: factor11\n\n", file = stdout())
  print(mod.surv.cvd.1)
  
  # Add Confidence Intervals
  mod.surv.cvd.1.df <- as.data.frame(summary(mod.surv.cvd.1)$coefficients)
  print(mod.surv.cvd.1.df)
  mod.surv.cvd.1.df$lower95 <- as.numeric(exp(confint(mod.surv.cvd.1)))[1]
  mod.surv.cvd.1.df$upper95 <- as.numeric(exp(confint(mod.surv.cvd.1)))[2]
  names(mod.surv.cvd.1.df) <- c("coef","HR","se(coef)","z","P", "lower95", "upper95")
  print(mod.surv.cvd.1.df)
  
  mod.surv.cvd.2 <- coxph(Surv(followup, cvd) ~
                            factor+
                            sex+
                            age,
                          data=factor.cov.dates)
  cat("\n\nResults for model 2: factor11 + sex + age\n\n", file = stdout())
  print(mod.surv.cvd.2)
  
  # Add Confidence Intervals
  mod.surv.cvd.2.df <- as.data.frame(summary(mod.surv.cvd.2)$coefficients)
  print(mod.surv.cvd.2.df)
  mod.surv.cvd.2.df$lower95 <- as.numeric(exp(confint(mod.surv.cvd.2)))[1:3]
  mod.surv.cvd.2.df$upper95 <- as.numeric(exp(confint(mod.surv.cvd.2)))[4:6]
  names(mod.surv.cvd.2.df) <- c("coef","HR","se(coef)","z","P", "lower95", "upper95")
  print(mod.surv.cvd.2.df)
  
  mod.surv.cvd.3 <- coxph(Surv(followup, cvd) ~
                            factor+
                            sex+
                            age+
                            tot_chol+
                            hdl_chol+
                            sbp+
                            dbp+
                            glucose+
                            smoke,
                          data=factor.cov.dates)
  
  cat("\n\nResults for model 3: factor + sex + age + tot_chol + hdl_chol + sbp + dbp + smoke + glucose\n\n", file = stdout())
  print(mod.surv.cvd.3)
  
  # Add Confidence Intervals
  mod.surv.cvd.3.df <- as.data.frame(summary(mod.surv.cvd.3)$coefficients)
  print(mod.surv.cvd.3.df)
  print(as.numeric(exp(confint(mod.surv.cvd.3))))
  mod.surv.cvd.3.df$lower95 <- as.numeric(exp(confint(mod.surv.cvd.3)))[1:11]
  mod.surv.cvd.3.df$upper95 <- as.numeric(exp(confint(mod.surv.cvd.3)))[12:22]
  names(mod.surv.cvd.3.df) <- c("coef","HR","se(coef)","z","P", "lower95", "upper95")
  print(mod.surv.cvd.3.df)
  
  return(mod.surv.cvd.3)

}

for (i in 1:30) {

  print(paste("ESTE ES EL FACTOR ",i,sep = ""))
  factor <- i
  biomarker <- factors.values$group1[,factor]
  #feature <- i
  #biomarker <- features.expression[,feature]
  cox.res <- coxRegression(biomarker = biomarker)

}

stop("till here")

### kaplan meier (can only be done for categorical data...)
## tres perfiles: 1-33, 33-66, 66-100

# summary(factor.cov.dates$factor)
# factor.quantiles <- quantile(factor.cov.dates$factor, probs = seq(0, 1, 0.33))
# factor.categorical <- ifelse(factor.cov.dates$factor < factor.quantiles[2], 1, factor.cov.dates$factor)
# factor.categorical <- ifelse(factor.categorical < factor.quantiles[3], 2, factor.categorical)
# factor.categorical <- ifelse(factor.categorical == 1 | factor.categorical == 2, factor.categorical, 3)
# factor.cov.dates$factor_categorical <- factor.categorical
# 
# survObject <- Surv(factor.cov.dates$followup, factor.cov.dates$cvd)
# #fit <- survfit(survObject ~ 1, data = factor.cov.dates)
# fit <- survfit(survObject ~ factor_categorical, data = factor.cov.dates)
# plot(fit)
# 
# # significativo?
# logRank <- survdiff(survObject ~ factor_categorical, data = factor.cov.dates)
# pval <- p.val <- 1 - pchisq(logRank$chisq, length(logRank$n) - 1)
# pval
# 
# library(survminer)
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/kaplan_meier_",factor,".png", sep = ""), height=3600, width=6000, res=600, units="px")
# ggsurvplot(fit, data=factor.cov.dates, pval = pval, legend.labs = c("Factor values tertile 1", "Factor values tertile 2", "Factor values tertile 3" ))
# dev.off()

###############################################################################################
############################################ Discriminaci?n ###################################
###############################################################################################

library(Hmisc)
#install.packages("data.table")
#Estad?stico C 
# Harrell's C for survival data (implemented in the rcorr.cens function of the Hmisc package).
# rcorr.cens Rank Correlation for Censored Data Computes the c index and the corresponding generalization of Somers' Dxy rank correlation for a censored response variable. Also works for uncensored and binary responses, 
# although its use of all possible pairings makes it slow for this purpose. Dxy and c are related by Dxy = 2(c ??? 0.5). rcorr.cens handles one predictor variable. rcorrcens computes rank correlation measures separately by a series of predictors. 
# In addition, rcorrcens has a rough way of handling categorical predictors. If a categorical (factor) predictor has two levels, it is coverted to a numeric having values 1 and 2. If it has more than 2 levels, 
# an indicator variable is formed for the most frequently level vs.all others, and another indicator for the second most frequent level and all others. The correlation is taken as the maximum of the two (in absolute value).

discr <- function(factor, sex = "All") {
  
  factor <- factor
  factor.cov.dates <- read.csv(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/factor_cov/factor",factor,"_cov_dates.csv", sep = ""), sep = ",", stringsAsFactors = FALSE)
  head(factor.cov.dates)
  dim(factor.cov.dates)
  
  ## Stratification by sex (only for factor 21)
  
  if (sex == "Male") {
    factor.cov.dates <- factor.cov.dates[factor.cov.dates$sex==1,]
    sex <- "Male" 
  }
  else if (sex == "Female") {
    factor.cov.dates <- factor.cov.dates[factor.cov.dates$sex==2,]
    sex <- "Female"
  }
  
  if (sex == "All") {
    model.cvd.exp <- "with(factor.cov.dates, coxph(Surv(followup, cvd) ~ log(age) + sex + log(tot_chol) + log(hdl_chol) + log(sbp) + log(dbp) + log(glucose) + smoke))"
    #model.cvd.exp <- "with(factor.cov.dates, coxph(Surv(followup, cvd) ~ log(age) + sex + log(tot_chol) + log(hdl_chol) + log(sbp) + log(dbp) + log(glucose) + smoke + log(NK) + log(CD8T) + log(CD4T) + log(Bcell) + log(Mono) + log(Gran) + log(sva1_METH))"
  }
  else if ((sex == "Male") | (sex == "Female")) {
    model.cvd.factor <- "with(factor.cov.dates, coxph(Surv(followup, cvd) ~ log(age) + log(tot_chol) + log(hdl_chol) + log(sbp) + log(dbp) + log(glucose) + smoke))"
    #model.cvd.exp <- "with(factor.cov.dates, coxph(Surv(followup, cvd) ~ log(age) + log(tot_chol) + log(hdl_chol) + log(sbp) + log(dbp) + log(glucose) + smoke + log(NK) + log(CD8T) + log(CD4T) + log(Bcell) + log(Mono) + log(Gran) + log(sva1_METH))"
  }
  
  # Without biomarker
  model.cvd <- eval(parse(text=model.cvd.exp))
  model.cvd
  
  Cstat <- rcorr.cens(-predict(model.cvd), model.cvd$y)
  Cstat <- as.data.frame(t(as.data.frame(Cstat)))
  Cstat$SE <- Cstat$S.D./2
  Cstat$Low95 <- Cstat$`C Index` - 1.96*Cstat$SE
  Cstat$Upper95 <- Cstat$`C Index` + 1.96*Cstat$SE
  
  # With biomarker
  model.cvd.factor <- eval(parse(text=model.cvd.exp))
  model.cvd.factor
  
  Cstat.factor <- rcorr.cens(-predict(model.cvd.factor), model.cvd.factor$y)
  Cstat.factor <- as.data.frame(t(as.data.frame(Cstat.factor)))
  Cstat.factor$SE <- Cstat.factor$S.D./2
  Cstat.factor$Low95 <- Cstat.factor$`C Index` - 1.96*Cstat$SE
  Cstat.factor$Upper95 <- Cstat.factor$`C Index` + 1.96*Cstat$SE
  
  Cstat.results <- rbind(Cstat, Cstat.factor)
  
  #Test de comparaci?n de C-index entre dos modelos (si p-valor < 0.05 es que los modelos discriminan diferente)
  result <- rcorrp.cens(predict(model.cvd), predict(model.cvd.factor), model.cvd$y, method=2)[1:2]
  p.val <- 2*(1-pnorm(as.numeric(abs(result[1]/result[2]))))
  print(p.val)
  
  Cstat.results$p_value <- p.val
  Cstat.results[1,"p_value"] <- NA
  
  write.table(Cstat.results,
              file=paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/Cstat_results_",sex,".csv", sep = ""),
              col.names = T,
              row.names = T,
              sep=",")
  
}

discr(factor = 9)



###############################################################################################
############################################ Reclassification #################################
###############################################################################################

library(nricens)

factor.cov.dates$sex <- as.factor(factor.cov.dates$sex)
factor.cov.dates$smoke <- as.factor(factor.cov.dates$smoke)
factor.cov.dates$followup <- as.numeric(factor.cov.dates$followup)
factor.cov.dates$factor <- as.numeric(factor.cov.dates$factor)

model <- with(factor.cov.dates, 
              coxph(Surv(followup, cvd) ~
                      #sex+
                      log(age)+
                      log(tot_chol)+
                      log(hdl_chol)+
                      log(sbp)+
                      log(dbp)+
                      log(glucose)+
                      smoke,
                    x = T))

model.factor <- with(factor.cov.dates, 
                     coxph(Surv(followup, cvd) ~
                             factor+
                             #sex+
                             log(age)+
                             log(tot_chol)+
                             log(hdl_chol)+
                             log(sbp)+
                             log(dbp)+
                             log(glucose)+
                             smoke),
                     x = T)

range(factor.cov.dates$followup)
length(which(is.na(factor.cov.dates$followup)))

# 5 years
t0 <- 365*5   #1825
p1 <- get.risk.coxph(mdl=model, t0=t0)
p2 <- get.risk.coxph(mdl=model.factor, t0=t0)

range(p1)
range(p2)

set.seed(123)
# nri <- nricens(time = factor.cov.dates$followup,
#                event = factor.cov.dates$cvd,
#                p.std = p1,
#                p.new = p2,
#                #mdl.std = model, 
#                #mdl.new = model.factor,
#                t0 = t0,
#                #updown = "category",
#                updown = "diff",
#                #cut = c(0.05, 0.1), 
#                cut = 0,
#                point.method = "km",
#                niter = 1000, 
#                alpha = 0.05, 
#                msg = T
# )

nri <- nricens(time = factor.cov.dates$followup,
               event = factor.cov.dates$cvd,
               p.std = p1,
               p.new = p2,
               #mdl.std = model,
               #mdl.new = model.factor,
               t0 = t0,
               updown = "category",
               cut = c(0.05, 0.1),
               point.method = "km",
               niter = 1000,
               alpha = 0.05,
               msg = T
)

nri_km <- nri$nri[c(1:3),]

nri_km

write.table(nri_km,
            file=paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/nri_km.csv", sep = ""),
            col.names = T,
            row.names = T,
            sep=",")


#clinical NRI (riesgo intermedio corregido por sesgo clínico)
nricatSurv<-function(resp,time.max,g1,g2,sel=NULL){ 
  cl <- match.call()
  keep<-!is.na(resp) & !is.na(g1) & !is.na(g2)
  resp<-resp[keep]
  g1<-g1[keep]
  g2<-g2[keep]
  if (is.null(sel)) sel<-1:NROW(resp)
  resp<-resp[sel]
  g1<-g1[sel]
  g2<-g2[sel]
  g1<-as.factor(g1)
  g2<-as.factor(g2)
  cens<-resp[,2]
  tt.cens<-table(g1,g2,cens)
  ans<-NULL
  for (i in levels(g1)){
    for (j in levels(g2)) {
      if (sum(g1==i & g2==j,na.rm=TRUE)>0){
        if (tt.cens[i,j,2]==0) # no hi ha morts
          ss<-1
        else{
          ss<-summary(survfit(resp~1,subset=g1==i & g2==j))
          if (all(ss$time>time.max)){
            ss<-1
          }else{
            ss<-ss$surv[ss$time<=time.max]
            ss<-ss[length(ss)]
          }
        }
      } else
        ss<-0
      ans<-c(ans,ss)
    }
  }
  ans<-matrix(ans,nrow=nlevels(g1),byrow=TRUE)
  rownames(ans)<-levels(g1)
  colnames(ans)<-levels(g2)
  tt<-table(g1,g2)
  # nri total
  ttcont<-tt*ans
  ttcas<-tt*(1-ans)
  n.cas<-sum(ttcas)
  n.cont<-sum(ttcont)
  p.up.cont<-sum(ttcont[upper.tri(ttcont)])/n.cont
  p.down.cont<-sum(ttcont[lower.tri(ttcont)])/n.cont
  p.up.cas<-sum(ttcas[upper.tri(ttcas)])/n.cas
  p.down.cas<-sum(ttcas[lower.tri(ttcas)])/n.cas
  beta.cas <- p.up.cas-p.down.cas
  beta.cont <- p.down.cont-p.up.cont
  # nri clinic
  ncateg<-nlevels(g1)
  # observed
  ww<-upper.tri(ttcas)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.up.cas<-sum(ttcas[ww])/sum(ttcas[-c(1,ncateg),])
  ww<-lower.tri(ttcas)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.down.cas<-sum(ttcas[ww])/sum(ttcas[-c(1,ncateg),])
  ww<-upper.tri(ttcont)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.up.cont<-sum(ttcont[ww])/sum(ttcont[-c(1,ncateg),])
  ww<-lower.tri(ttcont)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.down.cont<-sum(ttcont[ww])/sum(ttcont[-c(1,ncateg),])
  beta.int.cas <- p.up.cas-p.down.cas
  beta.int.cont <- p.down.cont-p.up.cont
  # expected
  ttcas.exp<-(t(ttcas)+ttcas)/2
  ttcont.exp<-(t(ttcont)+ttcont)/2
  ww<-upper.tri(ttcas.exp)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.up.cas<-sum(ttcas.exp[ww])/sum(ttcas.exp[-c(1,ncateg),])
  ww<-lower.tri(ttcas.exp)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.down.cas<-sum(ttcas.exp[ww])/sum(ttcas.exp[-c(1,ncateg),])
  ww<-upper.tri(ttcont.exp)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.up.cont<-sum(ttcont.exp[ww])/sum(ttcont.exp[-c(1,ncateg),])
  ww<-lower.tri(ttcont.exp)
  ww[1,]<-ww[ncateg,]<-FALSE
  p.down.cont<-sum(ttcont.exp[ww])/sum(ttcont.exp[-c(1,ncateg),])
  beta.int.exp.cas <- p.up.cas-p.down.cas
  beta.int.exp.cont <- p.down.cont-p.up.cont
  beta.int.corr.cas <- beta.int.cas-beta.int.exp.cas
  beta.int.corr.cont <- beta.int.cont-beta.int.exp.cont
  # out
  list(beta=beta.cas+beta.cont,beta.cas=beta.cas,beta.cont=beta.cont,
       beta.int=beta.int.cas+beta.int.cont,beta.int.cas=beta.int.cas,beta.int.cont=beta.int.cont,
       beta.int.corr=beta.int.corr.cas+beta.int.corr.cont,beta.int.corr.cas=beta.int.corr.cas,beta.int.corr.cont=beta.int.corr.cont,
       ttcont=ttcont,ttcas=ttcas,tt=tt,surv=ans,call=cl,n=NROW(resp))
  
}
# function to compute confidence intervals for NRIcat. To be applied to an object computed by nricatSurv
nricatSurvci<-function(x,B=300,digits=2, timesele, verbose=FALSE){
  z<-qnorm(1-0.05/2)
  n<-x$n  
  iii<-0  
  estim<-unlist(x[1:9])
  sim<-replicate(B,{
    if (verbose) print(iii<<-iii+1)
    unlist(update.default(x,sel=sample(1:n,replace=TRUE))[1:9])
  })
  se<-apply(sim,1,sd)
  pval<-2*(1-pnorm(abs(estim/se)))
  ans<-cbind(estim,estim-z*se,estim+z*se,pval)
  colnames(ans)<-c("estim","lower","upper","pval")
  rownames(ans)<-sub("beta","nri",rownames(ans))
  ans  
}

model <- with(factor.cov.dates, 
              coxph(Surv(followup, cvd) ~
                      #sex+
                      log(age)+
                      log(tot_chol)+
                      log(hdl_chol)+
                      log(sbp)+
                      log(dbp)+
                      log(glucose)+
                      smoke,
                    x = T))

model.factor <- with(factor.cov.dates, 
                     coxph(Surv(followup, cvd) ~
                             factor+
                             #sex+
                             log(age)+
                             log(tot_chol)+
                             log(hdl_chol)+
                             log(sbp)+
                             log(dbp)+
                             log(glucose)+
                             smoke),
                     x = T)

range(factor.cov.dates$followup)
length(which(is.na(factor.cov.dates$followup)))

# 5 years
t0 <- 365*5   #1825
p1 <- get.risk.coxph(mdl=model, t0=t0)
p2 <- get.risk.coxph(mdl=model.factor, t0=t0)

range(p1)
range(p2)

# puntos de corte de los riesgos estimados (5 a?os, 0.05 y 0.10)
talls <- c(0.05, 0.1)

# categorizamos los riesgos y creamos los grupos seg?n los puntos de corte determinados (talls)
g1 <- cut(p1, c(-Inf, talls, Inf))
g2 <- cut(p2, c(-Inf, talls, Inf))


set.seed(123)
# cálculo del NRI categórico (la función siguiente da el NRI total y también el intermedio sin sesgo).
nricat <- nricatSurvci(nricatSurv(Surv(factor.cov.dates$followup, factor.cov.dates$cvd), t0, g1, g2),
                       timesele=t0,
                       verbose=FALSE,
                       B=300)
nricat[,1:3] <- nricat[,1:3]*100
nricat

write.table(nricat,
            file=paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/nri_clinical.csv", sep = ""),
            col.names = T,
            row.names = T,
            sep=",")


date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: COX SURVIVAL SCRIPT ENDS ######################## \n\n", file = stdout())


