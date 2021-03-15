rm(list=ls())

library(survival)
library(survAUC)
library(Hmisc)

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

groups <- "nogroups"
samples <- "391"
# 9, 19, 21, 27
factor <- 21
factor.cov.dates <- read.csv(paste("MOFA/3.downstream_analysis/regicor/prediction/factor_cov/factor",factor,"_cov_dates.csv", sep = ""), sep = ",", stringsAsFactors = FALSE)
head(factor.cov.dates)
dim(factor.cov.dates)

## Stratification by sex

# Male
# factor.cov.dates <- factor.cov.dates[factor.cov.dates$sex==1,]
# sex <- "Male"
# # Female
# factor.cov.dates <- factor.cov.dates[factor.cov.dates$sex==2,]
# sex <- "Female"

##########################################################################################
############################## Mediana de seguimiento ####################################
##########################################################################################

# cvd
iqr_cvd <- (IQR(factor.cov.dates$followup))/365.25
median_cvd <- (median(factor.cov.dates$followup, na.rm=T))/365.25
median_cvd
iqr_cvd_low <- median_cvd-iqr_cvd
iqr_cvd_up <- median_cvd+iqr_cvd
iqr_cvd_low
iqr_cvd_up

# no cvd
#df.nocvd <- factor.cov.dates[factor.cov.dates$cvd == 0,]


##################################################################################################
### Association between biomarker (factor X) and incidence of CVD (we need the time to event) ###
##################################### COX REGRESSION #############################################
##################################################################################################

glm.res <- glm(cvd ~ factor, data = factor.cov.dates, family = "binomial")
summary(glm.res)
# El t-test debería dar lo mismo que el hecho en statistical_tests!
t.test(factor ~ cvd, data = factor.cov.dates)

mod.surv.cvd.1 <- coxph(Surv(followup, cvd) ~
                          factor,
                        data=factor.cov.dates)
cat("\n\nResults for model 1: factor\n\n", file = stdout())
mod.surv.cvd.1

# Add Confidence Intervals
mod.surv.cvd.1.df <- as.data.frame(summary(mod.surv.cvd.1)$coefficients)
mod.surv.cvd.1.df$lower95 <- as.numeric(exp(confint(mod.surv.cvd.1)))[1]
mod.surv.cvd.1.df$upper95 <- as.numeric(exp(confint(mod.surv.cvd.1)))[2]
names(mod.surv.cvd.1.df) <- c("coef","HR","se(coef)","z","P", "lower95", "upper95")
mod.surv.cvd.1.df$model <- "Model 1"
mod.surv.cvd.1.df$sex <- sex
mod.surv.cvd.1.df$factor <- paste("Factor ",factor,sep = "")
mod.surv.cvd.1.df

# Añadir interacción factor * sex en ls factores donde haya cierta correlación con sexo
mod.surv.cvd.2 <- coxph(Surv(followup, cvd) ~
                          factor+
                          sex+
                          age+
                          factor*sex
                          ,
                          data=factor.cov.dates)

cat("\n\nResults for model 2: factor + sex + age\n\n", file = stdout())
mod.surv.cvd.2

mod.surv.cvd.2.df <- as.data.frame(summary(mod.surv.cvd.2)$coefficients)
mod.surv.cvd.2.df
mod.surv.cvd.2.df$lower95 <- as.numeric(exp(confint(mod.surv.cvd.2)))[1:2] #[1:3]
mod.surv.cvd.2.df$upper95 <- as.numeric(exp(confint(mod.surv.cvd.2)))[3:4] #[4:6]
names(mod.surv.cvd.2.df) <- c("coef","HR","se(coef)","z","P", "lower95", "upper95")
mod.surv.cvd.2.df$model <- "Model 2"
mod.surv.cvd.2.df$sex <- sex
mod.surv.cvd.2.df$factor <- paste("Factor ",factor,sep = "")
mod.surv.cvd.2.df

mod.surv.cvd.3 <- coxph(Surv(followup, cvd) ~
                          factor+
                          sex+
                          age+
                          tot_chol+
                          hdl_chol+
                          sbp+
                          dbp+
                          glucose+
                          smoke+
                          factor*sex
                          ,
                        data=factor.cov.dates)

cat("\n\nResults for model 3: factor + sex + age + tot_chol + hdl_chol + sbp + dbp + smoke + glucose\n\n", file = stdout())
mod.surv.cvd.3

# Add Confidence Intervals
mod.surv.cvd.3.df <- as.data.frame(summary(mod.surv.cvd.3)$coefficients)
mod.surv.cvd.3.df
mod.surv.cvd.3.df$lower95 <- as.numeric(exp(confint(mod.surv.cvd.3)))[1:dim(mod.surv.cvd.3.df)[1]]
mod.surv.cvd.3.df$upper95 <- as.numeric(exp(confint(mod.surv.cvd.3)))[(dim(mod.surv.cvd.3.df)[1]+1):(dim(mod.surv.cvd.3.df)[1]*2)]
names(mod.surv.cvd.3.df) <- c("coef","HR","se(coef)","z","P", "lower95", "upper95")
mod.surv.cvd.3.df$model <- "Model 3"
mod.surv.cvd.3.df$sex <- sex
mod.surv.cvd.3.df$factor <- paste("Factor ",factor,sep = "")
mod.surv.cvd.3.df

#One for each factor
cox.results.f9 <- rbind(mod.surv.cvd.1.df,mod.surv.cvd.2.df[1,],mod.surv.cvd.3.df[1,])
cox.results.f19 <- rbind(mod.surv.cvd.1.df,mod.surv.cvd.2.df[1,],mod.surv.cvd.3.df[1,])
cox.results.f21.M <- rbind(mod.surv.cvd.1.df,mod.surv.cvd.2.df[1,],mod.surv.cvd.3.df[1,])
cox.results.f21.F <- rbind(mod.surv.cvd.1.df,mod.surv.cvd.2.df[1,],mod.surv.cvd.3.df[1,])
cox.results.f27 <- rbind(mod.surv.cvd.1.df,mod.surv.cvd.2.df[1,],mod.surv.cvd.3.df[1,])

cox.results <- rbind(cox.results.f9,cox.results.f19,cox.results.f21.M,cox.results.f21.F,cox.results.f27)
rownames(cox.results) <- NULL
cox.results$P_adjusted <- NA

# FDR adjust by model
cox.results[cox.results$model == "Model 1",]$P_adjusted <- p.adjust(cox.results[cox.results$model == "Model 1",]$P, method = "fdr")
cox.results[cox.results$model == "Model 2",]$P_adjusted <- p.adjust(cox.results[cox.results$model == "Model 2",]$P, method = "fdr")
cox.results[cox.results$model == "Model 3",]$P_adjusted <- p.adjust(cox.results[cox.results$model == "Model 3",]$P, method = "fdr")
cox.results

write.csv(cox.results, file = paste("MOFA/3.downstream_analysis/regicor/prediction/cox_results.csv", sep = ""),
          col.names = TRUE, row.names = FALSE, sep = ",")

mod.surv.cvd.4 <- coxph(Surv(followup, cvd) ~
                          pspline(factor)+
                          sex+
                          age+
                          tot_chol+
                          hdl_chol+
                          sbp+
                          dbp+
                          glucose+
                          smoke,
                        data=factor.cov.dates)

termplot(mod.surv.cvd.4, terms = "pspline(factor)", se = TRUE, rug = TRUE)
mod.surv.cvd.4

### kaplan meier (can only be done for categorical data...)
## tres perfiles: 1-33, 33-66, 66-100

summary(factor.cov.dates$factor)
factor.quantiles <- quantile(factor.cov.dates$factor, probs = seq(0, 1, 0.33))
factor.categorical <- ifelse(factor.cov.dates$factor < factor.quantiles[2], 1, factor.cov.dates$factor)
factor.categorical <- ifelse(factor.categorical < factor.quantiles[3], 2, factor.categorical)
factor.categorical <- ifelse(factor.categorical == 1 | factor.categorical == 2, factor.categorical, 3)
factor.cov.dates$factor_categorical <- factor.categorical

survObject <- Surv(factor.cov.dates$followup, factor.cov.dates$cvd)
#fit <- survfit(survObject ~ 1, data = factor.cov.dates)
fit <- survfit(survObject ~ factor_categorical, data = factor.cov.dates)
plot(fit)

# significativo?
logRank <- survdiff(survObject ~ factor_categorical, data = factor.cov.dates)
pval <- p.val <- 1 - pchisq(logRank$chisq, length(logRank$n) - 1) 
pval

library(survminer)
ggsurvplot(fit, data=factor.cov.dates, pval = pval)

###############################################################################################
############################################ Discriminaci?n ###################################
###############################################################################################

# Without biomarker
model.cvd <- with(factor.cov.dates, 
              coxph(Surv(followup, cvd) ~
                      sex+
                      log(age)+
                      log(tot_chol)+
                      log(hdl_chol)+
                      log(sbp)+
                      log(dbp)+
                      log(glucose)+
                      smoke))

Cstat <- rcorr.cens(-predict(model.cvd), model.cvd$y)
Cstat <- as.data.frame(t(as.data.frame(Cstat)))
Cstat$SE <- Cstat$S.D./2
Cstat$Low95 <- Cstat$`C Index` - 1.96*Cstat$SE 
Cstat$Upper95 <- Cstat$`C Index` + 1.96*Cstat$SE 

# With biomarker
model.cvd.factor <- with(factor.cov.dates, 
              coxph(Surv(followup, cvd) ~
                      sex+
                      log(age)+
                      log(tot_chol)+
                      log(hdl_chol)+
                      log(sbp)+
                      log(dbp)+
                      log(glucose)+
                      smoke+
                      factor)) ## PONGO LOG O NO? VALORES NEGATIVOS!

Cstat.factor <- rcorr.cens(-predict(model.cvd.factor), model.cvd.factor$y)
Cstat.factor <- as.data.frame(t(as.data.frame(Cstat.factor)))
Cstat.factor$SE <- Cstat.factor$S.D./2
Cstat.factor$Low95 <- Cstat.factor$`C Index` - 1.96*Cstat$SE 
Cstat.factor$Upper95 <- Cstat.factor$`C Index` + 1.96*Cstat$SE 

Cstat.results <- rbind(Cstat, Cstat.factor)

#Test de comparaci?n de C-index entre dos modelos (si p-valor < 0.05 es que los modelos discriminan diferente)
result <- rcorrp.cens(predict(model.cvd), predict(model.cvd.factor), model.cvd$y, method=2)[1:2]
p.val <- 2*(1-pnorm(as.numeric(abs(result[1]/result[2]))))
p.val

Cstat.results$p_value <- p.val
Cstat.results[1,"p_value"] <- NA

write.table(Cstat.results, 
            file=paste(groups,"_",samples,"/Cstat_results.csv", sep = ""),
            col.names = T,
            row.names = T,
            sep=",")

###############################################################################################
############################################ Reclasificaci?n ##################################
###############################################################################################

library(nricens)

factor.cov.dates$sex <- as.factor(factor.cov.dates$sex)
factor.cov.dates$smoke <- as.factor(factor.cov.dates$smoke)
factor.cov.dates$followup <- as.numeric(factor.cov.dates$followup)
factor.cov.dates$factor <- as.numeric(factor.cov.dates$factor)

model <- with(factor.cov.dates, 
              coxph(Surv(followup, cvd) ~
                      sex+
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
                               sex+
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

# con 5 a?os
t0 <- 365*5   #1825
p1 <- get.risk.coxph(mdl=model, t0=t0)
p2 <- get.risk.coxph(mdl=model.factor, t0=t0)

range(p1)
range(p2)

set.seed(123)
nri <- nricens(time = factor.cov.dates$followup,
               event = factor.cov.dates$cvd,
               p.std = p1,
               p.new = p2,
               #mdl.std = model, 
               #mdl.new = model.factor,
               t0 = t0,
               #updown = "category",
               updown = "diff",
               #cut = c(0.05, 0.1), 
               cut = 0,
               point.method = "km",
               niter = 1000, 
               alpha = 0.05, 
               msg = T
               )

nri_km <- nri$nri[c(1:3),]
nri_km
write.table(nri_km,
            file=paste(groups,"_",samples,"/nri_km.csv", sep = ""),
            col.names = T,
            row.names = T,
            sep=",")



#NRI cl?nico (riesgo intermedio corregido por sesgo cl?nico)
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
                      sex+
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
                               sex+
                               log(age)+
                               log(tot_chol)+
                               log(hdl_chol)+
                               log(sbp)+
                               log(dbp)+
                               log(glucose)+
                               smoke),
                       x = T)

range(db$cvd_date)
length(which(is.na(db$cvd_date)))
# [1] 0

# con 5 a?os
t0 <- 365.25*5   #1825
p1 <- get.risk.coxph(mdl=model1, t0=t0)
p2 <- get.risk.coxph(mdl=model2, t0=t0)

range(p1)
range(p2)

# p1 riesgos estimados con el modelo de referencia
# p2: riesgos estimados con el modelo con las variables gen?ticas
# times: variable tiempo hasta el evento
# cens: variable evento

# puntos de corte de los riesgos estimados (5 a?os, 0.05 y 0.10)
talls <- c(0.05, 0.1)

# categorizamos los riesgos y creamos los grupos seg?n los puntos de corte determinados (talls)
g1 <- cut(p1, c(-Inf, talls, Inf))
g2 <- cut(p2, c(-Inf, talls, Inf))

# c?lculo del NRI categ?rico (la funci?n siguiente da el NRI total y tambi?n el intermedio sin sesgo).
nricat <- nricatSurvci(nricatSurv(Surv(factor.cov.dates$followup, factor.cov.dates$cvd), t0, g1, g2),
                       timesele=t0,
                       verbose=FALSE,
                       B=300)
nricat[,1:3] <- nricat[,1:3]*100

nricat

