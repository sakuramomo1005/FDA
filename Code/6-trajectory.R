setwd(store_doc)
load('data_trans.RData')
head(data_trans)

### cluster 1 trajectory
cluster1 = data_trans[data_trans$cluster == 1, ]$subj
dat1 = dat[dat$subj %in% cluster1, ]

# draw the plots
placebo <- NULL
prozac <- NULL
imi <- NULL
for (jt in unique(dat1$trt)){ # fit lme for each arm
  dati <- dat1[dat1$trt==jt,]
  fit1 <- lmer(y ~ t1 + I(t1^2) + (t1+I(t1^2)|subj), data=dati)
  D <- as.matrix(VarCorr(fit1)$subj) # Covariance matrix for random effects
  D <- D[1:p,1:p]
  beta <- as.matrix(fixef(fit1))
  sigma <- attr(VarCorr(fit1), "sc")
  bis <- as.matrix(coef(fit1)$subj)%*%t(solve(A))
  
  responder <- NULL # record subjects that are responders or not
  age <- NULL
  BaselineCGI <- NULL
  for (isubj in unique(dati$subj)){
    datisubj <- dati[dati$subj ==isubj,]
    responder <- rbind(responder, datisubj$responder[1])
    age <- rbind(age, datisubj$age[1])
    BaselineCGI <- rbind(BaselineCGI, datisubj$BaselineCGI[1])
  }
  
  if (jt == 0){
    placebo$n <- length(unique(dati$subj))
    placebo$dat <- dati
    placebo$D <- D
    placebo$beta <- beta
    placebo$sigma <- sigma
    placebo$bis <- bis
    placebo$responder <- responder
    placebo$age <- age
    placebo$BaselineCGI <- BaselineCGI
    placebo$fit <- fit1
  }
  if (jt == 1){
    prozac$n <- length(unique(dati$subj))
    prozac$dat <- dati
    prozac$D <- D
    prozac$beta <- beta
    prozac$sigma <- sigma
    prozac$bis <- bis
    prozac$responder <- responder
    prozac$age <- age
    prozac$BaselineCGI <- BaselineCGI
    prozac$fit <- fit1
  }
  if (jt == 2){
    imi$n <- length(unique(dati$subj))
    imi$dat <- dati
    imi$D <- D
    imi$beta <- beta
    imi$sigma <- sigma
    imi$bis <- bis
    imi$responder <- responder
    imi$age <- age
    imi$BaselineCGI <- BaselineCGI
    imi$fit <- fit1
  }
}

png('cluster1 trajectory plot.png')
par(mfrow = c(1,2))
nf <- layout(matrix(c(0,0,1,2,1,2,1,2,1,2,1,2,1,2,1,2,0,0),9,2,byrow=TRUE))
layout.show(nf)
tplot=seq(0,6, by=.1)

plot(t, t*5.7, type="n", main="Figure 1: Fluoxetine and Placebo Treated Subjects",
     xlab="Week", ylab="HRSD")
for (i in 1:prozac$n){
  bi=A%*%as.matrix(prozac$bis[i,1:3])
  lines(tplot, bi[1,1]+bi[2,1]*tplot+bi[3,1]*tplot^2, col=4, lwd=1.3)
}
legend("bottomleft", c("Drug", "Placebo"), col=c(4,2), lwd=c(2,2))

plot(t, t*5.7, type="n", main="Figure 1: Fluoxetine and Placebo Treated Subjects",
     xlab="Week", ylab="HRSD")
for (i in 1:placebo$n){
  bi=A%*%as.matrix(placebo$bis[i,])
  lines(tplot, bi[1,1]+bi[2,1]*tplot+bi[3,1]*tplot^2, col=2)
}
legend("bottomleft", c("Drug", "Placebo"), col=c(4,2), lwd=c(2,2))
dev.off()


### cluster 4 trajectory
cluster4 = data_trans[data_trans$cluster == 4, ]$subj
dat4 = dat[dat$subj %in% cluster4, ]

# draw the plots
placebo <- NULL
prozac <- NULL
imi <- NULL
for (jt in unique(dat4$trt)){ # fit lme for each arm
  dati <- dat4[dat4$trt==jt,]
  fit1 <- lmer(y ~ t1 + I(t1^2) + (t1+I(t1^2)|subj), data=dati)
  D <- as.matrix(VarCorr(fit1)$subj) # Covariance matrix for random effects
  D <- D[1:p,1:p]
  beta <- as.matrix(fixef(fit1))
  sigma <- attr(VarCorr(fit1), "sc")
  bis <- as.matrix(coef(fit1)$subj)%*%t(solve(A))
  
  responder <- NULL # record subjects that are responders or not
  age <- NULL
  BaselineCGI <- NULL
  for (isubj in unique(dati$subj)){
    datisubj <- dati[dati$subj ==isubj,]
    responder <- rbind(responder, datisubj$responder[1])
    age <- rbind(age, datisubj$age[1])
    BaselineCGI <- rbind(BaselineCGI, datisubj$BaselineCGI[1])
  }
  
  if (jt == 0){
    placebo$n <- length(unique(dati$subj))
    placebo$dat <- dati
    placebo$D <- D
    placebo$beta <- beta
    placebo$sigma <- sigma
    placebo$bis <- bis
    placebo$responder <- responder
    placebo$age <- age
    placebo$BaselineCGI <- BaselineCGI
    placebo$fit <- fit1
  }
  if (jt == 1){
    prozac$n <- length(unique(dati$subj))
    prozac$dat <- dati
    prozac$D <- D
    prozac$beta <- beta
    prozac$sigma <- sigma
    prozac$bis <- bis
    prozac$responder <- responder
    prozac$age <- age
    prozac$BaselineCGI <- BaselineCGI
    prozac$fit <- fit1
  }
  if (jt == 2){
    imi$n <- length(unique(dati$subj))
    imi$dat <- dati
    imi$D <- D
    imi$beta <- beta
    imi$sigma <- sigma
    imi$bis <- bis
    imi$responder <- responder
    imi$age <- age
    imi$BaselineCGI <- BaselineCGI
    imi$fit <- fit1
  }
}

png('cluster4 trajectory plot.png')
par(mfrow = c(1,2))
nf <- layout(matrix(c(0,0,1,2,1,2,1,2,1,2,1,2,1,2,1,2,0,0),9,2,byrow=TRUE))
layout.show(nf)
tplot=seq(0,6, by=.1)

plot(t, t*5.7, type="n", main="Figure 1: Fluoxetine and Placebo Treated Subjects",
     xlab="Week", ylab="HRSD")
for (i in 1:prozac$n){
  bi=A%*%as.matrix(prozac$bis[i,1:3])
  lines(tplot, bi[1,1]+bi[2,1]*tplot+bi[3,1]*tplot^2, col=4, lwd=1.3)
}
legend("bottomleft", c("Drug", "Placebo"), col=c(4,2), lwd=c(2,2))

plot(t, t*5.7, type="n", main="Figure 1: Fluoxetine and Placebo Treated Subjects",
     xlab="Week", ylab="HRSD")
for (i in 1:placebo$n){
  bi=A%*%as.matrix(placebo$bis[i,])
  lines(tplot, bi[1,1]+bi[2,1]*tplot+bi[3,1]*tplot^2, col=2)
}
legend("bottomleft", c("Drug", "Placebo"), col=c(4,2), lwd=c(2,2))
dev.off()

par(mfrow = c(1,1))

setwd(load_doc)