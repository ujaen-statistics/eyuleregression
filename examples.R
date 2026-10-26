# ------------------------------------------------------------------------------
# APPLICATION EXAMPLES OF THE PAPER
#-------------------------------------------------------------------------------

# Required packages and R. files
# ------------------------------

packages <- c("COUNT", "GWRM", "DGLMExtPois", "pracma", "Rchoice", "doSNOW", 
              "foreach", "xtable")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(packages, require, character.only = TRUE))

archivos <- c("utils", "glm_EY", "summary_EY", "lrtest", "plots",
              "distribution_EY", "residuals", "glm_EY_methods", "expected",
              "gof")

invisible(lapply(paste0(archivos, ".R"), source))

# ------------------------------------------------------------------------------
# Takeover bids data
# ------------------------------------------------------------------------------

Bids$size.sq <- Bids$size^2
formula <- numbids ~ leglrest + rearest + finrest + whtknght + bidprem + insthold + size + size.sq + regulatn

# EYule regression model
# ----------------------

# Fitted model (Table 3)
ey_bids <- glm.EY(formula, numbids ~ 1, data = Bids, maxiter_series = 1000000, tol = 0.00001)
summary(ey_bids)

# Goodness-of-fit measures (Table 3)

## AIC
ey_bids$aic
## BIC
ey_bids$bic
## Diff
expec_ey <- EY_expected(ey_bids)
(Diff <- sum(abs(expec_ey$frequencies - expec_ey$observed_freq)))
## chi-squared
expec_ey$chi2

# Diagnosis plots

## Figure 3
par(mfrow = c(1,2))
set.seed(1234)
plot(ey_bids)
## Figure 4
par(mfrow = c(1,1))
set.seed(1234)
r1 <- residuals(ey_bids, envelope = T)
envelope_down <- apply(r1$sim_residuals, 1, min)
envelope_up <- apply(r1$sim_residuals, 1, max)
weirds <- which(r1$residuals < envelope_down|r1$residuals > envelope_up)
score_weirds <- qnorm(weirds / 127)
points(score_weirds,r1$residuals[weirds], col = "red")

# Validation of EYule regression model (LOOCV)

num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

s_bids <- nrow(Bids)

pb <- txtProgressBar(max = s_bids, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

y_jack_list_bids <- foreach(k = 1:s_bids, .errorhandling = "pass", .options.snow = opts) %dopar% {
  datosk <- Bids[-k, ]
  fit_k <- glm.EY(formula, numbids ~ 1, data = datosk, maxiter_series = 1000000, tol = 0.00001)
  pred_k <- predict.glm_EY(fit_k, newdata = Bids, type = "response")[k]
  return(as.numeric(pred_k))
}

close(pb)
stopCluster(cl)

y_jack_bids <- sapply(y_jack_list_bids, function(x) if (is.numeric(x) && length(x) == 1) x else NA)
mse_jack_bids <- sum((Bids$numbids - y_jack_bids)^2) / length(Bids$numbids)
mse_jack_bids

MSPE(ey_bids)

#Number of under-dispersed cases
mu_bids <- means_EY(ey_bids$lambdas, ey_bids$gammas, maxiter_series = 1000000)
s2_bids <- variances_EY(ey_bids$lambdas, ey_bids$gammas, maxiter_series = 1000000)
sum(mu_bids > s2_bids)

#Table 4
under <- mu_bids > s2_bids
summary(Bids[under, c(3,5:12)])
porc_under <- apply(Bids[under, 8:12], 2, table) / 102
mean_under <- c(apply(Bids[under, c(3, 5:7)], 2, mean), porc_under[2,])
var_under <- c(apply(Bids[under, c(3, 5:7)], 2, var), rep("-",5))

over <- mu_bids < s2_bids
summary(Bids[over, c(3, 5:12)])
porc_over <- apply(Bids[over, 8:12], 2, table) / 24
mean_over <- c(apply(Bids[over, c(3, 5:7)], 2, mean), porc_over[2,])
var_over <- c(apply(Bids[over,c(3, 5:7)], 2, var), rep("-",5))

desc_ey <- data.frame(Mean_under = mean_under, Variances_under = var_under, 
                    Mean_over = mean_over, Variances_over = var_over)

xtable(desc_ey, caption="Descriptive statistics for under- and over- dispersed cases", align = c("r","c","c","c","c"), digits = 4)


# hP regression model
# -------------------

# Fitted model (Table 3)
hp_bids <- glm.hP(formula, numbids ~ 1, data = Bids)
summary(hp_bids)

# Goodness-of-fit measures (Table 3)

## AIC
hp_bids$aic
## BIC
hp_bids$bic
## Diff
expec_hp <- hP_expected(hp_bids)
(Diff <- sum(abs(expec_hp$frequencies - expec_hp$observed_freq)))
## chi-squared
expec_hp$chi2

# CMP regression model
# --------------------

# Fitted model (Table 3)
cmp_bids <- glm.CMP(formula, numbids ~ 1, data = Bids)
summary(cmp_bids)

# Goodness-of-fit measures (Table 3)

## AIC
cmp_bids$aic
## BIC
cmp_bids$bic
## Diff
expec_cmp <- CMP_expected(cmp_bids)
(Diff <- sum(abs(expec_cmp$frequencies - expec_cmp$observed_freq)))
## chi-squared
expec_cmp$chi2

# CTP regression model
# --------------------

# Download and source the relevant R scripts from GitHub: https://github.com/ujaen-statistics/ctpregression

# Fitted model (Table 3)
ctp_bids <- ctp.fit(formula, data = Bids, astart = -1.25, bstart = 1.51)
summary.ctp(ctp_bids)

# Number of under-dispersed cases
sum(ctp_bids$parameters[1] < (-ctp_bids$fitted.values-1) / 2)

# Goodness-of-fit measures (Table 3)

## AIC
ctp_bids$aic
## BIC
ctp_bids$bic
## Diff
expec_ctp <- expected.ctp(ctp_bids)
(Diff <- sum(abs(expec_ctp$frequencies - expec_ctp$observed_freq)))
## chi-squared
expec_ctp$chi2

# Figure 2
par(mfrow = c(1,1))
freq_agr <- matrix(c(expec_ey$observed_freq, expec_ey$frequencies, expec_hp$frequencies, expec_cmp$frequencies,
                   expec_ctp$frequencies), byrow=T, nrow=5)
colnames(freq_agr) <- 0:11
rownames(freq_agr) <- c("Observed", "EY expected", "hP expected", "CMP expected", "CTP expected")
color.names <- c("black","grey75","grey50","grey25","white")
barplot(freq_agr,beside = T, xlab = "Y", ylab = "Frequencies", col = color.names)
legend("topright",rownames(freq_agr), cex=0.9, fill = color.names, bty = "n")


# ------------------------------------------------------------------------------
# Articles Data
# ------------------------------------------------------------------------------

data(Articles)
Articles$fem <- as.factor(Articles$fem)
Articles$mar <- as.factor(Articles$mar)
formula <- art ~ fem + mar + kid5 + phd + ment

# EYule regression model
# ----------------------

# Fitted model (Table 5)
formula1 <- art ~ ment
ey_articles <- glm.EY(formula, formula1, data = Articles, maxiter_series = 1000000)
summary(ey_articles)

# Goodness-of-fit measures (Table 5)

## AIC
ey_articles$aic
## BIC
ey_articles$bic
## Diff
expec_ey <- EY_expected(ey_articles)
(Diff <- sum(abs(expec_ey$frequencies - expec_ey$observed_freq)))
## chi-squared
expec_ey$chi2

# Diagnosis plots

## Figure 5
par(mfrow = c(1,2))
set.seed(1234)
plot(ey_articles)
## Figure 7
par(mfrow = c(1,1))
set.seed(1234)
r1 <- residuals(ey_articles, envelope = T)
envelope_down <- apply(r1$sim_residuals, 1, min)
envelope_up <- apply(r1$sim_residuals, 1, max)
weirds <- which(r1$residuals < envelope_down|r1$residuals > envelope_up)
score_weirds <- qnorm(weirds / 916)
points(score_weirds,r1$residuals[weirds], col = "red")

# Validation of the EYule regression model (LOOCV)
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

s_art <- nrow(Articles)

pb <- txtProgressBar(max = s_art, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

y_jack_list_art <- foreach(k = 1:s_art, .errorhandling = "pass", .options.snow = opts) %dopar% {
  datosk <- Articles[-k, ]
  fit_k <- glm.EY(formula, formula1, data = datosk,
                  maxiter_series = 1000000, tol = 0.000001)
  pred_k <- predict.glm_EY(fit_k, newdata = Articles, type = "response")[k]
  return(as.numeric(pred_k))
}

close(pb)
stopCluster(cl)

y_jack_art <- sapply(y_jack_list_art, function(x) if (is.numeric(x) && length(x) == 1) x else NA)
mse_jack_art <- sum((Articles$art - y_jack)^2) / length(Articles$art)
mse_jack_art

MSPE(ey_articles)

# Figure 8
mu_art <- means_EY(ey_articles$lambdas, ey_articles$gammas, maxiter_series = 1000000)
s2_art <- variances_EY(ey_articles$lambdas, ey_articles$gammas, maxiter_series = 1000000)
ai <- s2_art / mu_art
par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(4.1,4.1,1.1,1.1))
plot(ai ~ Articles$ment, ylim=c(0,6), ylab = "AI", xlab = "ment")
abline(1, 0, lty=2)
plot(ai ~ Articles$phd, ylim=c(0,6), ylab = "AI", xlab = "phd")
abline(1, 0, lty = 2)
boxplot(ai ~ Articles$kid5, ylim = c(0,6), ylab = "AI",xlab = "kid5")
abline(1, 0, lty = 2)
boxplot(ai ~ Articles$fem, ylim = c(0,6), ylab = "AI", xlab = "fem")
abline(1, 0, lty = 2)


# hP regression model
# -------------------

# Fitted model (Table 5)
hp_articles <- glm.hP(formula, ~ phd, data = Articles)
summary(hp_articles)

# Goodness-of-fit measures (Table 5)

## AIC
hp_articles$aic
## BIC
hp_articles$bic
## Diff
expec_hp <- hP_expected(hp_articles)
(Diff <- sum(abs(expec_hp$frequencies - expec_hp$observed_freq)))
## chi-squared
expec_hp$chi2


# CMP regression model
# --------------------

# Fitted model (Table 5)
cmp_articles <- glm.CMP(formula, ~ phd, data = Articles)
summary(cmp_articles)

# Goodness-of-fit measures (Table 5)

## AIC
cmp_articles$aic
## BIC
cmp_articles$bic
## Diff
expec_cmp <- CMP_expected(cmp_articles)
(Diff <- sum(abs(expec_cmp$frequencies - expec_cmp$observed_freq)))
## chi-squared
expec_cmp$chi2


# CTP regression model
# --------------------

# Fitted model (Table 5)
ctp_articles <- ctp.fit(formula, astart = 0.5, bstart = 1, data = Articles)
summary.ctp(ctp_articles)

# Number of under-dispersed cases
sum(ctp_articles$parameters[1] < (-ctp_articles$fitted.values-1) / 2)

# Goodness-of-fit measures (Table 5)

## AIC
ctp_articles$aic
## BIC
ctp_articles$bic
## Diff
expec_ctp <- expected.ctp(ctp_articles)
(Diff <- sum(abs(expec_ctp$frequencies - expec_ctp$observed_freq)))
## chi-squared
expec_ctp$chi2

# Figure 6
par(mfrow = c(1,1))
freq_agr <- matrix(c(expec_ey$observed_freq, expec_ey$frequencies, expec_hp$frequencies, expec_cmp$frequencies, expec_ctp$frequencies), byrow = T, nrow = 5)
colnames(freq_agr) <- 0:20
rownames(freq_agr)<-c("Observed", "EY expected", "hP expected", "CMP expected", "CTP expected")
color.names <- c("black", "grey", "grey66", "grey33", "white")
barplot(freq_agr,beside = T, xlab = "Y", ylab = "Frequencies", col = color.names)
legend("topright",rownames(freq_agr), cex=0.9, fill = color.names, bty = "n")

# Figure 9
mu_hp<-DGLMExtPois:::means_hp(hp_articles$lambdas,hp_articles$gammas)
s2_hp<-DGLMExtPois:::variances_hp(hp_articles$lambdas,hp_articles$gammas)
ai_hp<-s2_hp/mu_hp
mu_cmp<-DGLMExtPois:::means_cmp(cmp_articles$lambdas,cmp_articles$nus)
s2_cmp<-DGLMExtPois:::variances_cmp(cmp_articles$lambdas,cmp_articles$nus)
ai_cmp<-s2_cmp/mu_cmp

par(mfrow = c(1,2), mar=c(4.1,4.1,4.1,1.1))
plot(ai_hp~Articles$phd,ylim=c(0,6),ylab = "AI",xlab = "phd", main = "hP model: AI versus phd")
abline(1,0,lty=2)
plot(ai_cmp~Articles$phd,ylim=c(0,6),ylab = "AI",xlab = "phd", main = "CMP model: AI versus phd")
abline(1,0,lty=2)
