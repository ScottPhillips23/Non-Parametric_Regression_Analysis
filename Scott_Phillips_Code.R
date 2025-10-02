data(package = .packages(all.available = TRUE))
#Data and preliminary statistics and plots
###########################
install.packages('NHANES')
library(NHANES)
library(pheatmap)
library(car)
library(bayesboot)
library(boot)
library(ggplot2)
library(gridExtra)
data(NHANES)
attach(NHANES)
names(NHANES)

#Checking which variables have enough info to analyse

sum(is.na(NHANES$Poverty))
sum(is.na(NHANES$HHIncomeMid))
sum(is.na(NHANES$BMI))
sum(is.na(NHANES$Pulse))
sum(is.na(NHANES$DaysMentHlthBad))
sum(is.na(NHANES$Testosterone))
sum(is.na(NHANES$SleepHrsNight))
sum(is.na(NHANES$AlcoholYear))
sum(is.na(NHANES$AlcoholDay))
sum(is.na(NHANES$SexNumPartnLife))
sum(is.na(NHANES$Weight))
sum(is.na(NHANES$TotChol))
sum(is.na(NHANES$PhysActiveDays))
sum(is.na(NHANES$SmokeAge))
HealthdataOG <- data.frame(DaysPhysHlthBad, DaysMentHlthBad, Age, Gender, SleepHrsNight,
                          Poverty, BMI, Weight, BPSysAve, BPDiaAve, Pulse, TotChol,
                          AlcoholYear, SmokeAge, Testosterone, SexNumPartnLife)
Healthdata <- na.omit(HealthdataOG)
Healthall <- Healthdata[, -4]
Healthmale <- Healthdata[Healthdata$Gender == 'male',]
Healthmale <- Healthmale[, -4]
Healthfemale <- Healthdata[Healthdata$Gender == 'female',]
Healthfemale <- Healthfemale[, -4]

#First thing we want to do is calculate some mean and variance statistics for our variables

par(mfrow = c(1,1))

Mentmean <- mean(Healthdata$DaysMentHlthBad)
Mentmean
Mentvar <- var(Healthdata$DaysMentHlthBad)
Mentvar

boot.mean <- function(data, indices){
  dt <- data[indices]
  return(mean(dt))
}
boot.var <- function(data, indices){
  dt <- data[indices]
  return(var(dt))
}

Mentmean.bb <- bayesboot(data = Healthdata$DaysMentHlthBad, statistic = mean, R = 10000, 
                      use.weights = FALSE)
Mentvar.bb <- bayesboot(data = Healthdata$DaysMentHlthBad, statistic = var, R = 10000,
                        use.weights = FALSE)
Mentmean.cb <- boot(Healthdata$DaysMentHlthBad, statistic = boot.mean, R = 10000)
Mentvar.cb <- boot(Healthdata$DaysMentHlthBad, statistic = var, R = 10000)

hist(Mentvar.cb$t)
hist(Mentvar.bb$V1)

#Plot we didn't use to showcase bootstrapped means of predictor variable

hist(Mentmean.cb$t, breaks = 20, 
     main = 'Bootstrap Estimate of mean of bad mental health days',
     col = rgb(0, 0, 1, 1/4), xlab = '')
hist(Mentmean.bb$V1, breaks = 20, 
     col = rgb(1, 0, 0, 1/4), xlab = '', add = TRUE)
legend('topright', legend = c('Classical bootstrap', 'Bayesian bootstrap', 'Sample Mean'),
       col = c(rgb(0,0,0,1/4),rgb(1,0,0,1/4), 'red'), pch = 15, bty = 'n')
points(Mentmean,0, col = 'red', pch = 22, bg = 'red')

#
par(mfrow = c(1,1))
hist(log(Healthdata$DaysMentHlthBad+1), breaks = 20, 
     xlab = 'Log of Bad Mental Health Days', 
     main = 'Histogram of log of Bad Mental Health days', probability = TRUE)
lines(density(log(Healthdata$DaysMentHlthBad)))

#Nice looking histograms of predictor variable
plot2.tot <- ggplot(data = data.frame(x = log(Healthdata$DaysMentHlthBad+1)), aes(x = x)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.25, fill = 'lightblue', color = 'black') +
  labs(title = 'Histogram of log of bad mental health days', x = 'Log of days', y = 'Frequency') +
  geom_density(color = 'red', size = 1) 

hist(Healthdata$DaysMentHlthBad, breaks = 20, 
     xlab = 'Bad Mental Health Days',
     main = 'Histogram of bad mental health days', probability = TRUE)
lines(density(Healthdata$DaysMentHlthBad))

plot1.tot <- ggplot(data = data.frame(x = (Healthdata$DaysMentHlthBad)), aes(x = x)) +
  geom_histogram(aes(y = ..density..), binwidth = 2, fill = 'lightblue', color = 'black') +
  labs(title = 'Histogram of bad mental health days', x = 'Days', y = 'Frequency') +
  geom_density(color = 'red', size = 1)

grid.arrange(plot1.tot, plot2.tot,
             layout_matrix = matrix(1:2, 1,2))







##############Preliminary bstatistics we will need
library(e1071)
mean.DNP <- function(data, indices){
  mean(data[indices])
}
var.DNP <- function(data, indices){
  var(data[indices])
}
skew.DNP <- function(data, indices){
  skewness(data[indices])
}

#Female data
library(boot)
Healthfemaledata <- as.numeric(Healthfemale$DaysMentHlthBad)

female.ment.mean.mean <- mean(boot(Healthfemale$DaysMentHlthBad, statistic = mean.DNP, R = 10000)$t)
female.ment.mean.var <- var(boot(Healthfemale$DaysMentHlthBad, statistic = mean.DNP, R = 10000)$t)
female.ment.var.mean <- mean(boot(Healthfemale$DaysMentHlthBad, statistic = var.DNP, R = 10000)$t)
female.ment.var.var <- var(boot(Healthfemale$DaysMentHlthBad, statistic = var.DNP, R = 10000)$t)
female.ment.skew.mean <- mean(boot(Healthfemale$DaysMentHlthBad, statistic = skew.DNP, R = 10000)$t)

female.sleep.mean.mean <- mean(boot(Healthfemale$SleepHrsNight, statistic = mean.DNP, R = 10000)$t)
female.sleep.mean.var <- var(boot(Healthfemale$SleepHrsNight, statistic = mean.DNP, R = 10000)$t)
female.sleep.var.mean <- mean(boot(Healthfemale$SleepHrsNight, statistic = var.DNP, R = 10000)$t)
female.sleep.var.var <- var(boot(Healthfemale$SleepHrsNight, statistic = var.DNP, R = 10000)$t)
female.pov.mean.mean <- mean(boot(Healthfemale$Poverty, statistic = mean.DNP, R = 10000)$t)
female.pov.mean.var <- var(boot(Healthfemale$Poverty, statistic = mean.DNP, R = 10000)$t)
female.pov.var.mean <- mean(boot(Healthfemale$Poverty, statistic = var.DNP, R = 10000)$t)
female.pov.var.var <- var(boot(Healthfemale$Poverty, statistic = var.DNP, R = 10000)$t)
female.pulse.mean.mean <- mean(boot(Healthfemale$Pulse, statistic = mean.DNP, R = 10000)$t)
female.pulse.mean.var <- var(boot(Healthfemale$Pulse, statistic = mean.DNP, R = 10000)$t)
female.pulse.var.mean <- mean(boot(Healthfemale$Pulse, statistic = var.DNP, R = 10000)$t)
female.pulse.var.var <- var(boot(Healthfemale$Pulse, statistic = var.DNP, R = 10000)$t)
female.sex.mean.mean <- mean(boot(Healthfemale$SexNumPartnLife, statistic = mean.DNP, R = 10000)$t)
female.sex.mean.var <- var(boot(Healthfemale$SexNumPartnLife, statistic = mean.DNP, R = 10000)$t)
female.sex.var.mean <- mean(boot(Healthfemale$SexNumPartnLife, statistic = var.DNP, R = 10000)$t)
female.sex.var.var <- var(boot(Healthfemale$SexNumPartnLife, statistic = var.DNP, R = 10000)$t)

#Male data

male.ment.mean.mean <- mean(boot(Healthmale$DaysMentHlthBad, mean.DNP, R = 10000)$t)
male.ment.mean.var <- var(boot(Healthmale$DaysMentHlthBad, mean.DNP, R = 10000)$t)
male.ment.var.mean <- mean(boot(Healthmale$DaysMentHlthBad, var.DNP, R = 10000)$t)
male.ment.var.var <- var(boot(Healthmale$DaysMentHlthBad, var.DNP, R = 10000)$t)
male.ment.skew.mean <- mean(boot(Healthmale$DaysMentHlthBad, skew.DNP, R = 10000)$t)

male.phys.mean.mean <- mean(boot(Healthmale$DaysPhysHlthBad, mean.DNP, R = 10000)$t)
male.phys.mean.var <- var(boot(Healthmale$DaysPhysHlthBad, mean.DNP, R = 10000)$t)
male.phys.var.mean <- mean(boot(Healthmale$DaysPhysHlthBad, var.DNP, R = 10000)$t)
male.phys.var.var <- var(boot(Healthmale$DaysPhysHlthBad, var.DNP, R = 10000)$t)
male.sleep.mean.mean <- mean(boot(Healthmale$SleepHrsNight, mean.DNP, R = 10000)$t)
male.sleep.mean.var <- var(boot(Healthmale$SleepHrsNight, mean.DNP, R = 10000)$t)
male.sleep.var.mean <- mean(boot(Healthmale$SleepHrsNight, var.DNP, R = 10000)$t)
male.sleep.var.var <- var(boot(Healthmale$SleepHrsNight, var.DNP, R = 10000)$t)
male.test.mean.mean <- mean(boot(Healthmale$Testosterone, mean.DNP, R = 10000)$t)
male.test.mean.var <- var(boot(Healthmale$Testosterone, mean.DNP, R = 10000)$t)
male.test.var.mean <- mean(boot(Healthmale$Testosterone, var.DNP, R = 10000)$t)
male.test.var.var <- var(boot(Healthmale$Testosterone, var.DNP, R = 10000)$t)


#Now lets make a table, one for men and one for women, rows correspond 
#to mean and variances of variables, columns are means and variances of the measn and variances

library(ggpubr)

# Create a data frame for men
male_table <- data.frame(
  Male_Variable = c("DaysMentHlthBad", "DaysPhysHlthBad", "SleepHrsNight", "Testosterone"),
  Mean_of_Means = c(male.ment.mean.mean, male.phys.mean.mean, male.sleep.mean.mean, male.test.mean.mean),
  Variance_of_Means = c(male.ment.mean.var, male.phys.mean.var, male.sleep.mean.var, male.test.mean.var),
  Mean_of_Variances = c(male.ment.var.mean, male.phys.var.mean, male.sleep.var.mean, male.test.var.mean),
  Variance_of_Variances = c(male.ment.var.var, male.phys.var.var, male.sleep.var.var, male.test.var.var)
)


# View the table for men
print(male_table)

# Create a data frame for the female dataset
female_table <- data.frame(
  Female_Variable = c("DaysMentHlthBad", "SleepHrsNight", "Poverty", "Pulse", "SexNumPartnLife"),
  Mean_of_Means = c(female.ment.mean.mean,female.sleep.mean.mean,female.pov.mean.mean,
                    female.pulse.mean.mean, female.sex.mean.mean),
  Variance_of_Means = c(female.ment.mean.var,female.sleep.mean.var,female.pov.mean.var,
                        female.pulse.mean.var,female.sex.mean.var),
  Mean_of_Variances = c(female.ment.var.mean,female.sleep.var.mean,female.pov.var.mean,
                        female.pulse.var.mean,female.sex.var.mean),
  Variance_of_Variances = c(female.ment.var.var,female.sleep.var.var,female.pov.var.var,
                            female.pulse.var.var,female.sex.var.var)
)

female_table$Gender <- "Female"
male_table$Gender <- "Male"

colnames(female_table)[1] <- "Variable"
colnames(male_table)[1] <- "Variable"

combined_table <- rbind(female_table, male_table)

table_plot <- ggtexttable(combined_table, rows = NULL)

table_plot
###########################
#Linear Regression male data
###########################

#Initial linear model
linear.model3 <- lm(log(DaysMentHlthBad+1)~ ., data = Healthmale)
summary(linear.model3)

Healthmale

#Scaling of linear models
apply(Healthmale[, -2], 2, range)
scaled.Healthmale <- Healthmale
scaled.Healthmale[, -2] <- scale(Healthmale[, -2])
linear.model4 <- lm(log(DaysMentHlthBad+1) ~., data = scaled.Healthmale)
summary(linear.model4)

#Linear model after accounting for multicollinearity
covHealthmale <- cov(scaled.Healthmale[, -2])
corHealthmale <- cov2cor(covHealthmale)
pheatmap(corHealthmale, display_numbers = TRUE, cluster_Rows = FALSE,
         cluster_cols = FALSE, fontsize_numbers = 10,
         main = 'Correlations between male predictors')
linear.model4 <- lm(log(DaysMentHlthBad+1) ~. -BMI, data = scaled.Healthmale)
vif(linear.model4)
summary(linear.model4)

#We work with the scale data to give the C.I.'s comparability later on 

#Standard Confidence intervals
Wald.male <- confint(linear.model4, level = 0.99)
Wald.male

#Bootstrap using NPB
set.seed(1)
NPB.male <- Boot(linear.model4, method = 'case', R = 10000)
summary(NPB.male)

normal.NPB.male <- confint(NPB.male, level = 0.99, type = 'norm')
basic.NPB.male <- confint(NPB.male, level = 0.99, type = 'basic')
percent.NPB.male <- confint(NPB.male, level = 0.99, type = 'perc')
bc.NPB.male <- confint(NPB.male, level = 0.99, type = 'bca')

#Bootstrap using SRB now
set.seed(1)
SRB.male <- Boot(linear.model4, method = 'residual', R = 10000)

normal.SRB.male <- confint(SRB.male, level = 0.99, type = 'norm')
basic.SRB.male <- confint(SRB.male, level = 0.99, type = 'basic')
percent.SRB.male <- confint(SRB.male, level = 0.99, type = 'perc')
bc.SRB.male <- confint(SRB.male, level = 0.99, type = 'bca')

#Finally using the Bayesian bootstrap
betas.male <- function(data){
  coef(lm(log(DaysMentHlthBad + 1) ~. - BMI, data = data))
}

set.seed(1)
BBR.male <- bayesboot(data = scaled.Healthmale, statistic = betas.male, R = 10000, 
                      use.weights = FALSE, .progress = 'text')
summary(BBR.male)
conf.boot.male <- apply(BBR.male, 2, function(x) quantile(x, probs = c(0.05, 0.995)))
mean.boot.male <- apply(BBR.male, 2, function(x) mean(x))
mean.boot.male[2]
conf.boot.male
bayes.male.Phys.Hlth <- conf.boot.male[,2]
bayes.male.Sleep <- conf.boot.male[,4]
bayes.male.Test <- conf.boot.male[,13]

############C.I.'s for Phys.Hlth

CI.male.Phys.Hlth <- rbind(Wald.male[2,], normal.NPB.male[2,], basic.NPB.male[2,], 
                           percent.NPB.male[2,], bc.NPB.male[2,], normal.SRB.male[2,],
                           basic.SRB.male[2,], percent.SRB.male[2,], bc.SRB.male[2,], 
                           conf.boot.male[,2])
rownames(CI.male.Phys.Hlth) <- c('Wald.PhysHlth', 'normal.NPB.PhysHlth', 'basic.NPB.PhysHlth', 
                                 'percent.NPB.PhysHlth', 'bc.NPB.PhysHlth','normal.SRB.PhysHlth', 
                                 'basic.SRB.PhysHlth', 'percent.SRB.PhysHlth', 
                                 'bc.SRB.PhysHlth','bayes.PhysHlth')

point.male.Phys.Hlth <- c(rep(coef(linear.model4)[2], 9), mean.boot.male[2])

library(plotrix)

plotCI(x = 1:10, y = point.male.Phys.Hlth, li = CI.male.Phys.Hlth[,1], ui = CI.male.Phys.Hlth[,2],
       lwd = 2, main = 'DaysBadPhysHlth', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n', ylim = c(0.1,0.5))
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)


###############C.I.'s for sleep

CI.male.Sleep <- rbind(Wald.male[4,], normal.NPB.male[4,], basic.NPB.male[4,], 
                       percent.NPB.male[4,], bc.NPB.male[4,], normal.SRB.male[4,],
                       basic.SRB.male[4,], percent.SRB.male[4,], bc.SRB.male[4,], 
                       conf.boot.male[,4])
rownames(CI.male.Sleep) <- c('Wald.Sleep', 'normal.NPB.Sleep', 'basic.NPB.Sleep', 
                                 'percent.NPB.Sleep', 'bc.NPB.Sleep','normal.SRB.Sleep', 
                                 'basic.SRB.Sleep', 'percent.SRB.Sleep', 
                                 'bc.SRB.Sleep','bayes.Sleep')

point.male.Sleep <- c(rep(coef(linear.model4)[4], 9), mean.boot.male[4])

library(plotrix)

plotCI(x = 1:10, y = point.male.Sleep, li = CI.male.Sleep[,1], ui = CI.male.Sleep[,2],
       lwd = 2, main = 'Hours of sleep', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#######################Now for Testosterone

CI.male.Test <- rbind(Wald.male[13,], normal.NPB.male[13,], basic.NPB.male[13,], 
                       percent.NPB.male[13,], bc.NPB.male[13,], normal.SRB.male[13,],
                       basic.SRB.male[13,], percent.SRB.male[13,], bc.SRB.male[13,], 
                       conf.boot.male[,13])
rownames(CI.male.Test) <- c('Wald.Test', 'normal.NPB.Test', 'basic.NPB.Test', 
                             'percent.NPB.Test', 'bc.NPB.Test','normal.SRB.Test', 
                             'basic.SRB.Test', 'percent.SRB.Test', 
                             'bc.SRB.Test','bayes.Test')

point.male.Test <- c(rep(coef(linear.model4)[13], 9), mean.boot.male[13])

library(plotrix)

plotCI(x = 1:10, y = point.male.Test, li = CI.male.Test[,1], ui = CI.male.Test[,2],
       lwd = 2, main = 'Testosterone', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n', ylim = c(-0.3, 0.00001))
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#Combined C.I.'s for data

par(mfrow = c(1,1))
male.statsig.point <- c(rep(coef(linear.model4)[2], 9), mean.boot.male[2],
                        rep(coef(linear.model4)[4], 9), mean.boot.male[4],
                        rep(coef(linear.model4)[13], 9), mean.boot.male[13])
male.statsig.CI <- rbind(CI.male.Phys.Hlth, CI.male.Sleep, CI.male.Test)

plotCI(x = 1:30, y = male.statsig.point, li = male.statsig.CI[,1], ui = male.statsig.CI[,2],
       lwd = 2, main = 'Male C.I.', 
       col = c(1,rep(2,4),rep(3,4),4,5,rep(6,4),rep(7,4),8,9,rep(10,4),rep(11,4),12), 
       bty = 'l', ylab = '', xlab = 'Methods', xaxt ='n')
axis(1, at = 1:30, labels = rep(c('Wald', 'Normal', 'Basic', 'Percent', 'BC',
                              'Normal', 'Basic', 'Percent', 'BC','Bayes'), 3), cex.axis = 0.8,
     las = 2)
legend('topright', legend = c('WaldPhys','NPBPhys', 'SRBPhys','BayesPhys',
                              'WaldSleep','NPBSleep', 'SRBSleep','BayesSleep',
                              'WaldTest','NPBTest', 'SRBTest','BayesTest'), 
       col = seq(1:15), pch = rep(15,12), bty = 'n', cex = 0.6, ncol = 2)
abline(h = 0, lty = 3)

#So we now have all the confidence intervals for all statistically significant results at the 99% level
#And can proceed to the next stage

###########################
#Linear Regression female data
###########################

#Initial linear model
linear.model5 <- lm(log(DaysMentHlthBad + 1) ~., data = Healthfemale)
summary(linear.model5)

#Scaled linear model
apply(Healthfemale[, -2], 2, range)
scaled.Healthfemale <- Healthfemale
scaled.Healthfemale[, -2] <- scale(Healthfemale[, -2])
linear.model6 <- lm(log(DaysMentHlthBad+1) ~., data = scaled.Healthfemale)
summary(linear.model6)

#Final linear model after accounting for multicollinearity
covHealthfemale <- cov(scaled.Healthfemale[, -2])
corHealthfemale <- cov2cor(covHealthfemale)
pheatmap(corHealthfemale, display_numbers = TRUE, cluster_Rows = FALSE,
         cluster_cols = FALSE, fontsize_numbers = 10,
         main = 'Correlations between female predictors')
linear.model6 <- lm(log(DaysMentHlthBad + 1) ~. -Weight, data = scaled.Healthfemale)
vif(linear.model6)
summary(linear.model6)

#Wald interval first

Wald.female <- confint(linear.model6, level = 0.99)

#Now first find confidence intervals for NPB
set.seed(2)
NPB.female <- Boot(linear.model6, method = 'case', R = 10000)
summary(NPB.female)

normal.NPB.female = confint(NPB.female, level = 0.99, type = "norm")
basic.NPB.female = confint(NPB.female, level = 0.99, type = "basic")
percent.NPB.female = confint(NPB.female, level = 0.99, type = "perc")
bc.NPB.female = confint(NPB.female, level = 0.99, type = "bca")


#SRB

set.seed(1)
SRB.female <- Boot(linear.model6, method = 'residual', R = 10000)

normal.SRB.female = confint(SRB.female, level = 0.99, type = "norm")
basic.SRB.female = confint(SRB.female, level = 0.99, type = "basic")
percent.SRB.female = confint(SRB.female, level = 0.99, type = "perc")
bc.SRB.female = confint(SRB.female, level = 0.99, type = "bca")

basic.SRB.female

basic.NPB.female

#Bayesian 

betas.female <- function(data){
  coef(lm(log(DaysMentHlthBad + 1) ~. - Weight, data = data))
}

set.seed(1)
BBR.female <- bayesboot(data = scaled.Healthfemale, statistic = betas.female, R = 10000, 
                        use.weights = FALSE, .progress = 'text')

conf.boot.female <- apply(BBR.female, 2, function(x) quantile(x, probs = c(0.005, 0.995)))
mean.boot.female <- apply(BBR.female, 2, function(x) mean(x))

#Starting with CI's for sleep

CI.female.sleep <- rbind(Wald.female[4,], normal.NPB.female[4,], basic.NPB.female[4,],
                         percent.NPB.female[4,], bc.NPB.female[4,], normal.SRB.female[4,],
                         basic.SRB.female[4,], percent.SRB.female[4,], bc.SRB.female[4,], 
                         conf.boot.female[,4])

rownames(CI.female.sleep) <- c('Wald.sleep', 'normal.NPB.sleep', 'basic.NPB.sleep', 
                               'percent.NPB.sleep', 'bc.NPB.sleep','normal.SRB.sleep', 
                               'basic.SRB.sleep', 'percent.SRB.sleep', 
                               'bc.SRB.sleep','bayes.sleep')

point.female.sleep <- c(rep(coef(linear.model6)[4], 9), mean.boot.female[4])

plotCI(x = 1:10, y = point.female.sleep,li = CI.female.sleep[,1], ui = CI.female.sleep[,2], 
       lwd = 2, main = 'Sleep', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#C.I.'s for poverty

CI.female.poverty <- rbind(Wald.female[5,], normal.NPB.female[5,], basic.NPB.female[5,],
                         percent.NPB.female[5,], bc.NPB.female[5,], normal.SRB.female[5,],
                         basic.SRB.female[5,], percent.SRB.female[5,], bc.SRB.female[5,], 
                         conf.boot.female[,5])

rownames(CI.female.poverty) <- c('Wald.poverty', 'normal.NPB.poverty', 'basic.NPB.poverty', 
                               'percent.NPB.poverty', 'bc.NPB.poverty','normal.SRB.poverty', 
                               'basic.SRB.poverty', 'percent.SRB.poverty', 
                               'bc.SRB.poverty','bayes.poverty')

point.female.poverty <- c(rep(coef(linear.model6)[5], 9), mean.boot.female[5])

plotCI(x = 1:10, y = point.female.poverty,li = CI.female.poverty[,1], ui = CI.female.poverty[,2], 
       lwd = 2, main = 'poverty', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#Pulse

CI.female.pulse <- rbind(Wald.female[9,], normal.NPB.female[9,], basic.NPB.female[9,],
                           percent.NPB.female[9,], bc.NPB.female[9,], normal.SRB.female[9,],
                           basic.SRB.female[9,], percent.SRB.female[9,], bc.SRB.female[9,], 
                           conf.boot.female[,9])

rownames(CI.female.pulse) <- c('Wald.pulse', 'normal.NPB.pulse', 'basic.NPB.pulse', 
                                 'percent.NPB.pulse', 'bc.NPB.pulse','normal.SRB.pulse', 
                                 'basic.SRB.pulse', 'percent.SRB.pulse', 
                                 'bc.SRB.pulse','bayes.pulse')

point.female.pulse <- c(rep(coef(linear.model6)[9], 9), mean.boot.female[9])

plotCI(x = 1:10, y = point.female.pulse,li = CI.female.pulse[,1], ui = CI.female.pulse[,2], 
       lwd = 2, main = 'pulse', col = c(1,rep(2,4),rep(3,4),4), bty = 'l', ylab = '',
       xlab = 'Methods', xaxt ='n', ylim = c(-0.0001, 0.35))
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#Number of Sexual Partners

CI.female.partn <- rbind(Wald.female[14,], normal.NPB.female[14,], basic.NPB.female[14,],
                         percent.NPB.female[14,], bc.NPB.female[14,], normal.SRB.female[14,],
                         basic.SRB.female[14,], percent.SRB.female[14,], bc.SRB.female[14,], 
                         conf.boot.female[,14])

rownames(CI.female.partn) <- c('Wald.partn', 'normal.NPB.partn', 'basic.NPB.partn', 
                               'percent.NPB.partn', 'bc.NPB.partn','normal.SRB.partn', 
                               'basic.SRB.partn', 'percent.SRB.partn', 
                               'bc.SRB.partn','bayes.partn')

point.female.partn <- c(rep(coef(linear.model6)[14], 9), mean.boot.female[14])

plotCI(x = 1:10, y = point.female.partn,li = CI.female.partn[,1], ui = CI.female.partn[,2], 
       lwd = 2, main = 'number of sexual partners', col = c(1,rep(2,4),rep(3,4),4), bty = 'l',
       ylab = '', xlab = 'Methods', xaxt ='n')
axis(1, at = 1:10, labels = c('Wald', 'normal', 'basic', 'percent', 'BC',
                              'normal', 'basic', 'percent', 'BC','Bayes'), cex.axis = 0.8)
legend('topright', legend = c('NPB','SRB'), col = c(2,3), pch = c(15,15), bty = 'n', cex = 0.8)
abline(h = 0, lty = 3)

#C.I.'s with all information

female.statsig.point <- c(rep(coef(linear.model6)[4], 9), mean.boot.female[4],
                          rep(coef(linear.model6)[5], 9), mean.boot.female[5],
                          rep(coef(linear.model6)[9], 9), mean.boot.female[9],
                          rep(coef(linear.model6)[14], 9), mean.boot.female[14])

female.statsig.CI <- rbind(CI.female.sleep, CI.female.poverty, CI.female.pulse, CI.female.partn)
female.statsig.point
female.statsig.CI  

plotCI(x = 1:40, y = female.statsig.point, li = female.statsig.CI[,1], ui = female.statsig.CI[,2],
       lwd = 2, main = 'Female C.I.', 
       col = c(1,rep(2,4),rep(3,4),4,5,rep(6,4),rep(7,4),8,9,rep(10,4),rep(11,4),12, 13, rep(14,4),
               rep(15,4), 16), 
       bty = 'l', ylab = '', xlab = 'Methods', xaxt ='n')
axis(1, at = 1:40, labels = rep(c('Wald', 'Normal', 'Basic', 'Percent', 'BC',
                                  'Normal', 'Basic', 'Percent', 'BC','Bayes'), 4), cex.axis = 0.8,
     las = 2)
legend('topleft', legend = c('WaldSleep','NPBSleep', 'SRBSleep','BayesSleep',
                             'WaldPov','NPBPov', 'SRBPov','BayesPov',
                             'WaldPulse','NPBPulse', 'SRBPulse','BayesPulse',
                             'WaldPartn', 'NPBPartn', 'SRBPartn', 'BayesPartn'), 
       col = seq(1:16), pch = rep(15,16), bty = 'n', cex = 0.6, ncol = 2)
abline(h = 0, lty = 3)

###########################
#Clustering Male
###########################
#Now applying to see if we have clusters for our multivariate male data

library(BNPmix)
library(gridExtra)

#Initialise male datasets and grids for evaluation of model
Ymale = (Healthmale[, c('DaysMentHlthBad', 'DaysPhysHlthBad', 'SleepHrsNight', 'Testosterone')])

n = nrow(Ymale)

grid.length = 20
grid1 <- seq(min(Ymale[,1]), max(Ymale[,1]), length.out = grid.length)
grid2 <- seq(min(Ymale[,2]), max(Ymale[,2]), length.out = grid.length)
grid3 <- seq(min(Ymale[,3]), max(Ymale[,3]), length.out = grid.length)
grid4 <- seq(min(Ymale[,4]), max(Ymale[,4]), length.out = grid.length)

grid <- expand.grid(grid1, grid2, grid3, grid4)
dim(grid)
pairs(Ymale)

#Running of Clustering model
library(BNPmix)
DPprior <- PYcalibrate(Ek = 3, n = n, discount = 0)
mcmc <- list(niter = 2000, nburn = 1000, method = 'SLI')
prior <- list(strength = DPprior$strength, discount = 0, hyper = TRUE,
              m0 = c(4.213, 3.587, 6.493, 419.1), k0 = 0.01,
              Sigma0 = diag(c(69.8, 57.2, 1.60, 1000)))
output <- list(grid = grid, out_param = TRUE)
set.seed(1)
mult.fit.male <- PYdensity(y = Ymale, mcmc = mcmc, prior = prior, output = output)
summary(mult.fit.male)
apply()

plot(mult.fit.male, dim = c(4,1), show_clust = TRUE)

#Produce cluster plots for densities of outputs

means <- mult.fit.male$mean

post.means1 = matrix(NA, nrow = 1000, ncol = 3)
post.means2 = matrix(NA, nrow = 1000, ncol = 3)
post.means3 = matrix(NA, nrow = 1000, ncol = 3)
post.means4 = matrix(NA, nrow = 1000, ncol = 3)

for(i in 1:1000){
  post.means1[i,] = means[[i]][,1][1:3]
  post.means2[i,] = means[[i]][,2][1:3]
  post.means3[i,] = means[[i]][,2][1:3]
  post.means4[i,] = means[[i]][,4][1:3]
}


all.means1.male <- data.frame(cluster = factor(rep(c('Cluster 1', 'Cluster2', 'Cluster3'),
                                                  each = 1000)),
                             samples = c(post.means1[,1], post.means1[,2], post.means1[,3]))
all.means2.male <- data.frame(cluster = factor(rep(c('Cluster 1', 'Cluster2', 'Cluster3'),
                                                   each = 1000)),
                              samples = c(post.means2[,1], post.means2[,2], post.means2[,3]))
all.means3.male <- data.frame(cluster = factor(rep(c('Cluster 1', 'Cluster2', 'Cluster3'),
                                                   each = 1000)),
                              samples = c(post.means3[,1], post.means3[,2], post.means3[,3]))
all.means4.male <- data.frame(cluster = factor(rep(c('Cluster 1', 'Cluster2', 'Cluster3'),
                                                   each = 1000)),
                              samples = c(post.means4[,1], post.means4[,2], post.means4[,3]))


library(ggplot2)

plot.m1.male <- ggplot(all.means1.male, aes(x = samples, color = cluster))+
  geom_density(linewidth = 1.3)+
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Days') + ylab('Posterior density') + 
  ggtitle('Bad Mental Health')+
  xlim(0,30)

plot.m1.male

plot.m2.male <- ggplot(all.means2.male, aes(x = samples, color = cluster))+
  geom_density(linewidth = 1.3)+
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Days') + ylab('Posterior density') + 
  ggtitle('Bad Physical Health')+
  xlim(0,30)

plot.m2.male

plot.m3.male <- ggplot(all.means3.male, aes(x = samples, color = cluster))+
  geom_density(linewidth = 1.3)+
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Hours') + ylab('Posterior density') + 
  ggtitle('Sleep')+
  xlim(0,10)

plot.m3.male

plot.m4.male <- ggplot(all.means4.male, aes(x = samples, color = cluster))+
  geom_density(linewidth = 1.3)+
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('ng/dL') + ylab('Posterior density') + 
  ggtitle('Testosterone')+
  xlim(0,1000)

plot.m4.male

#Remove legends from most plots
plot.m1.male <- plot.m1.male + theme(legend.position = "none")
plot.m2.male <- plot.m2.male + theme(legend.position = "none")
plot.m3.male <- plot.m3.male + theme(legend.position = "none")

grid.arrange(plot.m1.male, plot.m2.male, plot.m3.male, plot.m4.male,
             layout_matrix = matrix(1:4, 1, 4), widths = c(1,1,1,1.5))

#Now we look at the posterior summaries of how the data points are assigned to the clusters.
#Do not include in report

pi.male = mult.fit.male$probs
pi.male[[1]][3]

post.pi.male <- matrix(NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  post.pi.male[i,1] = pi.male[[i]][1]
  post.pi.male[i,2] = pi.male[[i]][2]
  post.pi.male[i,3] = pi.male[[i]][3]
}
post.pi.male <- na.omit(post.pi.male)
summary(post.pi)

par(mfrow = c(1,3))
plot(density(post.pi.male[,1]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[1]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.pi.male[,2]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[2]), bty = 'l', col = 2, cex.lab = 1.3)
plot(density(post.pi.male[,3]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[3]), bty = 'l', col = 2, cex.lab = 1.3)
mtext('3-component mixture weights', side = 3, line = -2, outer = TRUE,
      cex = 1.5)


###########################
#Clustering Female
###########################
#Now applying to see if we have clusters for our multivariate female data, first insitalising
#the data sets and constructing the grids for evaluation
Yfemale = Healthfemale[, c('DaysMentHlthBad', 'SleepHrsNight', 'Poverty', 'Pulse', 'SexNumPartnLife')]
n = nrow(Yfemale)
summary(Y)

grid.length = 10
grid1 <- seq(min((Yfemale[,1])), max((Yfemale[,1])), length.out = grid.length)
grid2 <- seq(min(Yfemale[,2]), max(Yfemale[,2]), length.out = grid.length)
grid3 <- seq(min(Yfemale[,3]), max(Yfemale[,3]), length.out = grid.length)
grid4 <- seq(min(Yfemale[,4]), max(Yfemale[,4]), length.out = grid.length)
grid5 <- seq(min(Yfemale[,5]), max(Yfemale[,5]), length.out = grid.length)

grid <- expand.grid(grid1, grid2, grid3, grid4, grid5)
dim(grid)

#Running clustering algorithm
DPprior <- PYcalibrate(Ek = 2, n = n, discount = 0)
mcmc <- list(niter = 2000, nburn = 1000)
prior <- list(strength = DPprior$strength, discount = 0, hyper = FALSE,
              m0 =c(4.770, 6.885, 3.001, 74.78, 14.42), k0 = 0.01, Sigma0 = diag(c(56.64, 1.766, 2.811, 123.5, 1293)))
output <- list(grid = grid, out_param = TRUE)
set.seed(1)
mult.fit.female <- PYdensity(y = Yfemale, mcmc = mcmc, prior = prior, output = output)
summary(mult.fit.female)

#Producing posterior densities for the clusters 

means <- mult.fit.female$mean

post.means1 = matrix(NA, nrow = 1000, ncol = 3)
post.means2 = matrix(NA, nrow = 1000, ncol = 3)
post.means3 = matrix(NA, nrow = 1000, ncol = 3)
post.means4 = matrix(NA, nrow = 1000, ncol = 3)
post.means5 = matrix(NA, nrow = 1000, ncol = 3)
for(i in 1:(1000)){
  post.means1[i,] = means[[i]][,1][1:3]
  post.means2[i,] = means[[i]][,2][1:3]
  post.means3[i,] = means[[i]][,2][1:3]
  post.means4[i,] = means[[i]][,4][1:3]
  post.means5[i,] = means[[i]][,5][1:3]
}
summary(post.means1)


library(ggplot2)
all.means1.female = data.frame(cluster = factor(rep(c('Cluster1','Cluster2','Cluster3'), each = 1000)), 
                        samples = c(post.means1[,1],post.means1[,2],post.means1[,3]))
all.means2.female = data.frame(cluster = factor(rep(c('Cluster1','Cluster2','Cluster3'), each = 1000)), 
                        samples = c(post.means2[,1],post.means2[,2],post.means1[,3]))
all.means3.female = data.frame(cluster = factor(rep(c('Cluster1','Cluster2','Cluster3'), each = 1000)), 
                        samples = c(post.means3[,1],post.means3[,2],post.means1[,3]))
all.means4.female = data.frame(cluster = factor(rep(c('Cluster1','Cluster2','Cluster3'), each = 1000)), 
                        samples = c(post.means4[,1],post.means4[,2],post.means1[,3]))
all.means5.female = data.frame(cluster = factor(rep(c('Cluster1','Cluster2','Cluster3'), each = 1000)), 
                        samples = c(post.means5[,1],post.means5[,2],post.means1[,3]))

plot.m1.female = ggplot(all.means1.female, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Days') + ylab('Posterior density') + xlim(0,30) +
  ggtitle('Bad Mental Health') 

plot.m2.female = ggplot(all.means2.female, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Hours') + ylab('Posterior density') + xlim(0,20) +
  ggtitle('Sleep') 

plot.m3.female = ggplot(all.means3.female, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Index') + ylab('Posterior density') + xlim(0,30) +
  ggtitle('Poverty') 

plot.m4.female = ggplot(all.means4.female, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('BPM') + ylab('Posterior density') + xlim(60,100) +
  ggtitle('Pulse') 

plot.m5.female = ggplot(all.means5.female, aes(x=samples, color = cluster)) +
  geom_density(size = 1.3) +
  theme(plot.title = element_text(size=18), axis.title=element_text(size=15), axis.text=element_text(size=10)) +
  xlab('Number') + ylab('Posterior density') + xlim(0, 50) +
  ggtitle('Sexual Partners') 

# Remove individual legends from plots
plot.m1.female <- plot.m1.female + theme(legend.position = "none")
plot.m2.female <- plot.m2.female + theme(legend.position = "none")
plot.m3.female <- plot.m3.female + theme(legend.position = "none")
plot.m4.female <- plot.m4.female + theme(legend.position = "none")

grid.arrange(plot.m1.female, plot.m2.female, plot.m3.female, plot.m4.female, plot.m5.female,
             layout_matrix = matrix(1:5, 1, 5), widths = c(1,1,1,1,1.4))


#Now we produce posterior plots of how data points are assigned to the clusters.
#Do not include in report

pi.female <- mult.fit.female$probs
post.pi.female <- matrix(NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  post.pi.female[i,1] <- pi.female[[i]][1]
  post.pi.female[i,2] <- pi.female[[i]][2]
  post.pi.female[i,3] <- pi.female[[i]][3]

}
post.pi.female <- na.omit(post.pi.female)
summary(post.pi.female)

par(mfrow = c(1,3))
plot(density(post.pi.female[,1]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[1]), bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.pi.female[,2]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[2]), bty = 'l', col = 2, cex.lab = 1.3)
plot(density(post.pi.female[,3]), main = '', lwd = 2, ylab = 'Posterior density',
     xlab = expression(pi[3]), bty = 'l', col = 1, cex.lab = 1.3)
mtext('3-component mixture weights', side = 3, line = -2, outer = TRUE,
      cex = 1.5)
###########################
#DPM Male
###########################
#Now we want to simply apply a DPM regression model to our data
#Implementation of regression algorithm

#Initialisation of data and grid sequences for evaluation of model

grid_y = seq(min(Healthmale$DaysMentHlthBad), 
             max(Healthmale$DaysMentHlthBad), 
             length.out = 100)

X = Healthmale[, c('DaysPhysHlthBad', 'SleepHrsNight', 'Testosterone')]

n = nrow(X)

grid_x1 <- quantile(X[,1], probs = c(0.5))
grid_x2 <- quantile(X[,2], probs = c(0.01, 0.5, 0.9))
grid_x3 <- quantile(X[,3], probs = c(0.1, 0.5, 0.9))

grid_X <- expand.grid(grid_x1, grid_x2, grid_x3)

#Running of regression model
DPprior <- PYcalibrate(Ek = 2, n = n, discount = 0)
prior <- list(strength = DPprior$strength, discount = 0, hyper = FALSE,
              m0 = c(4.213, 3.587, 6.493, 419.1), s0 = diag(c(15, 10, 1, 100)), a0 = 2, b0 = 1.5)
mcmc <- list(niter = 11000, nburn = 1000, method = 'SLI')
output <- list(grid_x = grid_X, grid_y = grid_y, out_param = TRUE, out_type = 'FULL')
set.seed(1)
fit.male <- PYregression(y = Healthmale$DaysMentHlthBad, x = X, prior = prior, mcmc = mcmc,
                         output = output)

summary(fit.male)

#Produce density plots for the associated quantiles we have evaluated with 
#our regression algorithm

#10% of sleep hours per night
regplot1 <- data.frame(
  dens = as.vector(apply(fit.male$density[, 1:3,], c(1,2), mean)),
  qlow = as.vector(apply(fit.male$density[, 1:3,], c(1,2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.male$density[, 1:3,], c(1,2), quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste("Testosterone = ", grid_x3[1:3]), each = length(grid_y)),
                 level = rep(paste("Testosterone = ", grid_x3)))
)

#50% of sleep hours per night
regplot2 <- data.frame(
  dens = as.vector(apply(fit.male$density[, 4:6,], c(1,2), mean)),
  qlow = as.vector(apply(fit.male$density[, 4:6,], c(1,2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.male$density[, 4:6,], c(1,2), quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste('Testosterone = ', grid_x3[1:3]), each = length(grid_y)),
                 level = rep(paste("Testosterone = ", grid_x3)))
)

#90% of sleep hours per night
regplot3 <- data.frame(
  dens = as.vector(apply(fit.male$density[, 7:9,], c(1,2), mean)),
  qlow = as.vector(apply(fit.male$density[, 7:9,], c(1,2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.male$density[, 7:9,], c(1,2), quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x3)),
  label = factor(rep(paste('Testosterone = ', grid_x3[1:3]), each = length(grid_y)),
                 level = rep(paste('Testosterone = ', grid_x3)))
)

library(ggplot2)
plot1 = ggplot(regplot1) + theme_bw() +
  geom_line(data = regplot1, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot1, map = aes(x = grid, ymin = qlow, ymax = qupp),
              fill = 'green', alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = 'Bad Mental Health Days', y = '3 hours of sleep')

plot2 = ggplot(regplot2) + theme_bw() +
  geom_line(data = regplot2, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot2, map = aes(x = grid, ymin = qlow, ymax = qupp),
              fill = 'red', alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = 'Bad Mental Health Days', y = '6 hours of sleep')

plot3 = ggplot(regplot3) + theme_bw() +
  geom_line(data = regplot3, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot3, map = aes(x = grid, ymin = qlow, ymax = qupp),
              fill = 'blue', alpha = 0.3) +
  facet_wrap(~label, ncol = 5, nrow = 1) +
  labs(x = 'Bad Mental Health Days', y = '8 hours of sleep')

grid.arrange(plot1, plot2, plot3, layout_matrix = matrix(1:3))


#Posterior Regression co-efficient summaries

betas <- fit.male$beta
post.male.betas1 <- matrix(NA, nrow = 10000, ncol = 4)
for(i in 1:10000){
  post.male.betas1[i,] <- betas[[i]][1,]
}
post.male.betas1

#Ugly plots

par(mfrow = c(2,2))
plot(density(post.male.betas1[,1]), main = '', lwd = 2, ylab = '', xlab = expression(beta[0]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.male.betas1[,2]), main = '', lwd = 2, ylab = '', xlab = expression(beta[1]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.male.betas1[,3]), main = '', lwd = 2, ylab = '', xlab = expression(beta[2]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.male.betas1[,4]), main = '', lwd = 2, ylab = '', xlab = expression(beta[3]),
     bty = 'l', col = 1, cex.lab = 1.3)
mtext('Posterior Regression Coefficients', side = 3, line = -2, outer = TRUE,
      cex = 1.5)

#Good plots

plot.male1 <- ggplot(as.data.frame(post.male.betas1[,1]), aes(x = post.male.betas1[, 1]))+
  geom_density(color = 'red', size = 1) +
  labs(x = expression(Beta[0]))
plot.male2 <- ggplot(as.data.frame(post.male.betas1[,2]), aes(x = post.male.betas1[, 2]))+
  geom_density(color = 'red', size = 1)+
  labs(x = expression(Beta[1]))
plot.male3 <- ggplot(as.data.frame(post.male.betas1[,3]), aes(x = post.male.betas1[, 3]))+
  geom_density(color = 'red', size = 1) +
  labs(x = expression(Beta[2]))
plot.male4 <- ggplot(as.data.frame(post.male.betas1[,4]), aes(x = post.male.betas1[, 4]))+
  geom_density(color = 'red', size = 1)+
  labs(x = expression(Beta[3]))
grid.arrange(plot.male1, plot.male2, plot.male3, plot.male4,
             layout_matrix = matrix(1:4, 2,2),
             top = 'Denity plots of regression roefficients for male')

###########################
#DPM Female
###########################
#Now repeating this for the female data

#Initalise data frame for evaluation
Yfemale = Healthfemale[, c('DaysMentHlthBad', 'SleepHrsNight', 'Poverty', 'Pulse', 'SexNumPartnLife')]
n = nrow(Yfemale)


grid_y = seq(min(Healthfemale$DaysMentHlthBad), 
             max(Healthfemale$DaysMentHlthBad), 
             length.out = 100)

X = Healthfemale[, c('SleepHrsNight', 'Poverty', 'Pulse', 'SexNumPartnLife')]

grid_x1 <- quantile(X[,1], probs = c(0.1, 0.5, 0.9))
grid_x2 <- quantile(X[,2], probs = c(0.1, 0.5, 0.9))
grid_x3 <- quantile(X[,3], probs = c(0.5))
grid_x4 <- quantile(X[,4], probs = c(0.1, 0.5, 0.9))
grid_X <- expand.grid(grid_x1, grid_x2, grid_x3, grid_x4)

dim(X)

#Running of regression model
DPprior <- PYcalibrate(Ek = 3, n = n, discount = 0)
prior <- list(strength = DPprior$strength, discount = 0, hyper = FALSE,
              m0 = c(4.770, 6.885, 3.001, 74.78, 14.42), s0 = diag(c(10, 1, 1, 20, 50)),
              a0 = 2, b0 = 1.5)
mcmc <- list(niter = 11000, nburn = 1000)
output <- list(grid_x = grid_X, grid_y = grid_y, out_param = TRUE, out_type = 'FULL')
set.seed(1)
fit.female <- PYregression(y = Yfemale[,1], x = X, prior = prior, 
                           mcmc = mcmc, output = output)

summary(fit.female)
?PYregression
# Extract grid indices for plotting
num_poverty = length(grid_x2) 
num_partners = length(grid_x4) 

# Low sleep hours (10th percentile)
regplot1 <- data.frame(
  dens = as.vector(apply(fit.female$density[, 1:(num_poverty * num_partners), ], c(1, 2), mean)),
  qlow = as.vector(apply(fit.female$density[, 1:(num_poverty * num_partners), ], c(1, 2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.female$density[, 1:(num_poverty * num_partners), ], c(1, 2), quantile, probs = 0.975)),
  grid = rep(grid_y, num_poverty * num_partners),
  label = factor(rep(paste("Poverty = ", grid_x2, ", SexNumPartnLife = ", grid_x4), each = length(grid_y)),
                 level = rep(paste("Poverty = ", grid_x2, ", SexNumPartnLife = ", grid_x4)))
)

# Medium sleep hours (50th percentile)
regplot2 <- data.frame(
  dens = as.vector(apply(fit.female$density[, (num_poverty * num_partners + 1):(2 * num_poverty * num_partners), ], c(1, 2), mean)),
  qlow = as.vector(apply(fit.female$density[, (num_poverty * num_partners + 1):(2 * num_poverty * num_partners), ], c(1, 2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.female$density[, (num_poverty * num_partners + 1):(2 * num_poverty * num_partners), ], c(1, 2), quantile, probs = 0.975)),
  grid = rep(grid_y, num_poverty * num_partners),
  label = factor(rep(paste("Poverty = ", grid_x2, ", SexNumPartnLife = ", grid_x4), each = length(grid_y)),
                 level = rep(paste("Poverty = ", grid_x2, ", SexNumPartnLife = ", grid_x4)))
)

# High sleep hours (90th percentile) 
regplot3 <- data.frame(
  dens = as.vector(apply(fit.female$density[, (2 * num_poverty * num_partners + 1):(3 * num_poverty * num_partners), ], c(1, 2), mean)),
  qlow = as.vector(apply(fit.female$density[, (2 * num_poverty * num_partners + 1):(3 * num_poverty * num_partners), ], c(1, 2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.female$density[, (2 * num_poverty * num_partners + 1):(3 * num_poverty * num_partners), ], c(1, 2), quantile, probs = 0.975)),
  grid = rep(grid_y, num_poverty * num_partners),
  label = factor(rep(paste("Poverty = ", grid_x2, ", SexNumPartnLife = ", grid_x4), each = length(grid_y)),
                 level = rep(paste("Poverty = ", grid_x2, ", SexNumPartnLife = ", grid_x4)))
)

#This is a fix I attempted to get the plots to work
labels <- expand.grid(
  Poverty = paste("Poverty =", grid_x2, "   "),
  Partners = paste("SexPartners =", grid_x4, "   ")
)

regplot1$label <- factor(
  rep(apply(labels, 1, paste, collapse = ",     "), each = length(grid_y)),
  levels = apply(labels, 1, paste, collapse = ",     ")
)

regplot2$label <- regplot1$label
regplot3$label <- regplot1$label

# Low sleep hours
plot1 <- ggplot(regplot1) + 
  theme_bw() +
  geom_line(aes(x = grid, y = dens)) +
  geom_ribbon(aes(x = grid, ymin = qlow, ymax = qupp), fill = 'green', alpha = 0.3) +
  facet_wrap(~label, ncol = 3, nrow = 3, labeller = label_wrap_gen()) + 
  labs(x = 'Bad Mental Health Days', y = 'Density', title = 'Low Sleep Hours') +
  xlim(0,10)

# Medium sleep hours
plot2 <- ggplot(regplot2) + 
  theme_bw() +
  geom_line(aes(x = grid, y = dens)) +
  geom_ribbon(aes(x = grid, ymin = qlow, ymax = qupp), fill = 'red', alpha = 0.3) +
  facet_wrap(~label, ncol = 3, nrow = 3, labeller = label_wrap_gen()) + 
  labs(x = 'Bad Mental Health Days', y = 'Density', title = 'Medium Sleep Hours') +
  xlim(0,10)

# High sleep hours
plot3 <- ggplot(regplot3) + 
  theme_bw() +
  geom_line(aes(x = grid, y = dens)) +
  geom_ribbon(aes(x = grid, ymin = qlow, ymax = qupp), fill = 'blue', alpha = 0.3) +
  facet_wrap(~label, ncol = 3, nrow = 3, labeller = label_wrap_gen()) + 
  labs(x = 'Bad Mental Health Days', y = 'Density', title = 'High Sleep Hours') +
  xlim(0,10)

plot1
plot2
plot3
grid.arrange(plot1, plot3, 
             layout_matrix = matrix(1:2, 1, 2))
#Now we produce posterior plots for the estimators

betas.female <- fit.female$beta
post.betas1 <- matrix(NA, nrow = 10000, ncol = 5)
for(i in 1:10000){
  post.betas1[i,] <- betas.female[[i]][1,]
}

par(mfrow = c(2,3))
plot(density(post.betas1[,1]), main = '', lwd = 2, ylab = '', xlab = expression(beta[0]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas1[,2]), main = '', lwd = 2, ylab = '', xlab = expression(beta[1]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas1[,3]), main = '', lwd = 2, ylab = '', xlab = expression(beta[2]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas1[,4]), main = '', lwd = 2, ylab = '', xlab = expression(beta[3]),
     bty = 'l', col = 1, cex.lab = 1.3)
plot(density(post.betas1[,5]), main = '', lwd = 2, ylab = '', xlab = expression(beta[4]),
     bty = 'l', col = 1, cex.lab = 1.3)
mtext('Posterior Regression Coefficients', side = 3, line = -2, outer = TRUE,
      cex = 1.5)

#Produce equivalent plots but with nicer illustration

plot.female1 <- ggplot(as.data.frame(post.betas1[,1]), aes(x = post.betas1[, 1]))+
  geom_density(color = 'red', size = 1) +
  labs(x = expression(Beta[0]))
plot.female2 <- ggplot(as.data.frame(post.betas1[,2]), aes(x = post.betas1[, 2]))+
  geom_density(color = 'red', size = 1)+
  labs(x = expression(Beta[1]))
plot.female3 <- ggplot(as.data.frame(post.betas1[,3]), aes(x = post.betas1[, 3]))+
  geom_density(color = 'red', size = 1) +
  labs(x = expression(Beta[2]))
plot.female4 <- ggplot(as.data.frame(post.betas1[,4]), aes(x = post.betas1[, 4]))+
  geom_density(color = 'red', size = 1)+
  labs(x = expression(Beta[3]))
plot.female5 <- ggplot(as.data.frame(post.betas1[,5]), aes(x = post.betas1[, 5]))+
  geom_density(color = 'red', size = 1)+
  labs(x = expression(Beta[4]))
grid.arrange(plot.female1, plot.female2, plot.female3, plot.female4, plot.female5,
             layout_matrix = matrix(1:5, 1,5),
             top = 'Denity plots of regression roefficients for female')


###########################
citation('NHANES')