# CDS501 Assignment 2

# Create distribution to check on the missing value
# setwd("~/R Project")
corona <- read.table('MN997409.1-4NY0T82X016-Alignment-HitTable.csv', sep=',')
names(corona) <- c('query_acc.ver', 'subject_acc.ver', '%_identity', 
                   'alignment_length', 'mismatches', 'gap_opens',
                   'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit_score')
summary(corona)

# Calculate the correlation
corona1 <- corona[, c(3,4,5,6,7,8,9,10,11,12)]
# m <- cor(corona1)
# view(m)

# #Plot histogram for each variale
# Uncomment for first time running code
# install.packages("ggplot2")
# install.packages("tidyr")

library(ggplot2)
library(tidyr)


# Check for missing values for each column
sapply(corona1, function(x) sum(is.na(x)))
# %_identity alignment_length       mismatches        gap_opens          q.start 
# 0                0                0                0                0 
# q.end          s.start            s.end           evalue        bit_score 
# 0                0                0                0                0 

# ---
# For bivariate linear regression, we will use alignment_length
# to predict the bit_score.

# Hereby showing two approaches to perform regression

# Approach 1: Without splitting

# 1.1: Bivariate regression
# Step 1 - Visualize two variables on scatter plot.
plot(corona1$alignment_length, corona1$bit_score) +
  title("bit_score vs alignment_length")

# Step 2 - Create a linear regression model
corona.bivar.lm <- lm(bit_score ~ alignment_length,
                      data = corona1)

# Step 3 - Summary of the regression report.
summary(corona.bivar.lm)
# Output::
# Call:
#   lm(formula = bit_score ~ alignment_length, data = corona1)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -11150  -1179   2215   3524   7878 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -4.233e+03  5.565e+02  -7.606 5.11e-13 ***
#   alignment_length  1.725e+00  3.708e-02  46.512  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6321 on 261 degrees of freedom
# Multiple R-squared:  0.8923,	Adjusted R-squared:  0.8919 
# F-statistic:  2163 on 1 and 261 DF,  p-value: < 2.2e-16

# Step 4 - Add regression line to the x/y scatter plot.
abline(corona.bivar.lm, col = "blue")

# ----------

# 1.2: Multivariate regression
corona.multivar.lm <- lm(bit_score ~ `%_identity` + alignment_length + gap_opens + q.start,
                         data = corona1)

summary(corona.multivar.lm)
# Output::
# Call:
#   lm(formula = bit_score ~ `%_identity` + alignment_length + gap_opens + 
#        q.start, data = corona1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8445.6   -44.5    85.2   205.1  3921.2 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       8.571e+03  1.395e+03   6.144 3.03e-09 ***
#   `%_identity`     -9.020e+01  1.718e+01  -5.249 3.20e-07 ***
#   alignment_length  1.854e+00  1.247e-02 148.678  < 2e-16 ***
#   gap_opens        -1.113e+02  1.502e+00 -74.065  < 2e-16 ***
#   q.start           8.088e-04  6.036e-03   0.134    0.894    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 752.1 on 258 degrees of freedom
# Multiple R-squared:  0.9985,	Adjusted R-squared:  0.9985 
# F-statistic: 4.275e+04 on 4 and 258 DF,  p-value: < 2.2e-16

# ---

# Approach 2: Splitting data into dTrain and dTest
# Step 0.1 - Split data into dTrain and dTest
set.seed(0)corona1$gp = runif(dim(corona1)[1])
head(corona1)

test.size<- 0.1

dTrain <- subset(corona1, corona1$gp > test.size)
dTest <- subset(corona1, corona1$gp < test.size)

# Step 0.2 - Define functions R-Squared and RMSE
rsq <- function(y.true, pred.y) { 
  1 - sum((y.true - pred.y)^2) / sum((y.true - mean(y.true))^2) 
}

rmse <- function(y.true, pred.y) { 
  sqrt(mean( (y.true - pred.y)^2 ))
}

# 2.1: Bivariate regression
# Step 1 - Visualize two variables on scatter plot.
plot(dTrain$alignment_length, dTrain$bit_score) +
  title("bit_score vs alignment_length")

# Step 2 - Create a linear regression model
model.bivar <- lm(bit_score ~ alignment_length,
                  data = dTrain)

# Step 3 - Summary of the regression report.
summary(model.bivar)
# Output::
# Call:
#   lm(formula = bit_score ~ alignment_length, data = dTrain)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -11700  -1349   2148   3463   6959 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -4.222e+03  5.654e+02  -7.467 1.62e-12 ***
#   alignment_length  1.755e+00  3.743e-02  46.894  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6129 on 234 degrees of freedom
# Multiple R-squared:  0.9038,	Adjusted R-squared:  0.9034 
# F-statistic:  2199 on 1 and 234 DF,  p-value: < 2.2e-16

# Step 4 - Add regression line to the x/y scatter plot.
abline(model.bivar, col = "blue")

# Step 5 - Prediction on test set
dTest$pred.bit_score <- predict(model.bivar, newdata = dTest)
dTrain$pred.bit_score <- predict(model.bivar, newdata = dTrain)

# Step 6 - Calculate R-Squared
rsq(dTrain$bit_score, dTrain$pred.bit_score)
rsq(dTest$bit_score, dTest$pred.bit_score)
# Output:: might be different across each run

# Step 7 - Calculate RMSE
rmse(dTrain$bit_score, dTrain$pred.bit_score)
rmse(dTest$bit_score, dTest$pred.bit_score)
# Output:: might be different across each run

# cor(dTest$bit_score, dTest$pred.bit_score)
ggplot(data = dTest, aes(x = pred.bit_score, y = bit_score)) +
  geom_point(alpha=0.2, color="black") +
  geom_smooth(aes(x = pred.bit_score, y = bit_score), color="black") +
  geom_line(aes(x = bit_score, y = bit_score), color="blue", linetype=2)


# 2.2 Multivariate regression
# Step 1
model.multivar <- lm(bit_score ~ `%_identity` + alignment_length + gap_opens + q.start,
                     data = dTrain)

# Step 2
summary(model.multivar)
# Output::
# Call:
#   lm(formula = bit_score ~ `%_identity` + alignment_length + gap_opens + 
#        q.start, data = dTrain)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8420.7   -40.7   106.8   222.2  3924.5 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       8.284e+03  1.585e+03   5.226 3.96e-07 ***
#   `%_identity`     -8.718e+01  1.952e+01  -4.467 1.26e-05 ***
#   alignment_length  1.853e+00  1.411e-02 131.284  < 2e-16 ***
#   gap_opens        -1.110e+02  1.708e+00 -64.977  < 2e-16 ***
#   q.start           1.115e-03  6.979e-03   0.160    0.873    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 802 on 225 degrees of freedom
# Multiple R-squared:  0.9984,	Adjusted R-squared:  0.9983 
# F-statistic: 3.437e+04 on 4 and 225 DF,  p-value: < 2.2e-16

# Step 3
dTest$pred.mul.bit_score <- predict(model.multivar, newdata = dTest)

dTrain$pred.mul.bit_score <- predict(model.multivar, newdata = dTrain)

# Step 4
rsq(dTrain$bit_score, dTrain$pred.mul.bit_score)
rsq(dTest$bit_score, dTest$pred.mul.bit_score)
# Output:: might be different across each run

# Step 5
rmse(dTrain$bit_score, dTrain$pred.mul.bit_score)
rmse(dTest$bit_score, dTest$pred.mul.bit_score)
# Output:: might be different across each run

