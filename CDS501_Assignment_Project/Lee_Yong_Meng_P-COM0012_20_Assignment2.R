# CDS501 Assignment 2
# Group 3 - Coronavirus

# Read data from CSV file
corona <- read.table('MN997409.1-4NY0T82X016-Alignment-HitTable.csv', sep=',')

# Add header to data frame
names(corona) <- c('X.query.acc.ver.', 'X.subject_acc.ver.', 'X.pct.identity.', 
                   'X.alignment.length.', 'X.mismatches.', 'X.gap.opens.',
                   'X.q.start.', 'X.q.end.', 'X.s.start.', 'X.s.end.', 'X.evalue.', 'X.bit.score.')
summary(corona)
# Output::
# X.query.acc.ver.  X.subject_acc.ver. X.pct.identity.  X.alignment.length. X.mismatches.     X.gap.opens.   
# MN997409.1:263     AP006561.1:  4      Min.   : 77.56   Min.   : 1603       Min.   :   0.0   Min.   :  0.00  
# AY278554.2:  4      1st Qu.: 80.05   1st Qu.: 1925       1st Qu.: 142.0   1st Qu.: 12.00  
# AY282752.2:  4      Median : 82.30   Median : 5417       Median : 359.0   Median : 35.00  
# AY283796.1:  4      Mean   : 86.06   Mean   :10711       Mean   : 919.2   Mean   : 57.82  
# AY291451.1:  4      3rd Qu.: 90.19   3rd Qu.:17716       3rd Qu.: 989.0   3rd Qu.: 68.00  
# AY304486.1:  4      Max.   :100.00   Max.   :29882       Max.   :2952.0   Max.   :172.00  
# (Other)   :239                                                                            
# X.q.start.       X.q.end.       X.s.start.       X.s.end.       X.evalue.  X.bit.score.  
# Min.   :    1   Min.   : 1923   Min.   :    1   Min.   : 1672   Min.   :0   Min.   : 1011  
# 1st Qu.:   16   1st Qu.:21577   1st Qu.:    1   1st Qu.:21489   1st Qu.:0   1st Qu.: 2101  
# Median : 3956   Median :27910   Median : 3875   Median :27783   Median :0   Median : 3936  
# Mean   :11296   Mean   :21971   Mean   :11213   Mean   :21889   Mean   :0   Mean   :14240  
# 3rd Qu.:22539   3rd Qu.:29876   3rd Qu.:22429   3rd Qu.:29729   3rd Qu.:0   3rd Qu.:15175  
# Max.   :28257   Max.   :29882   Max.   :28137   Max.   :30256   Max.   :0   Max.   :55182  

# Remove first and second columns
corona <- corona[-c(1, 2)]
head(corona)
# Output::
# X.pct.identity. X.alignment.length. X.mismatches. X.gap.opens. X.q.start. X.q.end. X.s.start. X.s.end. X.evalue.
# 1         100.000               29882             0            0          1    29882          1    29882         0
# 2          99.990               29882             3            0          1    29882          1    29882         0
# 3          99.990               29882             3            0          1    29882          1    29882         0
# 4          99.990               29882             3            0          1    29882          1    29882         0
# 5          99.990               29882             3            0          1    29882          1    29882         0
# 6          99.993               29878             2            0          4    29881          1    29878         0
# X.bit.score.
# 1        55182
# 2        55166
# 3        55166
# 4        55166
# 5        55166
# 6        55164

# ------------------------------

# - Part 1 - Cross Validation
# Import "caret" package for Cross Validation
install.packages("caret") 
library(caret)

# 1.1 K-fold Cross Validation
# - Waiting for Ameer's part.
# Define training control
set.seed(123)
train.control.kfold <- trainControl(method = "cv", number = 5)
model.kfold <- train(X.bit.score. ~ X.alignment.length., data = corona, method = "lm", trControl = train.control.kfold)

print(model.kfold)
# Output::
# Linear Regression 
# 
# 263 samples
# 1 predictor
# 
# No pre-processing
# Resampling: Cross-Validated (5 fold) 
# Summary of sample sizes: 209, 212, 209, 211, 211 
# Resampling results:
#   
#   RMSE      Rsquared  MAE     
# 6320.651  0.891921  5080.146
# 
# Tuning parameter 'intercept' was held constant at a value of TRUE


# 1.2 LOOCV (Leave-one-out Cross Validation)

# Linear regression using LOOCV (X.alignment.length.)
train.control.loo <- trainControl(method = "LOOCV")
model.loo <- train(X.bit.score. ~ X.alignment.length., data = corona, method = "lm", trControl = train.control.loo)

# Print the outcome of loocv
print(model.loo)
# Output::
# Linear Regression 
# 
# 263 samples
# 1 predictor
# 
# No pre-processing
# Resampling: Leave-One-Out Cross-Validation 
# Summary of sample sizes: 262, 262, 262, 262, 262, 262, ... 
# Resampling results:
#   
#   RMSE      Rsquared   MAE    
# 6350.328  0.8904952  5102.54
# 
# Tuning parameter 'intercept' was held constant at a value of TRUE

# ----------------------------------------

# - Part 2 - Regression

# Import libraries for data visualization - plotting 
install.packages("ggplot2")
install.packages("tidyr")
library(ggplot2)
library(tidyr)

# Step 0.1 - Split data into dTrain and dTest
# Set seed value - ensure the result is always reproducible
set.seed(0)

# Generate new column with random numbers from 0 - 1
corona$gp <- runif(dim(corona)[1])
head(corona)
# Output::
# X.pct.identity. X.alignment.length. X.mismatches. X.gap.opens. X.q.start. X.q.end. X.s.start. X.s.end. X.evalue.
# 1         100.000               29882             0            0          1    29882          1    29882         0
# 2          99.990               29882             3            0          1    29882          1    29882         0
# 3          99.990               29882             3            0          1    29882          1    29882         0
# 4          99.990               29882             3            0          1    29882          1    29882         0
# 5          99.990               29882             3            0          1    29882          1    29882         0
# 6          99.993               29878             2            0          4    29881          1    29878         0
# X.bit.score.        gp
# 1        55182 0.8966972
# 2        55166 0.2655087
# 3        55166 0.3721239
# 4        55166 0.5728534
# 5        55166 0.9082078
# 6        55164 0.2016819

# Set size of test set to 10%
test.size <- 0.1

# Split data into training and test sets based on test.size
dTrain <- subset(corona, corona$gp > test.size)
dTest <- subset(corona, corona$gp < test.size)

# Step 0.2 - Define functions R-Squared and RMSE
# - R-squared
rsq <- function(y.true, pred.y) { 
  1 - sum((y.true - pred.y)^2) / sum((y.true - mean(y.true))^2) 
}

# - RMSE
rmse <- function(y.true, pred.y) { 
  sqrt(mean( (y.true - pred.y)^2 ))
}

# 2.1: Bivariate regression
# Step 1 - Visualize two variables on scatter plot.
plot(dTrain$X.alignment.length., dTrain$X.bit.score., 
     xlab = "alignment_length", ylab = "bit_score") +
  title("bit_score vs alignment_length")

# Step 2 - Create bivariate linear regression model
model.bivar <- lm(X.bit.score. ~ X.alignment.length., data = dTrain)
model.bivar

# Step 3 - Summary report of the biviariate regression on training set
summary(model.bivar)
# Output::
# Call:
#   lm(formula = X.bit.score. ~ X.alignment.length., data = dTrain)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -11246  -1247   2148   3457   7760 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -4.169e+03  5.666e+02  -7.359 2.72e-12 ***
#   X.alignment.length.  1.726e+00  3.789e-02  45.565  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6306 on 248 degrees of freedom
# Multiple R-squared:  0.8933,	Adjusted R-squared:  0.8929 
# F-statistic:  2076 on 1 and 248 DF,  p-value: < 2.2e-16

# Step 4 - Add regression line to the x/y scatter plot.
abline(model.bivar, col = "blue")

# Step 5 - Prediction on test set
dTest$pred.bit.score. <- predict(model.bivar, newdata = dTest)
# dTrain$pred.bit.score. <- predict(model.bivar, newdata = dTrain)

# Step 6 - Calculate R-Squared
rsq(dTrain$X.bit.score., dTrain$pred.bit.score.)
# Output::
# [1] 0.8932972

# Step 7 - Calculate RMSE
rmse(dTest$X.bit.score., dTest$pred.bit.score.)
# Output::
# [1] 6597.397


# 2.2 Multivariate regression
# Step 1 - Visualize variables for each feature with target on scatter plots.
plot(dTrain$X.alignment.length., dTrain$X.bit.score.,
     xlab = "alignment_length", ylab = "bit_score") +
  title("bit_score vs alignment_length")
cor(dTrain$X.alignment.length., dTrain$X.bit.score.)
# Output::
# [1] 0.945144

plot(dTrain$X.alignment.length., dTrain$X.mismatches.,
     xlab = "mismatches", ylab = "bit_score") +
  title("bit_score vs mismatches")
cor(dTrain$X.alignment.length., dTrain$X.mismatches.)
# Output::
# [1] 0.2043152

# Step 1 - Create multivariate linear regression model
model.multivar <- lm(X.bit.score. ~ X.alignment.length. + X.mismatches., data = dTrain)
model.multivar

# Step 2 - Summary report of the multiviariate regression on training set
summary(model.multivar)
# Output::
# Call:
#   lm(formula = X.bit.score. ~ X.pct.identity. + X.alignment.length. + 
#        X.gap.opens. + X.q.start., data = dTrain)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8499.8   -50.8    60.4   141.8  3905.8 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          8.912e+03  1.248e+03   7.144 1.03e-11 ***
#   X.pct.identity.     -9.420e+01  1.537e+01  -6.130 3.49e-09 ***
#   X.alignment.length.  1.858e+00  1.114e-02 166.767  < 2e-16 ***
#   X.gap.opens.        -1.115e+02  1.346e+00 -82.875  < 2e-16 ***
#   X.q.start.           7.578e-04  5.417e-03   0.140    0.889    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 655.6 on 245 degrees of freedom
# Multiple R-squared:  0.9989,	Adjusted R-squared:  0.9988 
# F-statistic: 5.371e+04 on 4 and 245 DF,  p-value: < 2.2e-16

# Step 3 -  Prediction on test set
dTest$pred.mul.bit.score.<- predict(model.multivar, newdata = dTest)
# dTrain$pred.mul.bit.score. <- predict(model.multivar, newdata = dTrain)

# Step 4 - Calculate R-squared
rsq(dTest$X.bit.score., dTest$pred.mul.bit.score.)
# Output::
# [1] 0.9906639

# Step 5 - Calculate RMSE
rmse(dTest$X.bit.score., dTest$pred.mul.bit.score.)
# Output:: 
# [1] 1775.943
