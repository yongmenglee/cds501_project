getwd()
corona <- read.csv('corona.csv', sep=',', header=T) 


install.packages("caret") 
library(caret)
install.packages ("tidyverse")
library(tidyverse)

data("corona")
sample_n(X.alignment.length, 3)

set.seed(123)
training.samples <- corona %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- corona[training.samples, ]
test.data <- corona[-training.samples, ]
# Build the model
model <- lm( X.alignment.length ~., data = train.data)
# Make predictions and compute the R2, RMSE and MAE
predictions <- model %>% predict(test.data)
data.frame( R2 = R2(predictions, test.data$ X.alignment.length),
            RMSE = RMSE(predictions, test.data$ X.alignment.length),
            MAE = MAE(predictions, test.data$ X.alignment.length))

RMSE(predictions, test.data$)/mean(test.data$X.alignment.length)