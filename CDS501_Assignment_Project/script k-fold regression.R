getwd()
# corona <- read.csv('corona.csv', sep=',', header=T) 
corona <- read.csv('MN997409.1-4NY0T82X016-Alignment-HitTable.csv', sep=',', header=T) 
# MN997409.1-4NY0T82X016-Alignment-HitTable.csv

install.packages("caret") 
library(caret)
install.packages ("tidyverse")
library(tidyverse)

# Corona is not a built-in dataset in the R package lol.
# data("corona")
# data("corona")
names(corona) <- c('X.query.acc.ver.', 'X.subject_acc.ver.', 'X.pct.identity.', 
                   'X.alignment.length.', 'X.mismatches.', 'X.gap.opens.',
                   'X.q.start.', 'X.q.end.', 'X.s.start.', 'X.s.end.', 'X.evalue.', 'X.bit.score.')

sample_n(corona$X.alignment.length., 3)

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