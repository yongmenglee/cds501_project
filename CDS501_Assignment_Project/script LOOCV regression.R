getwd()
corona <- read.csv('corona.csv', sep=',', header=T) 

install.packages("caret") 
library(caret)

#Linear regression using LOOCV (X.alignment.length.)
train.control <- trainControl(method = "LOOCV")
modellm <- train(X.bit.score. ~ X.alignment.length., data = corona, method = "lm", trControl = train.control)
print(modellm)