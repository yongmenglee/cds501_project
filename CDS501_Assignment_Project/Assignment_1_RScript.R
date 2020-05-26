#Create distribution to check on the missing value
setwd("~/R Project")
corona <- read.table('MN997409.1-4NY0T82X016-Alignment-HitTable.csv', sep=',')
names(corona) <- c('query_acc.ver', 'subject_acc.ver', '%_identity', 
                   'alignment_length', 'mismatches', 'gap_opens',
                   'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit_score')
summary(corona)

#Calculate the correlation
corona1 <- corona[, c(3,4,5,6,7,8,9,10,11,12)]
m <- cor(corona1)
view(m)

# #Plot histogram for each variale
# install.packages("ggplot2")
# install.packages("tidyr")
library(ggplot2)
library(tidyr)

#specify size here
corona1 %>% gather() %>% head()
histogram <- ggplot(gather(corona1), aes(value)) + geom_histogram(bins = 10) + facet_wrap(~key, scales = 'free_x')

#Plot scatter matrix
scatter <- pairs(corona1)

#Plot individual scatter 
install.packages("car")
library("car") 

scatterplot(mismatches ~ gap_opens, data = corona) 
scatterplot(s.start ~ q.start, data = corona) 
scatterplot(s.end ~ q.end, data = corona) 
scatterplot(corona$alignment_length ~ corona$bit_score, data = corona)




