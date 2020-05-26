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
m <- cor(corona1)
# view(m)

# #Plot histogram for each variale
# install.packages("ggplot2")
# install.packages("tidyr")

# Check for missing values for each column
sapply(corona1, function(x) sum(is.na(x)))
# %_identity alignment_length       mismatches        gap_opens          q.start 
# 0                0                0                0                0 
# q.end          s.start            s.end           evalue        bit_score 
# 0                0                0                0                0 

# ---
# For bivariate linear regression, we will use alignment_length
# to predict the bit_score.

# Step 1 - Create new data frame
corona.bivar <- data.frame(
  alignment_length = corona1$alignment_length,
  bit_score = corona1$bit_score
)

# Check if new data frame is implemented correctly
head(corona.bivar)

# Plot new data frame
plot(corona.bivar$alignment_length, corona.bivar$bit_score)


# Step 2 - Create a linear regression model
corona.bivar.regression <- lm(bit_score ~ alignment_length,
                              data = corona.bivar)

# Step 3 - Summary of the regression report.
summary(corona.bivar.regression)

# Step 4 - Add regression line to the x/y scatter plot.
abline(corona.bivar.regression, col = "blue")

# ---

corona.multivar <- data.frame(
  pct_identity = corona1[,"%_identity"],
  alignment_length = corona1$alignment_length,
  gap_opens = corona1$gap_opens,
  q.start = corona1$q.start,
  bit_score = corona1$bit_score
)

colnames(corona.multivar)
# corona.multivar

corona.multivar.regression <- lm(bit_score ~ pct_identity + alignment_length + gap_opens + q.start,
                              data = corona.multivar)

summary(corona.multivar.regression)


