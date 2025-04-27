
library(edgeR)
library(ggplot2)

# setting up the simulation requirements
set.seed(911)
tags <- 1000 #1000 tags simulated
n<- 4 # number of libraries
disp<- 0.42 #dispersion (φ)
delta <- disp/ (1 + disp) # to be used in the condLogLikDerDelta later
r<- 1/disp   # NB size parameter

# simulating the counts using rnbinom,set mu of 10 for each tag
counts <- matrix(rnbinom(G*n, mu = 10, size = r),ncol = n) # n column is 4 because there are 4 libraries
# total number of counts, going to use this later for the total tag count that will be plotted
z_total <- rowSums(counts)

#looping over each tag and determining the observed info
for (g in 1:tags) {
  y1 <- counts[g, 1:2] # z1
  y2 <- counts[g, 3:4] # z2
  # making a matrix of the two vectors to be used in the condLogLikeDerDelta function
  y <- cbind(y1,y2)
  # computes the second derivative which is the observed info to be plotted
  Jg <- -condLogLikDerDelta(counts, delta, der = 2) 
  Jg_vec_two <- Jg # storing Jg in a vector so that it can later be plotted by each count
}

# plotting
df <- data.frame(z_total, Jg = Jg_vec_two)
ggplot(df, aes(x = z_total, y = Jg)) +
geom_point() + geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, color = "blue") + theme(panel.background = element_blank()) + labs(x = "Total tag count",y = "Observed information (Jg)",) + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + scale_y_continuous(breaks = seq(0,100,by=20), limits = c(0,100)) + scale_x_continuous(breaks = seq(0,100,by=20))
