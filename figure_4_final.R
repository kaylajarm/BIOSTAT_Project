### This code is used to generate the false discovery plots from figure 4. I manually changed the b and n values since I didn't really use functions in this code ###

set.seed(911)
library(edgeR)
library(MASS)
source("C:/Users/megan/msage.R")
total_tags <- 10000 #there are total of 10,000 tags in the simulation
de_genes_n <- total_tags * 0.1 # 10% are DE genes
# making a vector of trues and falses regarding if they are differentially expressed
gene_types <- c(rep(TRUE, de_genes_n), rep(FALSE, total_tags - de_genes_n))
b <- 4  # #implanted difference for the small test
n <- 2  # replicates per group

# mean gene expressions for both groups
mu <- rep(10, total_tags) # #starting average gene expression
mu1 <- mu; mu1[gene_types] <- mu[gene_types] / sqrt(b) #DE genes group 1
mu2 <- mu; mu2[gene_types] <- mu[gene_types] * sqrt(b) #DE genes group 2

# creating the random library sizes
library_sizes <- matrix(runif(2 * n, 30000, 90000), nrow = 2) # paper says the range of library sizes simulated were from 30,000 to 90,000

# actually simulating the data
counts <- matrix(0, nrow = total_tags, ncol = 2 * n)
# creates random dispersion values since they wouldn't be the same between groups
disp1 <- rgamma(total_tags, shape = 0.85, scale = 0.5)
disp2 <- rgamma(total_tags, shape = 0.85, scale = 0.5)
# for loop to fill the counts matrix that will hold all the simulated data for the two groups, dividing by 1e6 to scale the data so that it would be counts per library, which uses counts per million scaling which takes into account the multiplied by the library size condition of the simulation. this just makes the simulation proportional to library size, and it does not work without scaling it like this
for (tag in 1:total_tags) {
  counts[tag, 1:n] <- rnbinom(n, mu = mu1[tag] * library_sizes[1, 1:n] / 1e6, size = 1 / disp1[tag])
  counts[tag, (n+1):(2*n)] <- rnbinom(n, mu = mu2[tag] * library_sizes[2, 1:n] / 1e6, size = 1 / disp2[tag])
}
# running edgeR to actually determine the differentially expressed genes
groups <- rep(1:2, each = n)
dge <- DGEList(counts = counts, group = groups) # Assembles a DGEList object from its components
dge <- normLibSizes(dge) #Calculate scaling factors to convert the raw library sizes for a set of sequenced samples into normalized effective library sizes.
dge <- estimateDisp(dge,trend.method = "loess") # Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags -> needed to run the actual model later

# Exact.EB
et <- exactTest(dge) # runs the moderated test from the paper
pvals_exact <- et$table$PValue # extracts the p-values from the exact test table

# Wald.EB
fit <- glmFit(dge) #Fit a negative binomial generalized log-linear model to the read counts for each gene.
lrt <- glmLRT(fit) #Given genewise generalized linear model fits, conduct likelihood ratio tests for a given coefficient or coefficient contrast.
pvals_wald <- lrt$table$PValue # extracts the pvalues from the Wald test

#function to replace a defunct package that estimates the dispersion of a Poisson model
luEstDisp <- function(y, mu, df.residual) {
  res <- (y - mu) / sqrt(mu)
  dispersion <- sum(res^2) / df.residual
  return(max(dispersion, 1e-6))
}

# Wald.PL, from the msage package but edited because it was trying to use defunct packages
runLu <- function(x, other.r = NULL) {
  wald <- rep(NA, nrow(x$data))
  if (!is.null(other.r)) wald2 <- rep(NA, nrow(x$data))
  phi <- rep(NA, nrow(x$data))
  o <- log(x$lib.size)
  for (i in 1:nrow(x$data)) {
    y <- as.numeric(x$data[i, ])
    f <- glm(y ~ offset(o) + x$group, family = poisson(link = "log"))
    mu <- f$fitted.values
    phi_hat <- luEstDisp(y, mu, f$df.residual)
    phi[i] <- phi_hat
    if (phi_hat < 1e-5) {
      ff <- f
    } else {
      ff <- glm(y ~ offset(o) + x$group, family = negative.binomial(theta = 1 / phi_hat))
    }
    wald[i] <- summary(ff)$coef[2, 3]
    if (!is.null(other.r)) {
      ff2 <- glm(y ~ offset(o) + x$group, family = negative.binomial(theta = other.r[i]))
      wald2[i] <- summary(ff2)$coef[2, 3]
    }
    if (i %% 100 == 0) cat("Lu - done:", i, "\n")
  }
  if (!is.null(other.r)) {
    return(list(wald = wald, wald2 = wald2, phi = phi))
  } else {
    return(list(wald = wald, phi = phi))
  }
}

#reorganizing the data to run the Lu test, since it does not take a DGE object
sim_data <- list(data = counts, group = factor(groups),lib.size = colSums(counts))

lu_result <- runLu(sim_data)
pvals_lu <- 2 * pnorm(-abs(lu_result$wald))  #converts Wald stats to p-values

#preparing and plotting the false discovery curves. order allows us to determine what has the least pvalue and the greatest to determine which is the most differentially expressed. cum sums to give false discoveries from the top genes
fd_exact <- cumsum(!gene_types[order(pvals_exact)])
fd_wald  <- cumsum(!gene_types[order(pvals_wald)])
fd_lu    <- cumsum(!gene_types[order(pvals_lu)])

#plotting: each fd is from range 1:1000 since we just want to see the 1000 top ranked genes, log = "y" makes it a logarithmic scale, ylim finds the max of the different fd curves because when I tried to set it to 500, it wouldn't populate correctly
plot(1:1000, fd_exact[1:1000], type = 'l', col = 'blue', lwd = 2,
     xlim = c(0, 1000), ylim = c(1, max(fd_exact[1:1000], fd_wald[1:1000], fd_lu[1:1000])),
     log = "y", xlab = "Top ranked genes", ylab = "False discoveries",
     main = paste0("n1=n2=", n, " per group, b = ", b))

lines(1:1000, fd_wald[1:1000], col = 'firebrick1', lwd = 2)
lines(1:1000, fd_lu[1:1000], col = 'forestgreen', lwd = 2)

legend("topleft", legend = c("Exact.EB", "Wald.EB", "Wald.PL (Lu)"),
       col = c("blue", "firebrick1", "forestgreen"), lwd = 2, bty = "n")

