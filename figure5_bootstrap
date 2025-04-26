# analysis of zhang dataset
# uses zhang commands
# assumes you unzip zhang.data.tar.gz to zhang directory

source("~/Desktop/BIOSTAT/msage/R/msage.R")
# use this as the source for the msage code

# don't library msage, it's gone :/
f <- dir("~/Desktop/BIOSTAT/bioinformatics_23_21_2881_s1/zhang", "txt", full.names = TRUE)
a <- read.sage(f[c(6,7,11,12)])
print(a)
k <- rowSums(a$data) > 2
zhang <- list(data = a$data[k,], groups = c(1,1,2,2), lib.size = a$lib.size)

alpha <- alpha.approxeb(zhang)
ms <- msage(zhang, alpha = alpha$alpha)
lu <- runLu(zhang, other.r = ms$r)


set.seed(123)
B <- 100
bootstrap_results <- matrix(NA, nrow = nrow(zhang$data), ncol = B)

for (i in 1:nrow(zhang$data)) {
  counts_tag <- zhang$data[i, ]
  
  for (b in 1:B) {
    idx <- sample(1:4, replace = TRUE)
    counts_boot <- counts_tag[idx]
    lib_size_boot <- zhang$lib.size[idx]
    group_boot <- zhang$groups[idx]
    
    zhang_boot <- list(
      data = matrix(counts_boot, nrow = 1),
      groups = group_boot,
      lib.size = lib_size_boot
    )
    
    # NEW: robust validity checks
    valid <- TRUE
    if (any(is.na(counts_boot)) || any(counts_boot < 0)) valid <- FALSE
    if (length(unique(group_boot[group_boot == 1])) < 1 || 
        length(unique(group_boot[group_boot == 2])) < 1) valid <- FALSE
    if (sum(zhang_boot$data[1, group_boot == 1]) <= 1 || 
        sum(zhang_boot$data[1, group_boot == 2]) <= 1) valid <- FALSE
    
    if (valid) {
      tryCatch({
        alpha_boot <- alpha.approxeb(zhang_boot)$alpha
        ms_boot <- msage(zhang_boot, alpha = alpha_boot)
        bootstrap_results[i, b] <- 1 / ms_boot$r
      }, error = function(e) {
        bootstrap_results[i, b] <- NA  # catch any crash silently
      })
    } else {
      bootstrap_results[i, b] <- NA
    }
  }
}

bootstrap_means <- rowMeans(bootstrap_results, na.rm = TRUE)
# BOOTSTRAPPING ADDITION ENDS HERE

# adjust p-values
adj.p <- p.adjust(ms$exact, method = "fdr")
k <- (adj.p < 0.05)

all <- cbind(zhang$data, adj.p, sign(rowSums(zhang$data[,1:2]) - rowSums(zhang$data[,3:4])))
write.table(all[k,], "all.csv", sep = ",")

# now plot
pdf("figure5bootstrap.pdf", width = 12, height = 6)
par(mfrow = c(1,2))

# Plot A
plot(ms$ps$p, 1/ms$r, log = "x", xlab = "lambda", ylab = "phi", 
     pch = 19, cex = 0.5, cex.axis = 1.5, cex.lab = 1.5)
# optional: overlay bootstrap means
points(ms$ps$p, bootstrap_means, col = "pink", pch = 4, cex = 0.6)
mtext("A", adj = 0, cex = 2, padj = -0.75)

# Plot B
plot(log2(ms$ps$p1) + log2(ms$ps$p2), log2(ms$ps$p2/ms$ps$p1), 
     ylim = c(-7,7), xlim = c(-35,-10), pch = "+", 
     xlab = "log proportion", ylab = "log Fold Change", 
     cex.axis = 1.5, cex.lab = 1.5)
mtext("B", adj = 0, cex = 2, padj = -0.75)
points(log2(ms$ps$p1)[k] + log2(ms$ps$p2)[k], log2(ms$ps$p2/ms$ps$p1)[k], 
       col = "grey", pch = 19)

dev.off()

