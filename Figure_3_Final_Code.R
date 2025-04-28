# Load necessary libraries
library(edgeR)
library(pROC)
library(MASS)
library(ggplot2)
library(gridExtra)

# Function to simulate count data
simulate_counts_adjusted <- function(total_tags = 10000, de_genes_n = 5000, b = 4, n = 5, 
                                     scale_factor = 1e6, disp_value = 0.2) {
  counts <- matrix(0, nrow = total_tags, ncol = 2 * n)
  gene_types <- c(rep(TRUE, de_genes_n), rep(FALSE, total_tags - de_genes_n))
  
  mu <- rep(10, total_tags)
  mu1 <- mu; mu1[gene_types] <- mu[gene_types] / sqrt(b)
  mu2 <- mu; mu2[gene_types] <- mu[gene_types] * sqrt(b)
  
  library_sizes <- matrix(runif(2 * n, 30000, 90000), nrow = 2)
  
  for (tag in 1:total_tags) {
    disp <- disp_value
    counts[tag, 1:n] <- rnbinom(n, mu = mu1[tag] * library_sizes[1, 1:n] / scale_factor, size = 1 / disp)
    counts[tag, (n + 1):(2 * n)] <- rnbinom(n, mu = mu2[tag] * library_sizes[2, 1:n] / scale_factor, size = 1 / disp)
  }
  
  list(counts = counts, groups = rep(1:2, each = n), gene_types = gene_types)
}

# Function to calculate p-values for different methods
edgeR_scores <- function(counts, groups, test_type = "Exact.EB") {
  dge <- DGEList(counts = counts, group = groups)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)
  
  if (test_type == "Exact.EB") {
    et <- exactTest(dge)
    return(et$table$PValue)
  } else if (test_type == "Wald.EB") {
    fit <- glmFit(dge)
    lrt <- glmLRT(fit)
    return(lrt$table$PValue)
  }
}

# Lu's Poisson Wald method
runLu <- function(counts, groups) {
  n_samples <- ncol(counts)
  lib_size <- colSums(counts)
  group_factor <- factor(c(rep(1, n_samples / 2), rep(2, n_samples / 2)))
  wald_stats <- rep(NA, nrow(counts))
  
  for (i in 1:nrow(counts)) {
    y <- counts[i, ]
    fit <- glm(y ~ offset(log(lib_size)) + group_factor, family = poisson(link = "log"))
    wald_stats[i] <- summary(fit)$coef[2, 3]
  }
  
  return(wald_stats)
}

# Function to generate a ggplot ROC curve for a given dispersion
generate_roc_plot <- function(disp_value, panel_label) {
  n_runs <- 40
  fpr_seq <- seq(0, 0.1, by = 0.002)
  
  tprs_exact <- matrix(0, nrow = n_runs, ncol = length(fpr_seq))
  tprs_wald <- matrix(0, nrow = n_runs, ncol = length(fpr_seq))
  tprs_lu <- matrix(0, nrow = n_runs, ncol = length(fpr_seq))
  
  for (i in 1:n_runs) {
    sim <- simulate_counts_adjusted(disp_value = disp_value)
    counts <- sim$counts
    groups <- sim$groups
    truth <- sim$gene_types
    
    scores_exact <- edgeR_scores(counts, groups, "Exact.EB")
    scores_wald <- edgeR_scores(counts, groups, "Wald.EB")
    scores_lu <- runLu(counts, groups)
    
    roc_exact <- suppressWarnings(roc(truth, scores_exact, direction = "<"))
    roc_wald <- suppressWarnings(roc(truth, scores_wald, direction = "<"))
    roc_lu <- suppressWarnings(roc(truth, scores_lu, direction = "<"))
    
    tprs_exact[i, ] <- approx(1 - roc_exact$specificities, roc_exact$sensitivities, 
                              xout = fpr_seq, yleft = 0, yright = 1)$y
    tprs_wald[i, ] <- approx(1 - roc_wald$specificities, roc_wald$sensitivities, 
                             xout = fpr_seq, yleft = 0, yright = 1)$y
    tprs_lu[i, ] <- approx(1 - roc_lu$specificities, roc_lu$sensitivities, 
                           xout = fpr_seq, yleft = 0, yright = 1)$y
  }
  
  avg_exact <- colMeans(tprs_exact)
  avg_wald <- colMeans(tprs_wald)
  avg_lu <- colMeans(tprs_lu)
  
  df <- data.frame(
    FPR = rep(fpr_seq, 3),
    TPR = c(avg_exact, avg_lu, avg_wald),
    Method = factor(rep(c("Exact.EB", "Wald.PL", "Wald.EB"), each = length(fpr_seq)),
                    levels = c("Exact.EB", "Wald.PL", "Wald.EB"))
  )
  
  p <- ggplot(df, aes(x = FPR, y = TPR, color = Method, linetype = Method)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("blue", "green", "red")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    labs(title = paste0(panel_label, " phi=", disp_value),
         x = "FPR", y = "TPR") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# Generate the three plots
set.seed(123) # for reproducibility

plotA <- generate_roc_plot(0.17, "A")
plotB <- generate_roc_plot(0.42, "B")
plotC <- generate_roc_plot(0.95, "C")

# Arrange them side by side
grid.arrange(plotA, plotB, plotC, ncol = 3)
