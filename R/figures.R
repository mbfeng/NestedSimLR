################################################################################
# figures.R
#
# Figures 1-3 of the paper, drawn with ggplot2 in the style of the original
# figures, from the summary files in a results folder.
################################################################################

# Figure 1(a): density of the outer scenarios S_tau.
# Figure 1(b): the conditional mean mu(theta). Closed form, no simulation.
figure1 <- function(out_dir) {
  library(ggplot2)
  params <- single_asset_params()
  probabilities <- seq(0.001, 0.999, by = 0.001)
  log_mean <- log(params$s0) + params$meanlog_out
  theta <- qlnorm(probabilities, log_mean, params$sdlog_out)
  curves <- data.frame(
    theta = theta,
    density = dlnorm(theta, log_mean, params$sdlog_out),
    mu = single_asset_true_value(theta, params)
  )
  x_label <- latex2exp::TeX("Projected stock prices $S_{\\tau}=\\theta$")

  density_plot <- ggplot(curves, aes(theta, density)) +
    geom_area(fill = "deeppink3", alpha = 0.5) +
    labs(x = x_label, y = "Outer scenario density function") +
    theme(axis.text = element_text(size = 16), text = element_text(size = 16))
  ggsave(file.path(out_dir, "Figure1a_OuterDist.pdf"), density_plot,
         width = 6, height = 4)

  mean_plot <- ggplot(curves, aes(theta, mu)) +
    geom_line(linewidth = 2, colour = "steelblue3") +
    labs(x = x_label, y = latex2exp::TeX("Conditional mean $\\mu(\\theta)$")) +
    theme(axis.text = element_text(size = 16), text = element_text(size = 16))
  ggsave(file.path(out_dir, "Figure1b_ConditionalMean.pdf"), mean_plot,
         width = 6, height = 4)
}

# Figure 2: 2.5% and 97.5% quantiles (over the macro runs) of the estimated
# mu(theta) at every scenario, for three designs, and the true mu(theta).
# Figure 2(b) zooms into 80 <= theta <= 120.
figure2 <- function(results_dir, out_dir) {
  library(ggplot2)
  bands <- read.csv(file.path(results_dir, "sec71_fig2_bands.csv"))
  bands$design <- factor(bands$design, levels = c("optimal", "SNS+", "regression"))
  legend_labels <- c("Opt. Design", "SNS+       ", "Regression ")

  band_plot <- ggplot(bands) +
    geom_ribbon(aes(x = theta, ymin = lower, ymax = upper,
                    colour = design, linetype = design),
                linewidth = 1.5, alpha = 0) +
    geom_line(aes(theta, true), linewidth = 1, colour = "black") +
    scale_linetype_manual(values = c(1, 2, 3), labels = legend_labels) +
    scale_colour_manual(values = c("#d7191c", "#7fc97f", "black"),
                        labels = legend_labels) +
    guides(colour = guide_legend(title = ""),
           linetype = guide_legend(title = "")) +
    theme_bw() +
    theme(legend.position = c(0.2, 0.8),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(1, "cm"),
          legend.key.width = unit(2, "cm"),
          text = element_text(size = 18)) +
    labs(x = latex2exp::TeX("outer scenarios $\\theta$"),
         y = latex2exp::TeX("conditional mean, $\\mu(\\theta)$"))
  ggsave(file.path(out_dir, "Figure2a_Bands.pdf"), band_plot,
         width = 7, height = 5)

  # expand = FALSE: the axes end exactly at the limits, as in the paper
  zoomed_plot <- band_plot +
    coord_cartesian(xlim = c(80, 120), ylim = c(30, 43), expand = FALSE)
  ggsave(file.path(out_dir, "Figure2b_Bands_Zoomed.pdf"), zoomed_plot,
         width = 7, height = 5)
}

# Figure 3(a): histogram of the variance ratios V(SNS) / V(optimal) at 1,000
# fixed scenarios (dashed line: their mean).
# Figure 3(b): the optimal design's budget against M (= N), averaged over the
# macro runs.
figure3 <- function(results_dir, out_dir) {
  library(ggplot2)
  ratios <- read.csv(file.path(results_dir, "sec73_fig3a_var_ratio.csv"))
  ratio_plot <- ggplot(ratios, aes(ratio)) +
    geom_histogram(bins = 50, colour = "black", fill = "red", alpha = 0.4) +
    geom_vline(xintercept = mean(ratios$ratio), colour = "blue",
               linetype = "dashed", linewidth = 1) +
    theme(text = element_text(size = 16, face = "bold")) +
    labs(x = "Variance Ratio", y = "Count")
  ggsave(file.path(out_dir, "Figure3a_VarRatio.pdf"), ratio_plot,
         width = 7, height = 5)

  budgets <- read.csv(file.path(results_dir, "sec73_fig3b_budgets.csv"))
  m_values <- sort(unique(budgets$M))
  average_budget <- numeric(length(m_values))
  for (i in seq_along(m_values)) {
    average_budget[i] <- mean(budgets$budget_optimal[budgets$M == m_values[i]])
  }
  growth <- data.frame(M = m_values, budget_optimal = average_budget)
  budget_plot <- ggplot(growth, aes(M, budget_optimal)) +
    geom_point(size = 5) +
    labs(x = "Number of Outer Scenarios, M", y = "",
         title = "Minimized Simulation Budget with M=N") +
    theme(plot.title = element_text(size = 14, face = "bold.italic",
                                    hjust = 0.5),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14))
  ggsave(file.path(out_dir, "Figure3b_SimBudgetGrowth.pdf"), budget_plot,
         width = 7, height = 5)
}
