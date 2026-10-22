# A cleaned up R file that builds on age_analysis.r to explore the relationship 
# between the proportion of A2 clade sequences and the proportion of 
#adult cases (>=18y) by year.

## Produces:
# 1. 

### First run age_analysis.r so you have 'allmeta' in your environment


# Store so we can manipulate, and take only data from 2013
d <- allmeta[allmeta$year >= 2013, ]

# proportion A2 per year (denominator: sequences with a clade call)
clade_d <- d[!is.na(d$subgenogroup), ]
a2_by_year <- aggregate(subgenogroup ~ year, data = clade_d,
                         FUN = function(x) mean(x == "A2"))
names(a2_by_year)[2] <- "prop_a2"
n_clade <- aggregate(subgenogroup ~ year, data = clade_d, FUN = length)
names(n_clade)[2] <- "n_clade"
a2_by_year <- merge(a2_by_year, n_clade, by = "year")

# proportion >=18y per year (denominator: sequences with age data)
age_d <- d[d$age_range3 %in% c("<18y", ">=18y"), ]
age_by_year <- aggregate(age_range3 ~ year, data = age_d,
                          FUN = function(x) mean(x == ">=18y"))
names(age_by_year)[2] <- "prop_adult"
n_age <- aggregate(age_range3 ~ year, data = age_d, FUN = length)
names(n_age)[2] <- "n_age"
age_by_year <- merge(age_by_year, n_age, by = "year")

# inner merge drops years with no age data (n_age == 0) automatically
plot_df <- merge(a2_by_year, age_by_year, by = "year")

cor_test <- cor.test(plot_df$prop_a2, plot_df$prop_adult, method = "spearman")

# scale the point size n_age (number of seqs)
cex_vals <- 1 + 3 * (plot_df$n_age - min(plot_df$n_age)) /
                    (max(plot_df$n_age) - min(plot_df$n_age))

plot(plot_df$prop_a2, plot_df$prop_adult,
     cex = cex_vals, pch = 19, col = "#D55E00",
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Proportion A2", ylab = "Proportion \u226518y",
     main = "A2 clade dominance vs. adult case proportion, by year")

abline(lm(prop_adult ~ prop_a2, data = plot_df, weights = n_age),
       lty = 2, col = "grey50")

text(plot_df$prop_a2, plot_df$prop_adult, labels = plot_df$year,
     pos = 3, offset = 0.6, cex = 0.8)

legend("bottomright",
       legend = sprintf("Spearman r = %.2f, p = %.3f",
                         cor_test$estimate, cor_test$p.value),
       bty = "n")

### Better plot -- gives a better legend and outlines circles
# --- reference sizes for the legend (pick round-ish numbers spanning the data) ---
size_breaks <- pretty(plot_df$n_age, n = 3)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max(plot_df$n_age)]

cex_from_n <- function(n) 0.6 + 1.8 * (n - min(plot_df$n_age)) /
                                      (max(plot_df$n_age) - min(plot_df$n_age))

png("a2_vs_adult.png", width = 4.2, height = 3.7, units = "in", res = 400)

par(mar = c(3.5, 3.5, 2, 1), family = "sans", mgp = c(2, 0.6, 0))

# 1. empty plot to set up axes/coordinate space only
plot(plot_df$prop_a2, plot_df$prop_adult, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Proportion variant A2", ylab = "Proportion adult cases", #\u226518y",
     main = "Proportion of A2 variant vs. adult cases, by year",
     bty = "l", las = 1, cex.lab = 0.85, cex.axis = 0.7, cex.main = 0.85, font.main = 2)

# 2. gridlines behind the data
grid(col = "grey90", lty = 1, lwd = 0.6)
box(bty = "l")

# 3. trend line
abline(lm(prop_adult ~ prop_a2, data = plot_df, weights = n_age),
       lty = 2, col = "grey50", lwd = 1)

# 4. points on top of everything
points(plot_df$prop_a2, plot_df$prop_adult,
       cex = cex_from_n(plot_df$n_age), pch = 21,
       bg = adjustcolor("#D55E00", alpha.f = 0.85), col = "grey30", lwd = 0.6)

text(plot_df$prop_a2, plot_df$prop_adult, labels = plot_df$year,
     pos = 3, offset = 0.5, cex = 0.55, col = "grey20")

legend("bottomright",
       legend = sprintf("Spearman r = %.2f, p = %.3f",
                         cor_test$estimate, cor_test$p.value),
       bty = "n", cex = 0.65)

legend("topleft",
       legend = size_breaks,
       pch = 21, pt.bg = adjustcolor("#D55E00", alpha.f = 0.85), col = "grey30",
       pt.cex = cex_from_n(size_breaks),
       title = "n (age data)", bty = "o", box.col = "grey60", bg = "white",
       cex = 0.6, y.intersp = 1.5)

dev.off()