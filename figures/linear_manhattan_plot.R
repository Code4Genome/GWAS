library(qqman)
head (results)

results$Chr = as.numeric(result$Chr)

max_y =30

par(las = 2)

png("manhattanplot.png", width=12, height=12, units="cm", res=300, type="cairo")

manhattan(results, chr = "Chr.y", bp = "Position", p = "pvalue", snp = "gene",
          col = c("gray10", "gray60"), chrlabs = NULL,
          suggestiveline = -log10(1e-04), genomewideline = -log10(5e-08),
          highlight = NULL, logp = TRUE, annotatePval = 1e-04,
          annotateTop = TRUE, ylim = c(0, max_y))

dev.off()
