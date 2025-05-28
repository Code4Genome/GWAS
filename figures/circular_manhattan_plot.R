library(CMplot)
library(ggplot2)
library(tidyverse)
library(devtools)

set.seed(123)

data = read.xlsx("circular_MP.xlsx")
data(circular)
 
#Dataset containing snps, Chr, Position and pvalue

CMplot(data,
       type = "p",                    # Point plot
       plot.type = "c",               # Circular Manhattan plot
       col = c("grey30", "grey60"),        # Alternating colors for chromosomes
       chr.labels = paste("Chr", 1:22),  # Chromosome labels
       threshold = c(1e-4, 1e-2),      # Significance thresholds
       cir.chr.h=1.5,
       threshold.col = c("green", "orange"),  # Threshold line colors
       threshold.lty = c(1, 2),       # Threshold line types,   # Threshold line widths
       amplify = TRUE,                # Highlight significant SNPs
       signal.col = c("red","blue"),  # Color for significant SNPs
       chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,
       outward=FALSE,
       signal.cex = 1.5,              # Size of significant SNP points
       signal.line=1,
       main = "Circular Manhattan Plot",  # Title of the plot
       file = "jpg",                  # Output file format
       dpi = 300,                     # Resolution
       file.output = TRUE,            # Save output to a file
       verbose = TRUE,
       width=10,
       height=10,
       r= 0.4)  


