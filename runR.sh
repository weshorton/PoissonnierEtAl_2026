#!/bin/sh

# Set to "/" for CodeOcean and "./" for github
baseDir="./"

mkdir -p $baseDir/results/fig2 $baseDir/results/fig4 $baseDir/results/fig7

### Figure 2 scripts
Rscript $baseDir/code/figure2_survivalGroups.R
Rscript $baseDir/code/figure2_survival.R
Rscript $baseDir/code/figure2_box.R
Rscript $baseDir/code/figure2_line.R

### Figure 4 scripts
Rscript $baseDir/code/figure4.R

### Figure 7 scripts
Rscript $baseDir/code/figure7_wrangle.R
Rscript $baseDir/code/figure7.R

### Extracting legend from ggplot in figure 4 code creates
### a blank Rplots.pdf file
rm Rplots.pdf
