###########################################################
######################## FIGURE 7 #########################
###########################################################

### Make line plots for 
### main figure 7 (panels A, B, C, D, E, F, G, H, I, and J)
### supplemental figure 7 (panels A, B, C, D, E, F, G, and H)

### Some mainuscript figures were made in prism, with the same data

source("./code/figure7_sourceRef.R")
library(data.table)
library(ggplot2)
library(immunarch)

###
### Read in ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Dirs
inDir_v <- "./data/fig7/rds"
plotDir_v <- "./fig7"

### Data
diversityMetrics_lsdt <- readRDS(file.path(inDir_v, "diversityMetrics_lsdt.rds"))
freqGroups_lsdt <- readRDS(file.path(inDir_v, "freqGroups_lsdt.rds"))
topClones_lsdt <- readRDS(file.path(inDir_v, "topClones_lsdt.rds"))
jaccard_lsdt <- readRDS(file.path(inDir_v, "jaccard_lsdt.rds"))
trackCloneCorpus_lsdt <- readRDS(file.path(inDir_v, "trackCloneCorpus_lsdt.rds"))
trackCloneQuery_lsdt <- readRDS(file.path(inDir_v, "trackCloneQuery_lsdt.rds"))
plotInfo_lsv <- readRDS(file.path(inDir_v, "plotInfo_lsv.rds"))

### Output
fig7_lsdt <- list()
fig7Titles_lsv <- list()
suppFig7_lsdt <- list()
suppFig7Titles_lsv <- list()

###
### Main Figure Panels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

###
### 7A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
fig7A_dt <- diversityMetrics_lsdt$Tumor[Experiment != "d95", mget(c("geoSample", "Treatment", "Title", "Clonality"))]

### Labels
fig7A_dt[, yPos := min(Clonality)]
fig7A_counts <- as.data.table(table(fig7A_dt$Treatment))
fig7A_dt <- merge(fig7A_dt, fig7A_counts, by.x = "Treatment", by.y = "V1", sort = F)

### Title
fig7ATitle_v <- "Fig 7a - Treatwise Clonal Index of end-stage (d105) tumors\nas well as d80 WT and Untreated"

### Plot
fig7A_gg <- ggplot(data = fig7A_dt, aes(x = Title, y = Clonality)) +
  geom_violin(aes(color = Treatment)) +
  geom_jitter(shape = 15, size = 2, position = position_jitter(width = 0.05, height = 0)) +
  geom_text(aes(y = yPos, label = N), size = 3, vjust = 2) +
  scale_color_manual(values = plotInfo_lsv$treatColors, breaks = names(plotInfo_lsv$treatColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  guides(color = "none") +
  ggtitle(fig7ATitle_v)

### Output
fig7_lsdt[["7a"]] <- fig7A_dt
fig7Titles_lsv[["7a"]] <- gsub("\\\n", ";", fig7ATitle_v)

pdf(file = file.path(plotDir_v, "7a.pdf"), height = 6, width = 8)
fig7A_gg
dev.off()


###
### 7B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
fig7B_dt <- diversityMetrics_lsdt$Tumor[Experiment != "d95", mget(c("geoSample", "Treatment", "Title", "Shannon Entropy"))]

### Labels
fig7B_dt[, yPos := min(`Shannon Entropy`)]
fig7B_counts <- as.data.table(table(fig7B_dt$Treatment))
fig7B_dt <- merge(fig7B_dt, fig7B_counts, by.x = "Treatment", by.y = "V1", sort = F)

### Title
fig7BTitle_v <- "Fig 7b - Treatwise Shannon Entropy of end-stage (d105) tumors\nas well as d80 WT and Untreated"

### Plot
fig7B_gg <- ggplot(data = fig7B_dt, aes(x = Title, y = `Shannon Entropy`)) +
  geom_violin(aes(color = Treatment)) +
  geom_jitter(shape = 15, size = 2, position = position_jitter(width = 0.05, height = 0)) +
  geom_text(aes(y = yPos, label = N), size = 3, vjust = 2) +
  scale_color_manual(values = plotInfo_lsv$treatColors, breaks = names(plotInfo_lsv$treatColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  guides(color = "none") +
  ggtitle(fig7BTitle_v)

### Output
fig7_lsdt[["7b"]] <- fig7B_dt
fig7Titles_lsv[["7b"]] <- gsub("\\\n", ";", fig7BTitle_v)

pdf(file = file.path(plotDir_v, "7b.pdf"), height = 6, width = 8)
fig7B_gg
dev.off()

###
### 7C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
fig7C_dt <- topClones_lsdt$Tumor[Experiment != "d95",]
fig7C_dt$Sample <- NULL; fig7C_dt$Experiment <- NULL

### Title
fig7CTitle_v <- "Fig 7c - Binned Clonal Frequencies of Primary Tumors"

### Prep Plot
tmpFig7C_dt <- melt(copy(fig7C_dt), measure.vars = plotInfo_lsv$topClones)
tmpFig7C_dt$variable <- factor(tmpFig7C_dt$variable, levels = rev(plotInfo_lsv$topClones))
tmpFig7C_dt[, mean := mean(value, na.rm = T), by = c("Title", "variable")]
tmpFig7C_dt <- unique(tmpFig7C_dt[,mget(c("Title", "Treatment", "mean", "variable"))])

### Labels
tmpFig7C_dt <- merge(tmpFig7C_dt, unique(fig7B_dt[,mget(c("Treatment", "N"))]), by = "Treatment", sort = F)

### Plot
fig7C_gg <- ggplot(data = tmpFig7C_dt, aes(x = Title, y = mean, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(y = 0, label = N), size = 3, vjust = 2) +
  scale_fill_manual(values = plotInfo_lsv$topCloneColors, breaks = names(plotInfo_lsv$topCloneColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  ggtitle(fig7CTitle_v)

### Output
fig7_lsdt[["7c"]] <- fig7C_dt
fig7Titles_lsv[["7c"]] <- gsub("\\\n", ";", fig7CTitle_v)

pdf(file = file.path(plotDir_v, "7c.pdf"), height = 6, width = 8)
fig7C_gg
dev.off()

###
### 7D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
fig7D_dt <- freqGroups_lsdt$Tumor[Experiment != "d95",]
fig7D_dt$Sample <- NULL; fig7D_dt$Experiment <- NULL

### Title
fig7DTitle_v <- "Fig 7d - Binned Clonal Frequencies of Primary Tumors"

### Prep Plot
tmpFig7D_dt <- melt(copy(fig7D_dt), measure.vars = plotInfo_lsv$freqGroups)
tmpFig7D_dt$variable <- factor(tmpFig7D_dt$variable, levels = plotInfo_lsv$freqGroups)
tmpFig7D_dt[, mean := mean(value, na.rm = T), by = c("Title", "variable")]
tmpFig7D_dt <- unique(tmpFig7D_dt[,mget(c("Title", "Treatment", "mean", "variable"))])

### Labels
tmpFig7D_dt <- merge(tmpFig7D_dt, unique(fig7B_dt[,mget(c("Treatment", "N"))]), by = "Treatment", sort = F)

### Plot
fig7D_gg <- ggplot(data = tmpFig7D_dt, aes(x = Title, y = mean, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(y = 0, label = N), size = 3, vjust = 2) +
  scale_fill_manual(values = plotInfo_lsv$freqGroupColors, breaks = names(plotInfo_lsv$freqGroupColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  ggtitle(fig7DTitle_v)

### Output
fig7_lsdt[["7d"]] <- fig7D_dt
fig7Titles_lsv[["7d"]] <- gsub("\\\n", ";", fig7DTitle_v)

pdf(file = file.path(plotDir_v, "7d.pdf"), height = 6, width = 8)
fig7D_gg
dev.off()

###
### 7E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
fig7E_dt <- jaccard_lsdt$TumorTop100

### Title
fig7ETitle_v <- "Fig 7e - Jaccard Index Between Top 100 Clones of Primary Tumor Samples"

### Plot
fig7E_gg <- ggplot(fig7E_dt, aes(x = treatmentComparison, y = Jaccard)) +
  geom_violin(aes(color = treatmentComparison)) +
  geom_jitter(shape = 15, size = 0.5, position = position_jitter(width = 0.05, height = 0)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  guides(color = "none") +
  ggtitle(fig7ETitle_v)

### Output
fig7_lsdt[["7e"]] <- fig7E_dt
fig7Titles_lsv[["7e"]] <- gsub("\\\n", ";", fig7ETitle_v)

pdf(file = file.path(plotDir_v, "7e.pdf"), height = 6, width = 8)
fig7E_gg
dev.off()

###
### 7F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
fig7F_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$TumorTop50, .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")

### Title
fig7FTitle_v <- "Fig 7f - Top 50 3x aPD-1R clones in primary tumors in top 50 of other 3x treatments"
fig7FaltTitle_v <- "Fig 7f - Top 50 3x aPD-1R clones in primary tumors in all clones of other 3x treatments"

### Plot
fig7F_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$TumorTop50, 
                                                 .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")))
fig7F_gg <- fig7F_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7FTitle_v)

fig7Falt_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$TumorAll,
                                                    .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")))
fig7Falt_gg <- fig7Falt_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7FaltTitle_v)

### Output
fig7_lsdt[["7f"]] <- fig7F_dt
fig7Titles_lsv[["7f"]] <- fig7FTitle_v

pdf(file = file.path(plotDir_v, "7f.pdf"), height = 6, width = 8)
fig7F_gg
fig7Falt_gg
dev.off()

###
### 7G ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
fig7G_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$TumorTop50, .which = trackCloneQuery_lsdt$Tumor3xaPDL1, .col = "aa")

### Title
fig7GTitle_v <- "Fig 7g - Top 50 3x aPD-L1 clones in primary tumors in top 50 of other 3x treatments"
fig7GaltTitle_v <- "Fig 7g - Top 50 3x aPD-L1 clones in primary tumors in all clones of other 3x treatments"

### Plot
fig7G_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$TumorTop50, 
                                                 .which = trackCloneQuery_lsdt$Tumor3xaPDL1, .col = "aa")))
fig7G_gg <- fig7G_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7GTitle_v)

fig7Galt_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$TumorAll, 
                                                    .which = trackCloneQuery_lsdt$Tumor3xaPDL1, .col = "aa")))
fig7Galt_gg <- fig7Galt_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7GaltTitle_v)

### Output
fig7_lsdt[["7g"]] <- fig7G_dt
fig7Titles_lsv[["7g"]] <- fig7GTitle_v

pdf(file = file.path(plotDir_v, "7g.pdf"), height = 6, width = 8)
fig7G_gg
fig7Galt_gg
dev.off()

###
### 7H ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
fig7H_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")

### Title
fig7HTitle_v <- "Fig 7h - Top 50 3x aPD-1R clones in primary tumors in top 50 of other 3x treatments in lung mets"

### Plot
fig7H_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, 
                                                 .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")))
fig7H_gg <- fig7H_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7HTitle_v)

### Output
fig7_lsdt[["7h"]] <- fig7H_dt
fig7Titles_lsv[["7h"]] <- fig7HTitle_v

pdf(file = file.path(plotDir_v, "7h.pdf"), height = 6, width = 8)
fig7H_gg
dev.off()

###
### 7I left ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
fig7Il_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$BloodAll, .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")

### Title
fig7IlTitle_v <- "Fig 7i left - Top 50 3x aPD-1R clones in primary tumors in all clones of other 3x treatments in blood"

### Plot
fig7Il_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$BloodAll, 
                                                  .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")))
fig7Il_gg <- fig7Il_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7IlTitle_v)

### Output
fig7_lsdt[["7il"]] <- fig7Il_dt
fig7Titles_lsv[["7il"]] <- fig7IlTitle_v

pdf(file = file.path(plotDir_v, "7il.pdf"), height = 6, width = 8)
fig7Il_gg
dev.off()

###
### 7I right ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
fig7Ir_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$BloodTop50, .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")

### Title
fig7IrTitle_v <- "Fig 7i right - Top 50 3x aPD-1R clones in primary tumors in top 50 of other 3x treatments in blood"

### Plot
fig7Ir_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$BloodTop50, 
                                                  .which = trackCloneQuery_lsdt$Tumor3xR, .col = "aa")))
fig7Ir_gg <- fig7Ir_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7IrTitle_v)

### Output
fig7_lsdt[["7ir"]] <- fig7Ir_dt
fig7Titles_lsv[["7ir"]] <- fig7IrTitle_v

pdf(file = file.path(plotDir_v, "7ir.pdf"), height = 6, width = 8)
fig7Ir_gg
dev.off()

###
### 7J ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
fig7J_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, .which = trackCloneQuery_lsdt$Tumor3xaPDL1, .col = "aa")

### Title
fig7JTitle_v <- "Fig 7h - Top 50 3x aPD-L1 clones in primary tumors in top 50 of other 3x treatments in lung mets"

### Plot
fig7J_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, 
                                                 .which = trackCloneQuery_lsdt$Tumor3xaPDL1, .col = "aa")))
fig7J_gg <- fig7J_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(fig7JTitle_v)

### Output
fig7_lsdt[["7j"]] <- fig7J_dt
fig7Titles_lsv[["7j"]] <- fig7JTitle_v

pdf(file = file.path(plotDir_v, "7j.pdf"), height = 6, width = 8)
fig7J_gg
dev.off()

###
### Supplemental Figure Panels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

###
### S7A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
figS7A_dt <- jaccard_lsdt$Tumor

### Labels
figS7A_counts <- as.data.table(table(figS7A_dt$treatmentComparison))
figS7A_dt <- merge(figS7A_dt, figS7A_counts, by.x = "treatmentComparison", by.y = "V1", sort = F)

### Title
figS7ATitle_v <- "Fig S7a - Jaccard Index Between Primary Tumor Samples"

### Plot
figS7A_gg <- ggplot(figS7A_dt, aes(x = treatmentComparison, y = Jaccard)) +
  geom_violin(aes(color = treatmentComparison)) +
  geom_jitter(shape = 15, size = 0.5, position = position_jitter(width = 0.05, height = 0)) +
  geom_text(aes(y = 0, label = N), size = 3, vjust = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  guides(color = "none") +
  ggtitle(figS7ATitle_v)

### Output
suppFig7_lsdt[["S7a"]] <- figS7A_dt
suppFig7Titles_lsv[["S7a"]] <- gsub("\\\n", ";", figS7ATitle_v)

pdf(file = file.path(plotDir_v, "S7a.pdf"), height = 6, width = 8)
figS7A_gg
dev.off()

###
### S7B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
figS7B_dt <- diversityMetrics_lsdt$Lung[Experiment != "d95", mget(c("geoSample", "Treatment", "Title", "Shannon Entropy"))]

### Labels
figS7B_dt[, yPos := min(`Shannon Entropy`)]
figS7B_counts <- as.data.table(table(figS7B_dt$Treatment))
figS7B_dt <- merge(figS7B_dt, figS7B_counts, by.x = "Treatment", by.y = "V1", sort = F)

### Title
figS7BTitle_v <- "Fig S7b - Treatwise Shannon Entropy of end-stage (d105) lung mets\nas well as d80 WT and Untreated"

### Plot
figS7B_gg <- ggplot(data = figS7B_dt, aes(x = Title, y = `Shannon Entropy`)) +
  geom_violin(aes(color = Treatment)) +
  geom_jitter(shape = 15, size = 2, position = position_jitter(width = 0.05, height = 0)) +
  geom_text(aes(y = yPos, label = N), size = 3, vjust = 2) +
  scale_color_manual(values = plotInfo_lsv$treatColors, breaks = names(plotInfo_lsv$treatColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  guides(color = "none") +
  ggtitle(figS7BTitle_v)

### Output
suppFig7_lsdt[["S7b"]] <- figS7B_dt
suppFig7Titles_lsv[["7b"]] <- gsub("\\\n", ";", figS7BTitle_v)

pdf(file = file.path(plotDir_v, "S7b.pdf"), height = 6, width = 8)
figS7B_gg
dev.off()

###
### S7C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
figS7C_dt <- diversityMetrics_lsdt$Lung[Experiment != "d95", mget(c("geoSample", "Treatment", "Title", "Clonality"))]

### Labels
figS7C_dt[, yPos := min(Clonality)]
figS7C_counts <- as.data.table(table(figS7C_dt$Treatment))
figS7C_dt <- merge(figS7C_dt, figS7C_counts, by.x = "Treatment", by.y = "V1", sort = F)

### Title
figS7CTitle_v <- "Fig S7c - Treatwise Clonal Index of end-stage (d105) lung mets\nas well as d80 WT and Untreated"

### Plot
figS7C_gg <- ggplot(data = figS7C_dt, aes(x = Title, y = Clonality)) +
  geom_violin(aes(color = Treatment)) +
  geom_jitter(shape = 15, size = 2, position = position_jitter(width = 0.05, height = 0)) +
  geom_text(aes(y = yPos, label = N), size = 3, vjust = 2) +
  scale_color_manual(values = plotInfo_lsv$treatColors, breaks = names(plotInfo_lsv$treatColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  guides(color = "none") +
  ggtitle(figS7CTitle_v)

### Output
suppFig7_lsdt[["S7b"]] <- figS7C_dt
suppFig7Titles_lsv[["S7b"]] <- figS7CTitle_v

pdf(file = file.path(plotDir_v, "S7c.pdf"), height = 6, width = 8)
figS7C_gg
dev.off()

###
### S7D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
figS7D_dt <- topClones_lsdt$Lung[Experiment != "d95",]
figS7D_dt$Sample <- NULL; figS7D_dt$Experiment <- NULL

### Title
figS7DTitle_v <- "Fig S7d - Binned Clonal Frequencies of Lung Mets"

### Prep Plot
tmpFigS7D_dt <- melt(copy(figS7D_dt), measure.vars = plotInfo_lsv$topClones)
tmpFigS7D_dt$variable <- factor(tmpFigS7D_dt$variable, levels = rev(plotInfo_lsv$topClones))
tmpFigS7D_dt[, mean := mean(value, na.rm = T), by = c("Title", "variable")]
tmpFigS7D_dt <- unique(tmpFigS7D_dt[,mget(c("Title", "Treatment", "mean", "variable"))])

### Labels
tmpFigS7D_dt <- merge(tmpFigS7D_dt, unique(figS7C_dt[,mget(c("Treatment", "N"))]), by = "Treatment", sort = F)

### Plot
figS7D_gg <- ggplot(data = tmpFigS7D_dt, aes(x = Title, y = mean, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(y = 0, label = N), size = 3, vjust = 2) +
  scale_fill_manual(values = plotInfo_lsv$topCloneColors, breaks = names(plotInfo_lsv$topCloneColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  ggtitle(figS7DTitle_v)

### Output
suppFig7_lsdt[["S7d"]] <- figS7D_dt
suppFig7Titles_lsv[["S7d"]] <- gsub("\\\n", ";", figS7DTitle_v)

pdf(file = file.path(plotDir_v, "S7d.pdf"), height = 6, width = 8)
figS7D_gg
dev.off()

###
### S7E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
figS7E_dt <- freqGroups_lsdt$Tumor[Experiment != "d95",]
figS7E_dt$Sample <- NULL; figS7E_dt$Experiment <- NULL

### Title
figS7ETitle_v <- "Fig S7e - Binned Clonal Frequencies of Lung Meets"

### Prep Plot
tmpFigS7E_dt <- melt(copy(figS7E_dt), measure.vars = plotInfo_lsv$freqGroups)
tmpFigS7E_dt$variable <- factor(tmpFigS7E_dt$variable, levels = plotInfo_lsv$freqGroups)
tmpFigS7E_dt[, mean := mean(value, na.rm = T), by = c("Title", "variable")]
tmpFigS7E_dt <- unique(tmpFigS7E_dt[,mget(c("Title", "Treatment", "mean", "variable"))])

### Labels
tmpFigS7E_dt <- merge(tmpFigS7E_dt, unique(figS7C_dt[,mget(c("Treatment", "N"))]), by = "Treatment", sort = F)

### Plot
figS7E_gg <- ggplot(data = tmpFigS7E_dt, aes(x = Title, y = mean, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(y = 0, label = N), size = 3, vjust = 2) +
  scale_fill_manual(values = plotInfo_lsv$freqGroupColors, breaks = names(plotInfo_lsv$freqGroupColors)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  ggtitle(figS7ETitle_v)

### Output
suppFig7_lsdt[["S7e"]] <- figS7E_dt
suppFig7Titles_lsv[["S7e"]] <- gsub("\\\n", ";", figS7ETitle_v)

pdf(file = file.path(plotDir_v, "S7e.pdf"), height = 6, width = 8)
figS7E_gg
dev.off()

###
### S7F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset
figS7F_dt <- jaccard_lsdt$Lung

### Labels
figS7F_counts <- as.data.table(table(figS7F_dt$treatmentComparison))
figS7F_dt <- merge(figS7F_dt, figS7F_counts, by.x = "treatmentComparison", by.y = "V1", sort = F)

### Title
figS7FTitle_v <- "Fig S7a - Jaccard Index Between Lung Met Samples"

### Plot
figS7F_gg <- ggplot(figS7F_dt, aes(x = treatmentComparison, y = Jaccard)) +
  geom_violin(aes(color = treatmentComparison)) +
  geom_jitter(shape = 15, size = 0.5, position = position_jitter(width = 0.05, height = 0)) +
  geom_text(aes(y = 0, label = N), size = 3, vjust = 2) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  guides(color = "none") +
  ggtitle(figS7FTitle_v)

### Output
suppFig7_lsdt[["S7f"]] <- figS7F_dt
suppFig7Titles_lsv[["S7f"]] <- gsub("\\\n", ";", figS7FTitle_v)

pdf(file = file.path(plotDir_v, "S7f.pdf"), height = 6, width = 8)
figS7F_gg
dev.off()

###
### S7G ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
figS7G_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, .which = trackCloneQuery_lsdt$Lung3xR, .col = "aa")

### Title
figS7GTitle_v <- "Fig S7g - Top 50 3x aPD-1R clones in lung mets in top 50 of other 3x treatments in lung mets"

### Plot
figS7G_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, 
                                                  .which = trackCloneQuery_lsdt$Lung3xR, .col = "aa")))
figS7G_gg <- figS7G_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(figS7GTitle_v)

### Output
suppFig7_lsdt[["S7g"]] <- figS7G_dt
suppFig7Titles_lsv[["S7g"]] <- figS7GTitle_v

pdf(file = file.path(plotDir_v, "S7g.pdf"), height = 6, width = 8)
figS7G_gg
dev.off()

###
### S7H ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Get data
figS7H_dt <- trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, .which = trackCloneQuery_lsdt$Lung3xaPDL1, .col = "aa")

### Title
figS7HTitle_v <- "Fig S7h - Top 50 3x aPD-1R clones in lung mets in top 50 of other 3x treatments in lung mets"

### Plot
figS7H_gg <- suppressWarnings(vis(trackClonotypes(.data = trackCloneCorpus_lsdt$LungTop50, 
                                                  .which = trackCloneQuery_lsdt$Lung3xaPDL1, .col = "aa")))
figS7H_gg <- figS7H_gg + scale_fill_manual(values = plotInfo_lsv$trackCloneColors, breaks = names(plotInfo_lsv$trackCloneColors)) +
  ggtitle(figS7HTitle_v)

### Output
suppFig7_lsdt[["7h"]] <- figS7H_dt
suppFig7Titles_lsv[["S7h"]] <- figS7HTitle_v

pdf(file = file.path(plotDir_v, "S7h.pdf"), height = 6, width = 8)
figS7H_gg
dev.off()

###
### Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Add titles
fig7_lsdt[["titles"]] <- as.data.frame(fig7Titles_lsv)
suppFig7_lsdt[["titles"]] <- as.data.frame(suppFig7Titles_lsv)

### Write - have to make this without muy package...
fig7_wb <- openxlsx::createWorkbook()
invisible(lapply(seq_along(fig7_lsdt), function(x) {
  openxlsx::addWorksheet(wb = fig7_wb, sheetName = names(fig7_lsdt)[x])
  openxlsx::writeData(wb = fig7_wb, sheet = x, x = fig7_lsdt[[x]])
}))
openxlsx::saveWorkbook(wb = fig7_wb, file = "./data/fig7/mainFigData.xlsx", overwrite = T)

suppFig7_wb <- openxlsx::createWorkbook()
invisible(lapply(seq_along(suppFig7_lsdt), function(x) {
  openxlsx::addWorksheet(wb = suppFig7_wb, sheetName = names(suppFig7_lsdt)[x])
  openxlsx::writeData(wb = suppFig7_wb, sheet = x, x = suppFig7_lsdt[[x]])
}))
openxlsx::saveWorkbook(wb = suppFig7_wb, file = "./data/fig7/suppFigData.xlsx", overwrite = T)