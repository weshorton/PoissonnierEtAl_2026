####################################################################
########################## FIGURE 4 WATERFALL ######################
####################################################################

### Make waterfall plot and run stats for main figure 4 panel C

source("./code/figure4_sourceRef.R")
library(data.table)
library(ggplot2)
library(ggpubr)

###
### Set Up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Read in Data
sheetNames_v <- c("d25 aPD-1 groups", "d25 aPD-L1 groups")
file_v <- "./data/fig4/tumorBurdenData.xlsx"
data_lsdt <- lapply(sheetNames_v, function(x) setDT(readxl::read_excel(file_v, sheet = x)))
names(data_lsdt) <- sheetNames_v

### Combine PDL1 combo
apdl1_v <- c("aPD-L1", "AIN_PTX_aPD-L1", "PLX_PTX_aPD-L1")
data_lsdt$allTreats <- rbind(data_lsdt$`d25 aPD-1 groups`,
                             data_lsdt$`d25 aPD-L1 groups`[Treatment %in% apdl1_v,])

### Change aCSF1R to PLX
data_lsdt <- sapply(data_lsdt, function(x) {
  x$Short <- gsub("CSF1Ri", "PLX", x$Short)
  return(x)}, simplify = F, USE.NAMES = T)

### Add "D" to each column name
metaCols_v <- c("Treatment", "Short", "Sample")
data_lsdt <- sapply(data_lsdt, function(x) {
  for (col_v in setdiff(colnames(x), metaCols_v)) {
    colnames(x)[colnames(x) == col_v] <- paste0("D", col_v)
  } # for
  return(x)}, simplify = F, USE.NAMES = T)

###
### Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

data_lsdt <- sapply(data_lsdt, function(x) {
  y <- calcChange(wide_dt = x, start_v = "D8", end_v = "D18")
  return(y)}, simplify = F, USE.NAMES = T)

data_lsdt <- sapply(data_lsdt, function(x) {
  x <- classifyChange(data_dt = x, cols_v = NULL, which_v = "pct")
  return(x)}, simplify = F, USE.NAMES = T)

###
### R/NR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### What is rate of change for aPD-1 R/NR between D8 and D18?
aPD1_check_dt <- as.data.table(table(data_lsdt$`d25 aPD-1 groups`[Short == "PLX+PTX+aPD-1", "D8_to_D18_pctRate"]))

### What about aPD-L1
aPDL1_check_dt <- as.data.table(table(data_lsdt$`d25 aPD-L1 groups`[Short == "PLX+PTX+aPD-L1", "D8_to_D18_pctRate"]))

### Get responders
aPDL1_responders_dt <- data_lsdt$`d25 aPD-L1 groups`[Short == "PLX+PTX+aPD-L1" & D8_to_D18_pctRate == "decrease: <0%", ]
aPD1_responders_dt <- data_lsdt$`d25 aPD-L1 groups`[Short == "PLX+PTX+aPD-1" & D8_to_D18_pctRate == "decrease: <0%", ]

###
### Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

plot_lsgg <- sapply(data_lsdt, function(x) {
  tmpColors_v <- colors_v[intersect(names(colors_v), unique(x$Short))]
  p_gg <- tumorWaterfall(wide_dt = x, xVar_v = "Sample", treatCol_v = "Short",
                         treatOrder_v = names(tmpColors_v), treatColors_v = tmpColors_v, addCounts_v = F,
                         yVar_v = "D8_to_D18_pctChange", groupVar_v = "D8_to_D18_pctRate",
                         addGaps_v = T, title_v = "D8 to D18 Tumor Burden Change", thinBars_v = T)
  return(p_gg$annotPlot)
}, simplify = F, USE.NAMES = T)

pdf(file = "./fig4/fig4c.pdf", width = 14, height = 10)
print(plot_lsgg$allTreats)
dev.off()

###
### Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

stats_KW_lsdt <- stats_Wilcox_lsdt <- list()
for (i in 1:length(data_lsdt)) {
  
  ### Get info
  currName_v <- names(data_lsdt)[i]
  currData_dt <- data_lsdt[[currName_v]]
  
  ### Make sure factor
  currData_dt$Short <- factor(currData_dt$Short, 
                              levels = intersect(names(colors_v), unique(currData_dt$Short)))
  
  ### Run Stats
  currStatsKW_dt <- simpleStats(data_dt = currData_dt, depVar_v = "D8_to_D18_pctChange",
                                indVar_v = "Short", indOrder_v = NULL, compType_v = "KW")
  currStatsWilcox_dt <- simpleStats(data_dt = currData_dt, depVar_v = "D8_to_D18_pctChange",
                                    indVar_v = "Short", indOrder_v = NULL, compType_v = "wilcox")
  
  stats_KW_lsdt[[currName_v]] <- currStatsKW_dt
  stats_Wilcox_lsdt[[currName_v]] <- currStatsWilcox_dt
  
} # for i

writexl::write_xlsx(x = stats_Wilcox_lsdt["allTreats"], path = "./fig4/wilcoxResults.xslx")
