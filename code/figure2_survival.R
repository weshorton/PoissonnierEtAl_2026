###################################################################
###################### FIGURE 2 SCATTER ###########################
###################################################################

### Set to "/" for CodeOcean and "./" for github
baseDir_v <- "./"

### Make line plots for 
  ### supplemental figure 2 panels S, T, U, and V

source(file.path(baseDir_v, "code/figure2_sourceRef.R"))
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

###
### Set Up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Read in data
data_dt <- as.data.table(readxl::read_xlsx(file.path(baseDir_v, "data/fig2/data.xlsx"), sheet = "all"))

### Read in survival
sheets_v <- readxl::excel_sheets(file.path(baseDir_v, "data/fig2/survivalGroups.xlsx"))
survGrp_lsdt <- lapply(sheets_v, function(x) {
  setDT(readxl::read_excel(path = file.path(baseDir_v, "data/fig2/survivalGroups.xlsx"), sheet = x))
})
names(survGrp_lsdt) <- sheets_v
survMed_dt <- setDT(readxl::read_excel(path = file.path(baseDir_v, "data/fig2/survivalMedians.xlsx"), sheet = "Sheet1"))

### Create metadata
metaCols_v <- c("PatientID", "CB", "BR", "progression (1 or 0)", "PFS_time", "TOT", "Dstart",
                "EOT Date", "Date_progression", "OFF Study", "EOT Reason")
meta_dt <- unique(data_dt[,mget(metaCols_v)])
#meta_dt <- unique(data_dt[,mget(c("PatientID", "BR", "CB"))])

### Remove non-evaluable patients
toRm_v <- meta_dt[BR == "NE", PatientID]
data_dt <- data_dt[!(PatientID %in% toRm_v),]

### Check
uniqPt_v <- unique(data_dt$PatientID)
if (length(uniqPt_v) != 35) stop("Bad patient substitution")

### Read in data type info
flow_dt <- fread(file.path(baseDir_v, "data/fig2/flowDataTypes.txt"))
luminex_dt <- fread(file.path(baseDir_v, "data/fig2/luminexDataTypes.txt"))
ref_dt <- rbind(flow_dt, luminex_dt)

### Change values to pct for those that require it
pct_dt <- ref_dt[Measure == "pct",]
for (i in 1:pct_dt[,.N]) {
  currPop_v <- pct_dt[i,Column]
  set(data_dt, j = currPop_v, value = data_dt[[currPop_v]]*100)
}

### Plot populations
kmPops_dt <- data.table("pop" = c("MyeloidDCs", "ClassicalMonocytes_CD14ppCD16n_PDL1"),
                        "measure" = c("S.C2D1", "C1D1.C2D1"),
                        "fig" = c("FigS2S", "FigS2U"))

scatterPops_dt <- data.table("pop" = c("MyeloidDCs", "ClassicalMonocytes_CD14ppCD16n_PDL1"),
                             "measure" = c("S.C2D1", "C1D1.C2D1"),
                             "fig" = c("FigS2T", "FigS2V"),
                             "xLab" = c("log2 Ratio of Percentage", 
                                       "log2 Ratio of PD-L1 Expression"),
                             "title" = c("Corr b/w L1-C2D1 Myeloid Freq and PFS",
                                         "Corr b/w C1D1-C2D1 Classical Monocyte PD-L1 and PFS"))

###
### Wrangle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Split data (a little unnecessary...)
data_lsdt <- komenFlowWrangle(data_dt = data_dt,
                              ref_dt = ref_dt,
                              grep_v = "CD8",
                              abs_v = "Screen",
                              num_v = c("C1D1", "C2D1", "C2D1"),
                              denom_v = c("Screen", "Screen", "C1D1"))

### Specify columns
metaCols_v <- colnames(data_lsdt[[1]])[1:(grep("CD8", colnames(data_lsdt[[1]]))[1]-1)]
meltCols_v <- c("PatientID", "Cycle", "CB", "progression (1 or 0)", "PFS_time")
measureCols_v <- setdiff(colnames(data_lsdt[[1]]), metaCols_v)
analytes_v <- luminex_dt$Column
pops_v <- setdiff(measureCols_v, analytes_v)
pctPops_v <- flow_dt[Measure == "pct",Column]
numPops_v <- flow_dt[Measure == "num",Column]

### Melt
melt_lsdt <- lapply(data_lsdt, function(x) {
  y <- melt(x[,mget(c(meltCols_v, measureCols_v))], id.vars = meltCols_v)
  y$PFS_time <- as.numeric(y$PFS_time)
  y$`progression (1 or 0)` <- as.numeric(y$`progression (1 or 0)`)
  return(y)
})

###
### Scatterplot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

for (i in 1:nrow(scatterPops_dt)) {
  
  ### Get info
  currPop_v <- scatterPops_dt$pop[i]
  currMeasure_v <- scatterPops_dt$measure[i]
  currFigName_v <- scatterPops_dt$fig[i]
  currXlab_v <- scatterPops_dt$xLab
  currTitle_v <- scatterPops_dt$title
  
  ### Subset Data
  currData_dt <- melt_lsdt[[currMeasure_v]][variable == currPop_v,]
  
  ### Subset NAs
  currData_dt <- currData_dt[!is.na(value),]
  
  ### Make plot
  curr_gg <- ggplot(data = currData_dt, aes(x = value, y = PFS_time)) +
    geom_point(size = 4) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 20), 
                            plot.subtitle = element_text(hjust = 0.5, size = 14), 
                            axis.text = element_text(size = 16), 
                            axis.title = element_text(size = 18), 
                            legend.text = element_text(size = 16), 
                            legend.title = element_text(size = 18)) +
    labs(x = currXlab_v, y = "PFS time (days)") +
    geom_smooth(method = 'lm', se = F) +
    stat_cor(method = "pearson") +
    ggtitle(currTitle_v)
  
  ### Cox Regression
  currSurv <- survival::Surv(time = currData_dt$PFS_time, event = currData_dt$`progression (1 or 0)`)
  
  ### Model Fit
  currFit <- survival::coxph(formula = currSurv ~ value, data = currData_dt)
  
  ### Extract p value
  currP_v <- summary(currFit)$coefficients[5]
  
  ### Fit test
  currTest <- survival::cox.zph(currFit)
  currCox_gg <- suppressWarnings(survminer::ggcoxzph(currTest)$`1`)
  currCox_gg@labels$title <- paste0(currCox_gg@labels$title, "\nbeta p: ", round(currP_v, digits = 5))
  
  ### Fit test 2
  currMartin_gg <- survminer::ggcoxdiagnostics(fit = currFit, type = "martingale", linear.predictions = F, ox.scale = "linear.predictions")
  
  ### Combine
  currCombo_gg <- ggpubr::ggarrange(curr_gg, currCox_gg, currMartin_gg, nrow = 1)
  currCombo_gg <- ggpubr::annotate_figure(currCombo_gg, top = text_grob(label = paste0("PFS Time vs. ", currPop_v, "\n", currMeasure_v), size = 24))
  
  pdf(file = file.path(baseDir_v, "results/fig2/", paste0(currFigName_v, "_scatterChecks_", currMeasure_v, "_", currPop_v, ".pdf")), width = 18)
  print(currCombo_gg)
  dev.off()
  
  pdf(file = file.path(baseDir_v, "results/fig2/", paste0(currFigName_v, "_scatter_", currMeasure_v, "_", currPop_v, ".pdf")))
  print(curr_gg)
  dev.off()
  
} # for i

###
### KM Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

timeCol_v <- "PFS_time"
eventCol_v <- "progression (1 or 0)"

### Specify columns
metaCols_v <- colnames(survGrp_lsdt[[1]])[1:(grep("CD8", colnames(survGrp_lsdt[[1]]))[1]-1)]
measureCols_v <- setdiff(colnames(survGrp_lsdt[[1]]), metaCols_v)

for (i in 1:nrow(kmPops_dt)) {
  
  ### Get info
  currPop_v <- kmPops_dt$pop[i]
  currMeasure_v <- kmPops_dt$measure[i]
  currFigName_v <- kmPops_dt$fig[i]
  currGrp_dt <- survGrp_lsdt[[currMeasure_v]]
  
  ### Subset Data
  currData_dt <- melt_lsdt[[currMeasure_v]][variable == currPop_v,]
  
  ### Subset NAs
  currData_dt <- currData_dt[!is.na(value),]
  
  ### Get Median
  currMed_v <- survMed_dt[Column == currPop_v, get(currMeasure_v)]
  
  ### Make Group
  currGrp_v <- paste0(currPop_v, "_grp")
  
  ### Merge
  currMerge_dt <- merge(currData_dt, currGrp_dt[,mget(c("PatientID", "Cycle", "CB", currGrp_v))],
                       by = c("PatientID", "Cycle", "CB"))
  
  ### Factorize
  currMerge_dt[[currGrp_v]] <- factor(currMerge_dt[[currGrp_v]], levels = c("lo", "hi"))
  
  ## Make survival object
  currSurv <- survival::Surv(time = currMerge_dt[[timeCol_v]], event = currMerge_dt[[eventCol_v]])
  
  ## Make formulas
  currFormula <- as.formula(paste0("currSurv ~ ", currGrp_v))
  
  ## Fit model
  currFit <- survminer::surv_fit(currFormula, data = currMerge_dt)
  
  ## Run log-rank
  currLR <- survival::survdiff(currFormula, data = currMerge_dt)
  currP_v <- currLR$pvalue
  
  ## Make plot
  curr_gg <- suppressMessages(suppressWarnings(
    survminer::ggsurvplot(fit = currFit,
                          data = currMerge_dt,
                          pval = T,
                          title = paste0(currMeasure_v, " - ", currPop_v, 
                                         "\nMedian: ", round(currMed_v, digits = 3)),
                          surv.median.line = "hv",
                          ggtheme = survTheme(),
                          risk.table = T,
                          xlab = timeCol_v,
                          legend.title = currPop_v,
                          legend.labs = c("below Med", "above Med"),
                          legend = c(0.8,0.8),
                          palette = c("blue", "red"))))
  
  pdf(file = file.path(baseDir_v, "results/fig2/", paste0(currFigName_v, "_KM_", currMeasure_v, "_", currPop_v, ".pdf")), 
      width = 7, height = 7, onefile = F)
  suppressMessages(suppressWarnings(print(curr_gg)))
  dev.off()
  
} # for i
