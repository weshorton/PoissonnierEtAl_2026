###################################################################
###################### FIGURE 2 SCATTER ###########################
###################################################################

### Make line plots for 
  ### supplemental figure 2 panels S, T, U, and V

source("./code/figure2_sourceRef.R")
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

###
### Set Up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Read in data
data_dt <- as.data.table(readxl::read_xlsx("./data/fig2/data.xlsx", sheet = "all"))

### Create metadata
meta_dt <- unique(data_dt[,mget(c("PatientID", "BR", "CB"))])

### Remove non-evaluable patients
toRm_v <- meta_dt[BR == "NE", PatientID]
data_dt <- data_dt[!(PatientID %in% toRm_v),]

### Check
uniqPt_v <- unique(data_dt$PatientID)
if (length(uniqPt_v) != 35) stop("Bad patient substitution")

### Read in data type info
flow_dt <- fread("./data/fig2/flowDataTypes.txt")
luminex_dt <- fread("./data/fig2/luminexDataTypes.txt")
ref_dt <- rbind(flow_dt, luminex_dt)

# ### Info for converting names to time
# time_lsv <- list("LIC1" = c("Screen", "C1D1"),
#                  "LIC2" = c("Screen", "C2D1"),
#                  "C1C2" = c("C1D1", "C2D1"))

### Change values to pct for those that require it
pct_dt <- ref_dt[Measure == "pct",]
for (i in 1:pct_dt[,.N]) {
  currPop_v <- pct_dt[i,Column]
  set(data_dt, j = currPop_v, value = data_dt[[currPop_v]]*100)
}

### Plot populations
figurePops_dt <- data.table("pop" = c("CD4_TemTh1_PD1", "MyeloidDCs", "Monocytes_PDL1", "ClassicalMonocytes_CD14ppCD16n_PDL1"),
                            "measure" = c("Screen", "S.C2D1", "C1D1.C2D1", "C1D1.C2D1"),
                            "fig" = c("FigS2S", "FigS2T", "FigS2U", "FigS2V"))

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

for (i in 1:nrow(figurePops_dt)) {
  
  ### Get info
  currPop_v <- figurePops_dt$pop[i]
  currMeasure_v <- figurePops_dt$measure[i]
  currFigName_v <- figurePops_dt$fig[i]
  
  ### Subset Data
  currData_dt <- melt_lsdt[[currMeasure_v]][variable == currPop_v,]
  
  ### Subset NAs
  currData_dt <- currData_dt[!is.na(value),]
  
  ### Column for shapes
  currData_dt[`progression (1 or 0)` == 0, Prog := "0"]
  currData_dt[`progression (1 or 0)` == 1, Prog := "1"]
  
  ### Make plot
  curr_gg <- ggplot(data = currData_dt, aes(x = value, y = PFS_time, shape = Prog, color = Prog)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("grey", "black"), breaks = c("0", "1")) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 20), 
                            plot.subtitle = element_text(hjust = 0.5, size = 14), 
                            axis.text = element_text(size = 16), 
                            axis.title = element_text(size = 18), 
                            legend.text = element_text(size = 16), 
                            legend.title = element_text(size = 18)) +
    labs(x = "Expression", y = "PFS time (days)", shape = "Prog") +
    ggtitle("PFS Time Scatter")
    #ggtitle(paste0("PFS Time vs. ", currPop_v, "\n", currMeasure_v, " Expression (", nrow(currData_dt), " Pts)"))
  
  #print(curr_gg)
    
  # pdf(file = file.path("./tmpSurv/", paste0(currMeasure_v, "_", currPop_v, ".pdf")))
  # print(curr_gg)
  # dev.off()
  
  ### Cox Regression
  currSurv <- survival::Surv(time = currData_dt$PFS_time, event = currData_dt$`progression (1 or 0)`)
  
  ### Model Fit
  currFit <- survival::coxph(formula = currSurv ~ value, data = currData_dt)
  
  ### Extract p value
  currP_v <- summary(currFit)$coefficients[5]
  print(summary(currFit))
  
  ### Fit test
  currTest <- survival::cox.zph(currFit)
  currCox_gg <- survminer::ggcoxzph(currTest)$`1`
  currCox_gg@labels$title <- paste0(currCox_gg@labels$title, "\nbeta p: ", round(currP_v, digits = 5))
  
  ### Fit test 2
  currMartin_gg <- survminer::ggcoxdiagnostics(fit = currFit, type = "martingale", linear.predictions = F, ox.scale = "linear.predictions")
  
  ### Combine
  currCombo_gg <- ggpubr::ggarrange(curr_gg, currCox_gg, currMartin_gg, nrow = 1)
  currCombo_gg <- ggpubr::annotate_figure(currCombo_gg, top = text_grob(label = paste0("PFS Time vs. ", currPop_v, "\n", currMeasure_v), size = 24))
  
  pdf(file = file.path("./fig2/surv/", paste0(currMeasure_v, "_", currPop_v, ".pdf")), width = 18)
  print(currCombo_gg)
  dev.off()
  

} # for i
