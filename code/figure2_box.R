###############################################################
######################## FIGURE 2 BOX #########################
###############################################################

### Make line plots for 
  ### main figure 2 (panels D, E, and J)
  ### supplemental figure 2 (panels N, O, and R)

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

### Change NE to NA
data_dt[BR == "NE", BR := NA]

# ### Read in stats results and split
# rawFile_v <- "./data/fig2/contrasts-for-135-outcomes_12JUN2025.xlsx"
# rawSheets_v <- grep("^Q[1-5]$", readxl::excel_sheets(rawFile_v), value = T)
# rawResults_lsdt <- lapply(rawSheets_v, function(x) {
#   setDT(readxl::read_excel(rawFile_v, sheet = x))
# })
# names(rawResults_lsdt) <- rawSheets_v
# 
# log2File_v <- "./data/fig2/contrasts-for-135-outcomes_log2_12JUN2025.xlsx"
# log2Sheets_v <- grep("^Q[1-5]_log2$", readxl::excel_sheets(log2File_v), value = T)
# log2Results_lsdt <- lapply(log2Sheets_v, function(x) {
#   setDT(readxl::read_excel(log2File_v, sheet = x))
# })
# names(log2Results_lsdt) <- log2Sheets_v
# 
# results_lslsdt <- list("raw" = rawResults_lsdt, "log2" = log2Results_lsdt)

### Read in data type info
flow_dt <- fread("./data/fig2/flowDataTypes.txt")
luminex_dt <- fread("./data/fig2/luminexDataTypes.txt")
ref_dt <- rbind(flow_dt, luminex_dt)

### Info for converting names to time
time_lsv <- list("LIC1" = c("Screen", "C1D1"),
                 "LIC2" = c("Screen", "C2D1"),
                 "C1C2" = c("C1D1", "C2D1"))

### Where to put N labels ("top", "bottom", or "none)
labelPosition_v <- "none"

###
### Wrangle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Split data (a little unnecessary...)
data_lsdt <- komenFlowWrangle(data_dt = data_dt,
                              ref_dt = ref_dt,
                              grep_v = "CD8",
                              num_v = NULL,
                              abs_v = c("Screen", "C1D1", "C2D1"))

### Specify columns
metaCols_v <- colnames(data_lsdt[[1]])[1:(grep("CD8", colnames(data_lsdt[[1]]))[1]-1)]
measureCols_v <- setdiff(colnames(data_lsdt[[1]]), metaCols_v)
analytes_v <- luminex_dt$Column
pops_v <- setdiff(measureCols_v, analytes_v)
pctPops_v <- flow_dt[Measure == "pct",Column]
numPops_v <- flow_dt[Measure == "num",Column]

### Combine all the data (should just be the same as input, without C4D1)
data_dt <- do.call(rbind, data_lsdt)

### Remove NA patients
data_dt <- data_dt[!is.na(BR),]
data_dt <- data_dt[BR != "NA",]

### Check
uniqPt_v <- unique(data_dt$PatientID)
if (length(uniqPt_v) != 31) stop("Bad patient substitution")

### Change values to pct for those that require it
pct_dt <- ref_dt[Measure == "pct",]
for (i in 1:pct_dt[,.N]) {
  currPop_v <- pct_dt[i,Column]
  set(data_dt, j = currPop_v, value = data_dt[[currPop_v]]*100)
}

###
### Plot Populations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

toPlot_lslsv <- list("fig2D" = list("Question" = "Q3", 
                                    "Class" = "plasma",
                                    "Scale" = "raw",
                                    "Time" = "C1C2",
                                    "Effect" = "Time",
                                    "Pops" = list(one = c("IFNG", "TNFSF13B"),
                                                  two = c("CXCL10", "CXCL11", "TNFRSF1B"),
                                                  three = "IL15",
                                                  four = c("CCL22", "CCL24"),
                                                  five = c("TNFRSF8"),
                                                  six = "IL2RA"),
                                    "Size" = c(6, 10), "LabSize" = c(6, 10)),
                     
                     
                     "fig2E" = list("Question" = "Q3", 
                                    "Class" = "bloodPct",
                                    "Scale" = "raw",
                                    "Time" = "C1C2",
                                    "Effect" = "Time",
                                    "Pops" = c("CD14nCD16n", "NKCells"),
                                    "Size" = c(6, 4), "LabSize" = c(10, 6)),
                     
                     
                     "fig2J" = list("Question" = "Q4", 
                                    "Class" = "bloodNum",
                                    "Scale" = "raw",
                                    "Time" = "C1C2",
                                    "Effect" = "Grp",
                                    "Pops" = c("CD4_TemTreg_PD1"),
                                    "Size" = c(6, 4), "LabSize" = c(10, 4)),
                     
                     
                     "figS2N" = list("Question" = "Q5", 
                                     "Class" = "bloodNum",
                                     "Scale" = "log2",
                                     "Time" = "C1C2",
                                     "Effect" = "Time",
                                     "Pops" = c("ClassicalMonocytes_CD14ppCD16n_PDL1", 
                                                "Non_classicalMonocyte_CD14pCD16pp_PDL1",
                                                "IntermediateMonocyte_CD14pCD16p_PDL1"),
                                     "Size" = c(6, 12), "LabSize" = c(8, 12)),
                     
                     
                     "figS2O" = list("Question" = "Q5",
                                     "Class" = "bloodNum",
                                     "Scale" = "log2",
                                     "Time" = "C1C2",
                                     "Effect" = "Time",
                                     "Pops" =  c("PlasmacytoidDCs_PDL1", "CD141p_mDCs_PDL1", "CD1cp_mDCs_PDL1"),
                                     "Size" = c(6, 12), "LabSize" = c(8, 12)),
                     
                     
                     "figS2R" = list("Question" = "Q4", 
                                     "Class" = "bloodNum",
                                     "Scale" = "raw",
                                     "Time" = "C1C2",
                                     "Effect" = "Grp",
                                     "Pops" = list(one = "CXCL9", two = "TNFSF13B", three = "TSLP"),
                                     "Size" = c(6, 6), "LabSize" = c(8, 6)))


###
### Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

counts_lsdt <- list()
for (i in 1:length(toPlot_lslsv)) {
  
  ### Get info
  currFigName_v <- names(toPlot_lslsv)[i]
  currInfo_lsv <- toPlot_lslsv[[currFigName_v]]
  
  currQ_v <- currInfo_lsv$Question
  currClass_v <- currInfo_lsv$Class
  currScale_v <- currInfo_lsv$Scale
  currTimeName_v <- currInfo_lsv$Time
  currEffect_v <- currInfo_lsv$Effect
  currPops_v <- currInfo_lsv$Pops
  currSize_v <- currInfo_lsv$Size
  currLabSize_v <- currInfo_lsv$LabSize
  
  # ### Get results
  # currRes_dt <- results_lslsdt[[currScale_v]][[grep(currQ_v, names(results_lslsdt[[currScale_v]]), value = T)]]
  # 
  # ### Subset
  # currRes_dt <- currRes_dt[Outcome %in% unlist(currPops_v),]
  # currRes_dt$Outcome <- factor(currRes_dt$Outcome, levels = unlist(currPops_v))
  
  ### Get effect and color
  if (currEffect_v == "Grp") {
    fillCol_v <- "CB"
  } else if (currEffect_v %in% c("Time", "Int")) {
    fillCol_v <- "Cycle"
  } # fi
  
  ### Get Times
  currTime_v <- time_lsv[[currTimeName_v]]
  
  ### Get y-scale
  if (currScale_v == "log2") {
    yVal_v <- "exp"
  } else {
    yVal_v <- "pct"
  } # fi
  
  ### Get data to plot, melt, subset for times
  currData_dt <- data_dt[,mget(c(metaCols_v, unlist(currPops_v)))]
  currMelt_dt <- melt(currData_dt, measure.vars = unlist(currPops_v))
  currMelt_dt <- currMelt_dt[Cycle %in% currTime_v,]
  
  ### Add facets
  if (class(currPops_v) == "list") {
    
    currMelt_dt$facet <- character()
    for (f_v in names(currPops_v)) {
      currMelt_dt[variable %in% currPops_v[[f_v]], facet := f_v]
    } # for f_v
    currMelt_dt$facet <- factor(currMelt_dt$facet, levels = names(currPops_v))
    
  } # fi list
  
  ### Change Screen to LI
  currMelt_dt[Cycle == "Screen", Cycle := "LI"]
  currNewTime_v <- gsub("Screen", "LI", currTime_v)
  names(timeColors_v) <- gsub("Screen", "LI", names(timeColors_v))

  ### Make sure factors are correct
  currMelt_dt$Cycle <- factor(currMelt_dt$Cycle, levels = currNewTime_v)
  currMelt_dt$CB <- factor(currMelt_dt$CB, levels = c("PD", "SD_PR"))
  currMelt_dt$variable <- factor(currMelt_dt$variable, levels = unlist(currPops_v))
  
  ### Subset colors to only be in time
  currTimeColors_v <- timeColors_v[names(timeColors_v) %in% currNewTime_v]
  
  ### Set X-variable
  if (currQ_v %in% c("Q4", "Q5")) {
    xVar_v <- "CB"
  } else {
    xVar_v <- "variable"
  } # fi
  
  ### Remove NA values
  currMelt_dt <- currMelt_dt[!is.na(value)]
  
  ### Get counts
  counts_dt <- as.data.table(table(currMelt_dt[,mget(unique(c(fillCol_v, xVar_v)))]))
  counts_lsdt[[currFigName_v]] <- counts_dt
  
  ### Add them
  currMelt_dt <- merge(currMelt_dt, counts_dt, by = unique(c(fillCol_v, xVar_v)), sort = F)
  
  vars_v <- unique(c(xVar_v, fillCol_v))
  if (labelPosition_v == "top") {
    ### Get 3rd quartile y value
    currMelt_dt[,yPos := quantile(value)[4], by = vars_v]
  } else if (labelPosition_v == "bottom") {
    ### Get minimum y value
    currMelt_dt[,yPos := min(value), by = vars_v]
  } else if (labelPosition_v != "none") {
    stop("Bad value for labelPosition_v")
  }
  
  ### Base Plot
  curr_gg <- ggplot(data = currMelt_dt, aes(x = !!sym(xVar_v), y = value, fill = !!sym(fillCol_v))) +
    ggtitle(paste0("Significant Effect of ", currEffect_v, "\nUsing ", gsub("_.*", "", currQ_v))) +
    geom_boxplot(fatten = 1, lwd = 0.3) +
    labs(x = NULL) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 18), 
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          axis.text.y = element_text(angle = 0, size = 14),
          axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 14),
          legend.key.size = unit(1.5, "lines"))
  
  ### Add count label
  if (labelPosition_v == "top") {
    curr_gg <- curr_gg +
      geom_text(aes(y = yPos, label = N), position = position_dodge(width = 0.75), vjust = -0.5, size = 3)
  } else if (labelPosition_v == "bottom") {
    curr_gg <- curr_gg +
      geom_text(aes(y = yPos, label = N), position = position_dodge(width = 0.75), vjust = 1.5, size = 3)
  }
  
  ### Y-axis scale
  if (currScale_v == "log2") {
    curr_gg <- curr_gg + scale_y_continuous(trans = "log2", labels = label_log2)
  }
  
  ### Add Color
  if (currEffect_v == "Grp") {
    curr_gg <- curr_gg + scale_fill_manual(values = grpColors_v, labels = names(grpColors_v))
  } else if (currEffect_v %in% c("Time", "Int")) {
    curr_gg <- curr_gg + scale_fill_manual(values = currTimeColors_v, labels = names(currTimeColors_v))
  } 
  
  ### Facet
  if (currFigName_v == "fig2D") {
    curr_gg <- curr_gg + facet_wrap("~facet", nrow = 1, scales = "free")
  }
  
  if (currFigName_v %in% c("figS2N", "figS2O", "figS2R")) {
    curr_gg <- curr_gg + facet_wrap("~variable", nrow = 1, scales = "free")
  }
  
  ### Make version with no label and no legend
  currNoLab_gg <- curr_gg + theme(axis.text.x = element_blank())
  currNoLab_gg <- currNoLab_gg + guides(fill = "none")
  
  ### Handle label position name
  if (labelPosition_v == "top") {
    currFigName_v <- paste0(currFigName_v, "_countsTop")
  } else if (labelPosition_v == "bottom") {
    currFigName_v <- paste0(currFigName_v, "_countsBottom")
  } else {
    currFigName_v <- paste0(currFigName_v, "_noCounts")
  }
  
  ### Save
  ggsave(filename = file.path("./fig2", paste0(currFigName_v, "_noLabel.pdf")), device = "pdf",
         plot = currNoLab_gg, height = currSize_v[1], width = currSize_v[2])
  
  ggsave(filename = file.path("./fig2", paste0(currFigName_v, ".pdf")), device = "pdf",
         plot = curr_gg, height = currLabSize_v[1], width = currLabSize_v[2])
  
} # for i

### Write
fig2box_wb <- openxlsx::createWorkbook()
invisible(lapply(seq_along(counts_lsdt), function(x) {
  openxlsx::addWorksheet(wb = fig2box_wb, sheetName = names(counts_lsdt)[x])
  openxlsx::writeData(wb = fig2box_wb, sheet = x, x = counts_lsdt[[x]])
}))
openxlsx::saveWorkbook(wb = fig2box_wb, file = "./fig2/boxplotSampleCounts.xlsx", overwrite = T)
