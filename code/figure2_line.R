###############################################################
######################## FIGURE 2 LINE ########################
###############################################################

### Set to "/" for CodeOcean and "./" for github
baseDir_v <- "./"

### Make line plots for 
  ### main figure 2 (panels A, B, C, and I)
  ### supplemental figure 2 (panels P (left), P (right), and Q)

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

### Create metadata
meta_dt <- unique(data_dt[,mget(c("PatientID", "BR", "CB"))])

### Remove non-evaluable patients
toRm_v <- meta_dt[is.na(BR) | BR == "NE", PatientID]
data_dt <- data_dt[!(PatientID %in% toRm_v),]

### Read in stats results and split
# results_dt <- as.data.table(readxl::read_xlsx("./data/fig2/contrasts-for-135-outcomes_12JUN2025.xlsx", sheet = "Q1toQ5"))
# results_lsdt <- sapply(paste0("Q", 1:5), function(x) results_dt[Label == x,], simplify = F, USE.NAMES = T)

### Read in data type info
flow_dt <- fread(file.path(baseDir_v, "data/fig2/flowDataTypes.txt"))
luminex_dt <- fread(file.path(baseDir_v, "data/fig2/luminexDataTypes.txt"))
ref_dt <- rbind(flow_dt, luminex_dt)

### Plot populations
figurePops_dt <- data.table("pop" = c("Non_classicalMonocyte_CD14pCD16pp", "CSF1", "CCL2", "CD45p_PDL1",
                                      "CD4_PD1", "CD4_TemTreg_PD1", "CD8_PD1"),
                            "fig" = c("Fig2A", "Fig2B", "Fig2C", "Fig2I", 
                                      "suppFig2P_left", "suppFig2P_right", "suppFig2Q"))

###
### Wrangle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Metadata columns are the first column through the "merger" column
allMetaCols_v <- colnames(data_dt)[1:which(colnames(data_dt) == "merger")]
meltCols_v <- c("PatientID", "Cycle", "CB")
measureCols_v <- setdiff(colnames(data_dt), allMetaCols_v)

### Melt
melt_dt <- melt(data_dt[,mget(c(meltCols_v, measureCols_v))], id.vars = meltCols_v)

### Total unique patients
uniqPt_v <- unique(melt_dt$PatientID)
if (length(uniqPt_v) != 31) stop("Bad patient substitution")

### Update Screen to LI
melt_dt[Cycle == "Screen", Cycle := "LI"]

### Mark changes (classify into >10%+, >10%-, <10%+/-)
change_lsdt <- lapply(names(timeCompare_lsv), function(x) {
  markChange(data_dt = melt_dt, var_lsv = timeCompare_lsv[x], col_v = "Cycle", mCol_v = "PatientID", metaCol_v = "CB", change_v = 10)})
change_dt <- do.call(rbind, change_lsdt)
colnames(change_dt)[colnames(change_dt) == "variable"] <- "measureVar"

## Set levels for plot order
melt_dt$Cycle <- factor(melt_dt$Cycle, levels = c("LI", "C1D1", "C2D1", "C4D1"))

## Get comparisons to run
comparisons_v <- unique(change_dt$Compare)

###
### Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

for (i in 1:nrow(figurePops_dt)) {
  
  ### Get current run info
  currPop_v <- figurePops_dt[i, pop]
  currFig_v <- figurePops_dt[i, fig]
  currMeasureType_v <- ref_dt[Column == currPop_v, Measure]
  
  ### Subset data
  currData_dt <- change_dt[measureVar == currPop_v,]
  rmCols_v <- c("Change", "pctChange")
  currData_dt[,(rmCols_v) := NULL]
  
  ### Melt again
  currMelt_dt <- melt(currData_dt, measure.vars = c("value_1", "value_2"))
  currMelt_dt$Compare <- factor(currMelt_dt$Compare, levels = c("LI.C1D1", "C1D1.C2D1", "C2D1.C4D1"))
  
  ### Multiply percents
  if (currMeasureType_v == "pct") currMelt_dt$value <- currMelt_dt$value * 100
  
  ### Make table of patients in each change group for each timepoint comparison
  currChanges_lsdt <- list()
  cols_v <- c("Dir", "PD", "SD_PR")
  for (j in 1:length(comparisons_v)) {
    
    ### Subset for time comparison
    comp_v <- comparisons_v[j]
    currCompData_dt <- currData_dt[Compare == comp_v,]
    
    ### Make table
    rec_dt <- tableChange(change_dt = currCompData_dt, cCol_v = "Direction", metaCol_v = "CB", var_v = currPop_v, out_v = "dt")
    if (is.null(rec_dt)) next
    
    ### Add empties
    for (c_v in setdiff(cols_v, colnames(rec_dt))) rec_dt[[c_v]] <- NA
    
    ### Add comparison label
    rec_dt$Compare <- comp_v
    
    ### Sort columns
    currChanges_lsdt[[comp_v]] <- rec_dt[,mget(c("Compare", cols_v))]
    
  } # for j
  
  ### Combine and reshape
  all_dt <- do.call(rbind, currChanges_lsdt)
  wide_dt <- reshape(all_dt, timevar = "Compare", idvar = "Dir", direction = "wide")
  wide_dt[is.na(wide_dt)] <- 0
  
  ### Turn table into grob
  currChanges_grob <- tableGrob(wide_dt)
  
  ## Make plot
  curr_gg <- ggplot(data = currMelt_dt[!is.na(value)], 
                    aes(x = variable, y = value, group = PatientID, color = Direction)) +
    geom_line(linewidth = 0.3) +
    facet_wrap(~ Compare, nrow = 1) +
    scale_x_discrete(limits = c("value_1", "value_2"), expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    ggtitle(currPop_v) +
    guides(color = "none") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 8), 
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          strip.background = element_blank(),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14),
          strip.text.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0, size = 11)) +
    scale_color_manual(values = lineGroupColors_v, labels = names(lineGroupColors_v))
  
  if (currPop_v == "CSF1") {
    curr_gg <- curr_gg + expand_limits(y = 0)
  }
  
  ### Save plot
  ggsave(filename = file.path(baseDir_v, "results/fig2", paste0(currFig_v, ".pdf")), plot = curr_gg, 
         device = "pdf", width = 3, height = 3)
  
  ### Save table
  pdf(file = file.path(baseDir_v, "results/fig2", paste0("table_", currFig_v, ".pdf")), width = 10, height = 10)
  grid.draw(currChanges_grob)
  dev.off()
  
} # for i

