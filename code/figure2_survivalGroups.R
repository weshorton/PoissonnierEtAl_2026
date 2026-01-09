#####################################################################
########################## SURVIVAL GROUPS ##########################
#####################################################################

### Read in data and calculate survival groups based on median expression

source("./code/figure2_sourceRef.R")
library(data.table)

###
### Set Up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Read in data
data_dt <- as.data.table(readxl::read_xlsx("./data/fig2/data.xlsx", sheet = "all"))

### Read in reference - indicates columns to use and their data type
ref_dt <- fread("./data/fig2/flowDataTypes.txt")

### Extract metadata
meta_dt <- unique(data_dt[,mget(metaCols_v)])

### Set survival columns
timeCol_v <- "PFS_time"
eventCol_v <- "progression (1 or 0)"

### Remove samples
### These patients are not able to be used for survival comparisons
meta_dt <- meta_dt[BR != "NE" | is.na(BR),]

### Change all percentages from decimal to percent
for (i in 1:ref_dt[,.N]) {
  currPop_v <- ref_dt[i,Column]
  if (ref_dt[i,Measure] == 'pct') {
    set(data_dt, j = currPop_v, value = data_dt[[currPop_v]]*100)
  } # fi
} # for i

### Calculate Ratios Between Timepoints
data_lsdt <- komenFlowWrangle(data_dt = data_dt,
                              ref_dt = ref_dt,
                              grep_v = "CD8",
                              abs_v = "Screen",
                              num_v = c("C1D1", "C2D1", "C2D1"),
                              denom_v = c("Screen", "Screen", "C1D1"))

### Get metadata and measurement columns
metaCols_v <- colnames(data_lsdt[[1]])[1:(grep("CD8", colnames(data_lsdt[[1]]))[1]-1)]
measureCols_v <- setdiff(colnames(data_lsdt[[1]]), metaCols_v)

###
### Calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Empty lists to hold output
groupData_lsdt <- list()
groupMedians_lsv <- list()

### Run calculations for each comparison
for (i in 1:length(data_lsdt)) {
  
  ### Get single data occurrence
  currMeasureName_v <- names(data_lsdt)[i]
  currData_dt <- data_lsdt[[currMeasureName_v]]
  
  ### Empty var to hold
  currMedians_v <- NULL
  
  ### Calculate medians and split into groups for each column
  for (j in 1:nrow(ref_dt)) {
    
    ### Get column
    currColumn_v <- ref_dt$Column[j]
    
    ### Make the group name
    currGroup_v <- paste0(currColumn_v, "_grp")
    
    ### Check
    if (!(currColumn_v %in% colnames(currData_dt))) {
      stop(sprintf("Missing column: %s\n", currColumn_v))
    }
    
    ### Calculate Median
    currMed_v <- median(currData_dt[[currColumn_v]], na.rm = T)
    
    ### Group patients using median
    currData_dt[,(currGroup_v) := ifelse(get(currColumn_v) > currMed_v, "hi", "lo")]
    
    ## Add to output
    currMedians_v <- c(currMedians_v, currMed_v)
    
  } # for j
  
  ### Extract group data for output
  groupData_lsdt[[currMeasureName_v]] <- currData_dt[,mget(c(metaCols_v, grep("_grp", colnames(currData_dt), value = T)))]
  
  ### Extract medians for output
  names(currMedians_v) <- ref_dt$Column
  groupMedians_lsv[[currMeasureName_v]] <- currMedians_v
  
} # for i

###
### Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

###
### Wrangle Medians List ~~~~
###

### Change medians vectors into a tables
groupMedians_lsdt <- sapply(groupMedians_lsv, function(x) {
  data.table("Column" = names(x), "Median" = x)
}, simplify = F, USE.NAMES = T)

### Merge
groupMedians_dt <- merge(groupMedians_lsdt$Screen, groupMedians_lsdt$S.C1D1, by = "Column",
                         suffixes = c("_Screen", "_S.C1D1"))
colnames(groupMedians_dt) <- gsub("Median_", "", colnames(groupMedians_dt))

groupMedians_dt <- merge(groupMedians_dt, groupMedians_lsdt$S.C2D1, by = "Column")
colnames(groupMedians_dt)[colnames(groupMedians_dt) == "Median"] <- "S.C2D1"

groupMedians_dt <- merge(groupMedians_dt, groupMedians_lsdt$C1D1.C2D1, by = "Column")
colnames(groupMedians_dt)[colnames(groupMedians_dt) == "Median"] <- "C1D1.C2D1"

###
### Write group data ~~~~~~~~
###

### Create workbook
wb <- openxlsx::createWorkbook()

### Load sheets
invisible(lapply(seq_along(groupData_lsdt), function(x) {
  openxlsx::addWorksheet(wb = wb, sheetName = names(groupData_lsdt)[x])
  openxlsx::writeData(wb = wb, sheet = x, x = groupData_lsdt[[x]])
}))

### Save workbook
openxlsx::saveWorkbook(wb = wb, file = "./data/fig2/survivalGroups.xlsx", overwrite = T)

###
### Write medians ~~~~~~~~~~~
###

writexl::write_xlsx(x = groupMedians_dt, path = "./data/fig2/survivalMedians.xlsx")
