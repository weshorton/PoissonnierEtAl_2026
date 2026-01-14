####################################################################
######################## FIG2 SOURCE REF ###########################
####################################################################

### Miscellaneous colors, column groups, functions, etc. used in the plotting/wrangling code for Figure 2
library(ggplot2)

###
### Metadata columns ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Columns from the data sheets that contain metadata information
metaCols_v <- c("PatientID",                  # deidentified patient ID
                "CB",                         # clinical benefit (PD and SD combined)
                "BR",                         # RECIST best response
                "Best Response_Park",         # Extra best response to for checks
                "progression (1 or 0)",       # Censor
                "PFS_time",                   # Progression-free survival time
                "TOT",                        # Time on trial
                "Dstart",                     # Date of trial start
                "EOT Date",                   # Date of trial end               
                "Progression",                # Progression record
                "Date_progression",           # Progression Date
                "OFF Study",                  # Reason for study removal
                "EOT Reason")                 # Reason for end of trial

### Time series comparisons
timeCompare_lsv <- list("LI.C1D1" = c("LI", "C1D1"),
                        "C1D1.C2D1" = c("C1D1", "C2D1"),
                        "C2D1.C4D1" = c("C2D1", "C4D1"))

### Colors to use for line plots
lineGroupColors_v <- c("#47AD48", "#D7222B", "#0A0B09")
names(lineGroupColors_v) <- c("UP", "DN", "NC")

### Colors for time designations (boxplot)
timeColors_v <- c("Screen" = "#33B44A", "C1D1" = "#EF3324", "C2D1" = "#006164", "C4D1" = "#9654A2")

### Colors for response group designations (boxplot)
grpColors_v <- RColorBrewer::brewer.pal(8, "Dark2")[c(2,4)]
names(grpColors_v) <- c("PD", "SD_PR")

survTheme <- function() {
  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 18), 
          plot.subtitle = element_text(hjust = 0.5, size = 14), 
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14))
}

###
### Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Log2 output
label_log2 <- function(x) parse(text = paste0('2^', log(x, 2)))

###################################################################################################

#' Komen Flow Wrangle
#' @description Remove extra columns and perform a few calculations
#' @param data_dt input data with all the patients and measures
#' @param ref_dt data.table with columns to keep. This is pre-determined based on what measurements you're interested.
#' @param grep_v pattern to use to search for meta columns. Meta columns are all columns BEFORE this grep is found
#' @param abs_v vector of timepoints to extract absolute data for
#' @param num_v vector of timepoints that are numerator log2 ratio data for
#' @param denom_v vector of timepoints that are the denominator for log2 ratio calculations
#' @param ratio_v two options. "log2" - log2 ratio; "pctChange" - percent change (num_v - denom_v) / denom_v * 100
#' @details
#' Input data has all data used in project - clinical information, flow output, and luminex output
#' We are only running survival analysis on flow output, so need to remove other stuff
#' We need to calculate measurements relative to something else to compare changes
#' @export
komenFlowWrangle <- function(data_dt, 
                             ref_dt, 
                             grep_v = "CD8", 
                             abs_v = NA,
                             num_v = c("C1D1", "C2D1", "C2D1"), 
                             denom_v = c("Screen", "Screen", "C1D1"), 
                             ratio_v = "log2") {
  
  ##
  ## SUBSET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  ## Identify meta columns
  ## The first measurement call is CD8_Tem_CD8, so grep on CD8 to split
  firstGrep_v <- grep(grep_v, colnames(data_dt))[1]
  metaCols_v <- colnames(data_dt)[1:(firstGrep_v - 1)]
  
  ## Identify measure columns
  ## Depending on the ref table you loaded, this could be
  ## every measurement in the table, or some subset.
  measureCols_v <- ref_dt$Column
  
  ### Clean table
  data_dt <- data_dt[,mget(c(metaCols_v, measureCols_v))]
  
  ## Make sure measure columns are numeric
  ## The flow data has NA's recorded as n/a, so need to fix those so we can make everything numeric
  data_dt <- data_dt[, lapply(.SD, function(x) gsub("n/a", NA, x))]
  for (col_v in measureCols_v) set(data_dt, j = col_v, value = as.numeric(data_dt[[col_v]]))
  
  ##
  ## ABSOLUTE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  abs_lsdt <- list()
  if (!is.na(abs_v[1])) {
    for (i in 1:length(abs_v)) {abs_lsdt[[ abs_v[i] ]] <- data_dt[Cycle == abs_v[i],]}
  } # fi
  
  ##
  ## RATIO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  if (!is.null(num_v)) {
    
    ## Check equivalency
    if (length(num_v) != length(denom_v)) stop("Numerator and Denominator lengths do not match.")
    
    ratio_lsdt <- list()
    for (i in 1:length(num_v)) {
      base_dt <- data_dt[Cycle == denom_v[i],]
      time_dt <- data_dt[Cycle == num_v[i],]
      ratio_dt <- calcRatio(base_dt = base_dt, time_dt = time_dt, mCol_v = measureCols_v, ratio_v = ratio_v)
      ratio_lsdt[[ i ]] <- ratio_dt
    } # for i
    names(ratio_lsdt) <- gsub("creen", "", paste(denom_v, num_v, sep = "."))
    
    out_lsdt <- c(abs_lsdt, ratio_lsdt)
    
  } else {
    
    out_lsdt <- abs_lsdt
    
  } # fi
  
  ##
  ## OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  return(out_lsdt)
  
} # komenFlowWrangle

###################################################################################################

#' Calculate Ratio
#' @description Calculate ratio (log2 or pct change) between timepoint and baseline
#' @param base_dt baseline data
#' @param time_dt timepoint data
#' @param mCol_v measure columns that need ratios calculated
#' @param ratio_v two options. "log2" - log2 ratio; "pctChange" - percent change (num_v - denom_v) / denom_v
#' @details IMPORTANT!!!! Values listed as 0 will be turned to 0.0001 so that the ratio can be taken.
#' @return input table along with extra columns for each ratio
#' @export
calcRatio <- function(base_dt, time_dt, mCol_v, ratio_v) {
  
  ## Find patients in each
  inBaseOutTime_v <- setdiff(base_dt$PatientID, time_dt$PatientID)
  inTimeOutBase_v <- setdiff(time_dt$PatientID, base_dt$PatientID)
  
  ## Remove them
  if (length(inBaseOutTime_v) > 0) base_dt <- base_dt[!(PatientID %in% inBaseOutTime_v),]
  if (length(inTimeOutBase_v) > 0) time_dt <- time_dt[!(PatientID %in% inTimeOutBase_v),]
  
  ## Turn 0 to 0.0001
  base_dt[, (mCol_v) := lapply(.SD, function(x) ifelse(x <= 0, 0.0001, x)), .SDcols = mCol_v]
  time_dt[, (mCol_v) := lapply(.SD, function(x) ifelse(x <= 0, 0.0001, x)), .SDcols = mCol_v]
  
  ## Make sure they're the same order
  base_dt <- base_dt[order(PatientID)]
  time_dt <- time_dt[order(PatientID)]
  
  ## Calculate Ratio
  if (ratio_v == "log2") {
    ratio_df <- as.data.frame(log2(as.matrix(time_dt[,mget(mCol_v)]) / as.matrix(base_dt[,mget(mCol_v)])))
  } else if (ratio_v == "pctChange") {
    ratio_df <- as.data.frame(((as.matrix(time_dt[,mget(mCol_v)]) - as.matrix(base_dt[,mget(mCol_v)])) / as.matrix(base_dt[,mget(mCol_v)]))*100)
  } else {
    stop("Bad value for ratio_v")
  }
  
  ## Add back meta
  out_dt <- cbind(base_dt[,mget(setdiff(colnames(base_dt), mCol_v))], ratio_df)
  
  ## Return
  return(out_dt)
  
} # calcRatio

###################################################################################################

#' Mark Changes
#' @description For sets of time-points, determine if data is increasing or decreasing between them
#' @param data_dt data in melted format
#' @param col_v Which column is used to distinguish time points
#' @param var_lsv which values of col_v to use (should be list of length 1, with 2 elements). Must be named list.
#' @param mCol_v which column should be used to merge, along with 'variable'
#' @param metaCol_v meta.data columns that should be kept through transformation
#' @param change_v percent increase/decrease threshold. (e.g. 10 means 10% increase or decrease required)
#' @details
#' Label patients into down ("DN"), up ("UP"), and no change ("NC") based on percent change
#' @return table of same dimensions as data_dt, but with one extra column indicating the direction of change
#' @export
markChange <- function(data_dt, col_v = "Time", var_lsv, mCol_v = "PatientID", metaCol_v = "CB", change_v = 10) {
  
  ## Merge data
  m_dt <- merge(data_dt[get(col_v) == var_lsv[[1]][1],],
                data_dt[get(col_v) == var_lsv[[1]][2],],
                by = c(mCol_v, 'variable', metaCol_v), 
                suffixes = c("_1", "_2"),
                sort = F)
  
  ## Calculate Change and pct change
  m_dt$Change <- m_dt$value_2 - m_dt$value_1
  m_dt$pctChange <- m_dt$Change / m_dt$value_1 * 100
  
  ## New timepoint column
  m_dt$Compare <- names(var_lsv)[1]
  
  ## Classify
  m_dt[, ("Direction") := ifelse(pctChange >= change_v, 'UP',
                                 ifelse(pctChange <= -change_v, "DN", "NC"))]
  
  ## Return
  return(m_dt)
  
} # markChange

###################################################################################################

#' Table of Changes
#' @description Count numbers of changes
#' @param change_dt change table
#' @param cCol_v column of change
#' @param metaCol_v meta columns to keep
#' @param var_v value of 'variable' column
#' @param out_v 'character'. 'v' : output as vector with "\n"; "dt" : output as data.table
#' @details
#' Get counts of how many patients are in each change group
#' at each time point comparison
#' @return table to include in plot
#' @export
tableChange <- function(change_dt, cCol_v, metaCol_v, varCol_v = "measureVar", var_v, out_v = "v") {
  
  ## Convert NA
  change_dt[is.na(get(metaCol_v)), (metaCol_v) := "NA"]
  
  ## Make table
  table_tab <- table(change_dt[get(varCol_v) == var_v, mget(c(cCol_v, metaCol_v))])
  
  ## Reformat
  if (nrow(table_tab) > 0) {
    ## Turn to data.table
    table_df <- NULL
    for (c_v in 1:ncol(table_tab)) table_df <- cbind(table_df, table_tab[,c_v])
    rownames(table_df) <- rownames(table_tab)
    colnames(table_df) <- colnames(table_tab)
    table_dt <- convertDFT(table_df, newName_v = " ")
    if (out_v == "dt") {
      colnames(table_dt)[1] <- "Dir"
      return(table_dt)
    } else if (out_v == "v") {
      ## Turn to vectors
      table_v <- apply(table_dt, 1, function(x) paste(x, collapse = " - "))
      title_v <- paste(colnames(table_dt)[-1], collapse = " - ")
      ## Add new lines
      table_v <- paste(table_v, collapse = "\n")
      table_v <- paste(title_v, table_v, sep = "\n")
      ## Output
      return(table_v)
    } else {
      stop("Please select either 'v' or 'dt' for the out_v argument.")
    }
  } # fi
}

###################################################################################################

#' Convert between data.table and data.frame
#' @description 
#' Change data.tables into data.frames with specified row.names or
#' data.frames/matrices into data.tables, copying over row.names as the first column.
#' If data.frame/matrix doesn't have row.names, then no columns will be added to data.table
#' @param data_dft data in either data.table or data.frame format (can also be a matrix)
#' @param col_v character or numeric vector. if converting from dt to df, column name or index of which column to use as row.names.
#' NA (default) will use 1st column; NULL will not add rownames. NEW can also provide multiple columns here. The names and values
#' will be pasted together with "_" to create one column to become the rownames.
#' @param newName_v character vector. if converting from df/mat to dt, what to name new column. (default is "V1")
#' if newName_v is already a column name, will paste "_2" to end of newName_v.
#' @param rmCol_v boolean value indicating whether to remove the column used to make the rownames from the output table (T) or to leave it (F)
#' @param split_v character vector. default = NULL. Used for converting TO data.table. If there are multiple columns' worth of data in the rownames, split on this and assign results to multiple columns.
#' both the column names and the values by the provided delimiter.
#' @description
#' TODO! Add description
#' 
#' @return either a data.table or data.frame (opposite class of input)
#' @examples 
#' # Data
#' my_df <- data.frame("A" = 1:10, "B" = LETTERS[1:10], "C" = letters[11:20])
#' my_df2 <- my_df; rownames(my_df2) <- paste0("Row", 1:10)
#' my_mat <- as.matrix(my_df); my_mat2 <- as.matrix(my_df2)
#' my_dt <- data.table("AA" = 10:1, "BB" = LETTERS[5:14], "CC" = letters[20:11])
#' convertDFT(data_dft = my_df)
#' convertDFT(data_dft = my_df2)
#' convertDFT(data_dft = my_df2, newName_v = "Test")
#' convertDFT(data_dft = my_mat2, newName_v = "MatTest")
#' convertDFT(data_dft = my_dt)
#' convertDFT(data_dft = my_dt, col_v = NULL)
#' convertDFT(data_dft = my_dt, col_v = "BB")
#' convertDFT(data_dft = my_dt, col_v = 3, rmCol_v = F)
#' @export
convertDFT <- function(data_dft, col_v = NA, newName_v = NULL, rmCol_v = T, split_v = NULL) {
  
  ## Row names function
  addRowNames <- function(data_dft, out_dft, newName_v, split_v = NULL) {
    
    ### For no split
    if (is.null(split_v)) {
      
      newName_v <- ifelse(newName_v %in% colnames(out_dft), paste0(newName_v, "_2"), newName_v)
      out_dft[[newName_v]] <- rownames(data_dft)
      out_dft <- out_dft[, c(ncol(out_dft), 1:(ncol(out_dft)-1)), with = F]
      
    } else {
      
      for (i in 1:length(newName_v)) {
        
        currNewName_v <- newName_v[i]
        currNewName_v <- ifelse(currNewName_v %in% colnames(out_dft), paste0(currNewName_v, "_2"), currNewName_v)
        currNewVals_v <- sapply(rownames(data_dft), function(x) strsplit(x, split = split_v)[[1]][i])
        out_dft[,(currNewName_v) := currNewVals_v]
        
      } # for i
      
    } # fi split_v
    
    return(out_dft)
    
  } # addRowNames
  
  ## Get class
  class_v <- class(data_dft)
  
  ## New for R 4.x - matrix now has two classes - matrix and arrary
  ## So if matrix is in class, just ignore the array
  if (length(grep("matrix", class_v)) > 0) class_v <- "matrix"
  
  ## Convert data.table to data.frame
  if ("data.table" %in% class_v) {
    
    ## Convert
    out_dft <- as.data.frame(data_dft)
    
    ## Special case for multiple columns. Have to paste together and then remove the others
    if (length(col_v) > 1) {
      cols_v <- col_v
      collapse_v <- ifelse(is.null(split_v), "-_-", split_v)
      col_v <- paste(cols_v, collapse = collapse_v)
      out_dft[[col_v]] <- apply(out_dft[,c(cols_v)], 1, function(x) paste(x, collapse = collapse_v))
      for (c_v in cols_v) out_dft[[c_v]] <- NULL
    } else if (length(col_v) == 1) {
      col_v <- ifelse(is.na(col_v), colnames(data_dft)[1], 
                      ifelse(is.numeric(col_v), colnames(data_dft)[col_v], col_v))
    } else if (is.null(col_v)) {
      col_v <- NULL
    } else {
      stop("Didn't expect this to get triggered...check your if statements!")
    } # fi
    
    ## Add row names and handle column that provided names
    if (length(col_v) > 0) {
      
      rownames(out_dft) <- out_dft[[col_v]]
      
      ## Remove column that provided rownames
      if (rmCol_v) {
        whichCol_v <- which(colnames(out_dft) == col_v)
        out_dft <- out_dft[,-whichCol_v, drop = F]
      } # fi
      
    } # fi
    
    ## Convert data.frame to data.table
  } else if (class_v == "data.frame"){
    
    ## If no newName provided, assign V1
    ## If split_v, will split the name.
    if (is.null(newName_v)) {
      newName_v <- "V1"
    } else {
      if (!is.null(split_v)) {
        newName_v <- strsplit(newName_v, split = split_v)[[1]]
      } # fi
    } # fi
    
    ## Convert
    out_dft <- as.data.table(data_dft)
    
    ## Handle row names
    if (!identical(rownames(data_dft), as.character(1:nrow(data_dft)))) {
      out_dft <- addRowNames(data_dft, out_dft, newName_v, split_v = split_v)
    } # fi
    
    ## Convert matrix to data.table
  } else if (class_v == "matrix") {
    
    if (length(col_v) > 1) stop("Can only handle 1 col_v if input is matrix")
    
    ## Convert
    out_dft <- as.data.table(data_dft)
    
    ## Handle row names
    if (!is.null(rownames(data_dft))) {
      out_dft <- addRowNames(data_dft, out_dft, newName_v)
    } # fi
    
  } else {
    stop("Neither 'data.table', 'data.frame', nor 'matrix' were in the class of data_dft. Please check your input data.")
  } # fi
  
  ## Return
  return(out_dft)
  
} # convertDFT