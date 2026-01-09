####################################################################
######################## FIG4 SOURCE REF ###########################
####################################################################

### Miscellaneous colors, column groups, functions, etc. used in the plotting/wrangling code for Figure 7

###
### Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Read Directory
#' @description Read in all of the files in a given directory
#' @param dir_v Path to directory containing files that should be read in
#' @param files_v Vector containing the files within dir_v to read. See details for more info.
#' @param pattern_v Regex pattern to search for files within dir_v to read. See details for more info.
#' @param names_v Vector containing names to assign to each element. If blank, will use file name minus extension.
#' @details There are multiple ways to determine which files will be read. First, files_v will be checked. If
#' this argument has been specified, then these are the files that will be read in. If it is blank (NA), then
#' pattern_v will be checked. If this has been specified, then it will be used to find matching files. If both
#' files_v and pattern_v are blank, then all of the files in the directory will be read.
#' @return List of data.tables
#' @export
readDir <- function(dir_v, files_v = NULL, pattern_v = NULL, names_v = NULL) {
  
  ### File determination
  if (is.null(files_v)) files_v <- list.files(dir_v, pattern = pattern_v)
  
  ### Read in
  data_lsdt <- lapply(files_v, function(x) {
    ext_v <- tools::file_ext(x)
    if (ext_v %in% c("xls", "xlsx")) {
      readAllExcel(file.path(dir_v, x))
    } else {
      fread(file.path(dir_v, x))
    } # fi
  })
  
  ### Add names
  if (is.null(names_v)) names_v <- tools::file_path_sans_ext(files_v)
  names(data_lsdt) <- names_v
  
  ### Return
  return(data_lsdt)
} # readDir

###################################################################################################

#' Get Shared Clones
#' @description Get shared clonotypes using V gene, J gene, and amino acid CDR3
#' @param clones_dt properly subset and rbound clone table
#' @return table of shared clones only
#' @export
getShared <- function(clones_dt) {
  
  ### New column combining V, CDR3, and J
  clones_dt$V_aaSeq_J <- paste(clones_dt$V.name, clones_dt$CDR3.aa, clones_dt$J.name, sep = "_")
  
  ### Table
  shared_dt <- as.data.table(table(clones_dt$V_aaSeq_J))[N > 1,]
  
  ### Re-split columns
  shared_dt$V.name <- sapply(shared_dt$V1, function(x) strsplit(x, split = "_")[[1]][1])
  shared_dt$J.name <- sapply(shared_dt$V1, function(x) strsplit(x, split = "_")[[1]][3])
  shared_dt$CDR3.aa <- sapply(shared_dt$V1, function(x) strsplit(x, split = "_")[[1]][2])
  shared_dt <- shared_dt[,mget(c("V.name", "J.name", "CDR3.aa"))]
  
  ### Output
  return(shared_dt)
  
} # getShared

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
#' @return either a data.table or data.frame (opposite class of input)
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
        #out_dft[[currNewName_v]] <- sapply(rownames(data_dft), function(x) strsplit(x, split = split_v)[[1]][i])
        
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