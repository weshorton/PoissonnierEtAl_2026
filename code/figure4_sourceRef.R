####################################################################
######################## FIG4 SOURCE REF ###########################
####################################################################

library(data.table)

### Miscellaneous colors, column groups, functions, etc. used in the plotting/wrangling code for Figure 4

meta_dt <- data.table(Short = c("Veh", "PTX", "PLX", "PLX+PTX", "aPD-1", "PTX+aPD-1", 
                                "PLX+PTX+aPD-1", "aPD-L1", "PTX+aPD-L1", "PLX+PTX+aPD-L1"), 
                      Treat = c("Veh_IgG2b", "IgG1_PTX_IgG2a", "PLX+IgG2a", "PLX_PTX_IgG2a", "aPD-1", 
                                "AIN_PTX_aPD-1", "PLX_PTX_aPD-1", "aPD-L1", "PTX+aPD-L1", "PLX_PTX_aPD-L1"), 
                      Color = c("#0B0B09", "#F5EA30", "#EA8A2C", "#DB332A", "#7FCCE7", "#4A82B6",
                                "#DB3682", "#61B345", "#888B39", "#208643"))

colors_v <- meta_dt$Color
names(colors_v) <- meta_dt$Short

pctChangeLevels_v <- c("decrease: <0%", "slow: 0-25%", "med: 25-50%", "fast: 50-75%", "vfast: 75-100%", "mtd: >100%")
pctChangeColors_v <- rev(RColorBrewer::brewer.pal(11, "RdYlGn"))
pctChangeColors_v <- pctChangeColors_v[c(1, 2, 7, 8, 10,11)]
names(pctChangeColors_v) <- pctChangeLevels_v

###
### Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Tumor Waterfall
#' @param wide_dt wide tumor burden table
#' @param xVar_v variable to put on x-axis. Should be 'mouse'
#' @param treatCol_v variable to annotate treatments. Should be 'treatment'
#' @param yVar_v column that has y-variable. Should be in format Dstart_to_Dend_pctChange
#' @param addCounts_v logical to add the number of observations of each treatment above the bars
#' @param groupVar_v column that has y-var groupings. If NULL, is assumed to be gsub("pctChange", "growthRate", yVar_v)
#' @param groupColors_v this should be the level colors that are already made!
#' @param treatOrder_v optional treatment order. If NULL, treatments are assumed to be in correct order
#' @param treatColors_v optional treatment colors. If not provided, will use default ggplot colors
#' @param annotCol_v optional column identifying another variable to annotate plot by
#' @param annotOrder_v optional order of secondary annotation (not sure if this will work or not)
#' @param annotColors_v optional colors for the secondary annotation
#' @param addGaps_v logical indicating to add a spacer between treatments
#' @param title_v optional title
#' @param thinBars_v logical indicating to adjust bar theme to make them thinner
#' @param maxLabel_v set the y-value for treatment N labels. 'global' or 'individual' Passed to getTxLabelPositions
#' @return plot
#' @export
tumorWaterfall <- function(wide_dt, xVar_v = "mouse", treatCol_v = "treatment",
                           yVar_v = "D5_to_D14_pctChange", addCounts_v = T,
                           groupVar_v = NULL, groupColors_v = NULL,
                           treatOrder_v = NULL, treatColors_v = NULL, 
                           annotCol_v = NULL, annotOrder_v = NULL, annotColors_v = NULL,
                           addGaps_v = T, title_v = NULL, thinBars_v = F, maxLabel_v = "global") {
  
  ### Determine type of measurement
  isPct_v <- grepl("pctChange", yVar_v)
  isLog2_v <- grepl("log2Change", yVar_v)
  
  ### Get groupVar
  if (is.null(groupVar_v)) {
    if (isPct_v) {
      groupVar_v <- gsub("pctChange", "growthRate", yVar_v)
    } else if (isLog2_v) {
      groupVar_v <- gsub("log2Change", "log2growthRate", yVar_v)
    } else {
      stop("Bad groupVar_v substitution.")
    } # fi
  } # fi null groupVar
  
  ### Get colors
  if (is.null(groupColors_v)) {
    if (isPct_v) {
      groupColors_v <- pctChangeColors_v
    } else if (isLog2_v) {
      groupColors_v <- log2ChangeColors_v
    } # fi
  } # fi null
  
  ### Subset Data
  sub_dt <- wide_dt[!is.na(get(groupVar_v)),]
  
  ### Make sure x-variable is categorical
  sub_dt[[xVar_v]] <- as.character(sub_dt[[xVar_v]])
  
  ### Handle treatments
  if (is.null(treatOrder_v)) {
    if (is.factor(sub_dt[[treatCol_v]])) {
      treatOrder_v <- levels(sub_dt[[treatCol_v]])
    } else {
      treatOrder_v <- unique(sub_dt[[treatCol_v]])
    } # fi
  } else {
    sub_dt[[treatCol_v]] <- factor(as.character(sub_dt[[treatCol_v]]), levels = treatOrder_v)
  }
  
  ### Get treat colors
  if (is.null(treatColors_v)) {
    treatColors_v <- scales::hue_pal()(length(treatOrder_v))
    names(treatColors_v) <- treatOrder_v
  }
  
  ### Handle secondary annotation
  if (!is.null(annotCol_v)) {
    
    ### Order
    if (is.null(annotOrder_v)) {
      if (is.factor(sub_dt[[annotCol_v]])) {
        annotOrder_v <- levels(sub_dt[[annotCol_v]])
      } else {
        annotOrder_v <- unique(sub_dt[[annotCol_v]])
        sub_dt[[annotCol_v]] <- as.character(sub_dt[[annotCol_v]]) # make sure it's character
      } # fi is factor
    } else {
      sub_dt[[annotCol_v]] <- factor(as.character(sub_dt[[annotCol_v]]), levels = annotOrder_v)
    } # fi is null order
    
    ### Colors
    if (is.null(annotColors_v)) {
      annotColors_v <- scales::hue_pal(direction = -1)(length(annotOrder_v))
      names(annotColors_v) <- annotOrder_v
    } # fi is null colors
    
  } # fi is null annotCol
  
  ### Sort by treatment and value
  if (is.null(annotCol_v)) {
    setorderv(sub_dt, c(treatCol_v, yVar_v), order = c(1, -1))
  } else {
    setorderv(sub_dt, c(treatCol_v, annotCol_v, yVar_v), order = c(1, 1, -1))
  } # fi null annot
  
  ### Get table of treatments
  treats_dt <- as.data.table(table(sub_dt[[treatCol_v]]))
  
  ### Get count of each treatments and the x/y coordinates to put them
  sampleCounts_dt <- getTxLabelPositions(data_dt = sub_dt, treatCol_v = treatCol_v,
                                         treatOrder_v = treatOrder_v, xVar_v = xVar_v,
                                         yVar_v = yVar_v, groupVar_v = groupVar_v, 
                                         maxLabel_v = maxLabel_v)
  
  if (addGaps_v) {
    sub_dt <- addGaps(data_dt = sub_dt, treats_v = treatOrder_v,
                      treatCol_v = treatCol_v, xVar_v = xVar_v,
                      yVar_v = yVar_v, groupVar_v = groupVar_v,
                      annotCol_v = annotCol_v)
    ### Get white colors
    treatColors_v <- c(treatColors_v, "gap" = "white")
    sub_dt[grep("gap", get(xVar_v)), (treatCol_v) := "gap"]
    if (!is.null(annotCol_v)) annotColors_v <- c(annotColors_v, "gap" = "white")
  } # fi
  
  ### Factorize xvar to maintain order
  sub_dt[[xVar_v]] <- factor(sub_dt[[xVar_v]], levels = sub_dt[[xVar_v]])
  
  ### Make plot
  plot_gg <- ggplot(sub_dt, aes(x = !!sym(xVar_v), y = !!sym(yVar_v), fill = !!sym(groupVar_v))) +
    scale_fill_manual(values = groupColors_v, breaks = names(groupColors_v)) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 18), 
          plot.subtitle = element_text(hjust = 0.5, size = 14), 
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14))
    labs(x = xVar_v, y = "Percent Change", title = gsub("_pctChange", "", yVar_v))
  
  ### Make bars
  if (thinBars_v) {
    plot_gg <- plot_gg + geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 1))
  } else {
    plot_gg <- plot_gg + geom_bar(stat = "identity")
  } # fi thinBars
  
  ### Add Counts
  if (addCounts_v) plot_gg <- plot_gg + geom_text(data = sampleCounts_dt, aes(x = !!sym(xVar_v), y = !!sym(yVar_v), label = N))
  
  ### Add title
  if (!is.null(title_v)) plot_gg <- plot_gg + ggtitle(title_v)
  
  ### Annotate prep
  annot_lsv <- list("x" = treatCol_v); names(annot_lsv) <- treatCol_v
  annotColors_lsv <- list("x" = treatColors_v); names(annotColors_lsv) <- treatCol_v
  
  ### Add optional annotations
  if (!is.null(annotCol_v)) {
    annot_lsv[[annotCol_v]] <- annotCol_v
    annotColors_lsv[[annotCol_v]] <- annotColors_v
  }
  
  ### Run
  plot_lsgg <- annotatedBar(plot_gg = plot_gg, x_v = xVar_v,
                            annot_lsv = annot_lsv,
                            annotColors_lsv = annotColors_lsv)

  ### Output
  return(plot_lsgg)
  
} # tumorWaterfall

###################################################################################################

#' Calculate Change (percent or log)
#' @description
#' Calculate change between two columns
#' @param wide_dt wide tumor burden table
#' @param start_v name of first timepoint column
#' @param end_v name of second timepoint column
#' @param log_v logical indicating if we should calculate log-transformed change in addition to standard percent change.
#' @return wide_dt with new column of change
#' @details
#' Formula for regular percent change is ((end - start)/start) * 100
#' Formula for log change is log2(end/start)
#' Can give multiple start and end times, 1st start will go with first end, 2nd with second, etc.
#' They have to match lengths, unless one of them is just 1 value. e.g start_v = c(0, 3, 7), end_v = 14
#' will use Days 0, 3, and 7 as start dates and 14 as the end date for each.
#' start_v = c(0, 7) and end_v = c(7, 14) will do 0-7 and 7-14
#' start_v = c(0, 3, 7) and end_v = c(7, 14) will throw an error.
#' @export
calcChange <- function(wide_dt, start_v, end_v, log_v = F) {
  
  ### Handle lengths (repeat the length 1 variable to match longer variable)
  if (length(start_v) != length(end_v)) {
    if (length(start_v) == 1) {
      start_v <- rep(start_v, length(end_v))
    } else if (length(end_v) == 1) {
      end_v <- rep(end_v, length(start_v))
    } else {
      stop("start_v and end_v are of different lengths, but neither is length 1.\n")
    } # fi
  } # fi
  
  for (i in 1:length(start_v)) {
    
    currStart_v <- start_v[i]
    currEnd_v <- end_v[i]
    
    if (!currStart_v %in% colnames(wide_dt)) {
      warning(sprintf("%s was provided as a start day, but not found in data. Please check!\n", currStart_v))
      next
    }
    
    if (!currEnd_v %in% colnames(wide_dt)) {
      warning(sprintf("%s was provided as an end day, but not found in data. Please check!\n", currEnd_v))
      next
    }
    
    ### Calculate percent change
    wide_dt[, pctChange := (get(currEnd_v) - get(currStart_v))/get(currStart_v) * 100]
    
    ### Change the name - pctChange
    newName_v <- paste(currStart_v, "to", currEnd_v, "pctChange", sep = "_")
    colnames(wide_dt)[colnames(wide_dt) == "pctChange"] <- newName_v
    
    if (log_v) {
      
      ### Calculate log change
      wide_dt[, log2Change := log2(get(currEnd_v) / get(currStart_v))]
      
      ### Change the name - logChange
      newName_v <- paste(currStart_v, "to", currEnd_v, "log2Change", sep = "_")
      colnames(wide_dt)[colnames(wide_dt) == "log2Change"] <- newName_v
      
    } # fi log_v
    
  } # for i
  
  ### Output
  return(wide_dt)
  
} # calcChange

###################################################################################################

#' Classify Change
#' @description
#' Bin specified columns into different groups
#' @param data_dt table with various columns to classify
#' @param cols_v optional vector of column names to classify. If NULL (default) will run on all
#' @param which_v either 'pct' or 'log2' indicating which change columns to use. SEE DETAILS
#' @param groups_lsv list indicating how to bin. See details.
#' @details Each element name of groups_lsv is a distinct bin. The values
#' for each element are the lower and upper limits, except for the first and last
#' The upper limit is inclusive (i.e. <=), while the lower is exclusive (i.e. <)
#' have to provide the which_v argument because the groups_lsv info is only valid for pct OR log2, not both
#' potentially update this to handle either at the same time, but not right now.
#' @return data.table with same number of rows as input and one more column per either cols_v value or columns
#' of D#_to_D# format
#' @export
classifyChange <- function(data_dt,
                           cols_v = NULL,
                           which_v = "pct",
                           groups_lsv = list("decrease: <0%" = 0,
                                             "slow: 0-25%" = c(0, 25),
                                             "med: 25-50%" = c(25, 50),
                                             "fast: 50-75%" = c(50, 75),
                                             "vfast: 75-100%" = c(75, 100),
                                             "mtd: >100%" = 100)) {
  
  ### Get columns
  if (is.null(cols_v)) {
    cols_v <- grep("D[0-9]+_to_D[0-9]+", colnames(data_dt), value = T)
    if (length(cols_v) == 0) stop("cols_v argument is NULL, but no columns of format 'D[0-9]+_to_D[0-9]+' found. Please check!\n")
  } # fi
  cols_v <- grep(which_v, cols_v, value = T)
  
  ### Classify
  for (col_v in cols_v) {
    
    new_v <- gsub("Change", "Rate", col_v)
    data_dt[[new_v]] <- classifyColumn(rate_v = data_dt[[col_v]],
                                       groups_lsv = groups_lsv)
  } # for i
  
  return(data_dt)
  
} # classifyChange

###################################################################################################

#' Classify Column
#' @description
#' Bin percent changes into different groups
#' @param rate_v vector of percent changes to bin
#' @param groups_lsv list indicating how to bin. See details.
#' @details Each element name of groups_lsv is a distinct bin. The values
#' for each element are the lower and upper limits, except for the first and last
#' The upper limit is inclusive (i.e. <=), while the lower is exclusive (i.e. <)
#' @return vector of same length as rate_v, but as categories
#' @export
classifyColumn <- function(rate_v, 
                           groups_lsv = list("decrease: <0%" = 0,
                                             "slow: 0-25%" = c(0, 25),
                                             "med: 25-50%" = c(25, 50),
                                             "fast: 50-75%" = c(50, 75),
                                             "vfast: 75-100%" = c(75, 100),
                                             "mtd: >100%" = 100)) {
  
  tmp_dt <- data.table("rate" = rate_v, "bin" = character(length = length(rate_v)))
  for (i in 1:length(groups_lsv)) {
    
    binName_v <- names(groups_lsv)[i]
    binVals_v <- groups_lsv[[binName_v]]
    
    if (i == 1) {
      tmp_dt[rate <= binVals_v[1], bin := binName_v]
    } else if (i < length(groups_lsv)) {
      tmp_dt[(rate > binVals_v[1] & rate <= binVals_v[2]), bin := binName_v]
    } else {
      tmp_dt[rate > binVals_v[1], bin := binName_v]
    } # fi
    
  } # for i
  
  out_v <- tmp_dt$bin
  out_v <- factor(out_v, levels = names(groups_lsv))
  return(out_v)
  
} # classifyColumn

###################################################################################################

#' Extract ggplot legend
#' @description Extract legend as separate gtable from ggplot object
#' @param a.gplot a ggplot object with a legend
#' @return a gtable object of the legend
#' @export
g_legend <- function(a.gplot){
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
    return(legend)
  } else {
    message("No legend grob in provided ggplot")
  }
} # g_legend

#' Annotated Bar Plot
#' @description Standard bar-plot with heatmap-like x-axis annotations
#' @param data_dt melted data table for plotting (see details. Must supply this or plot_gg)
#' @param plot_gg ggplot to annotate (see details. Must supply this or data_dt)
#' @param x_v name of x-axis variable (usually Sample_ID, Treatment, etc.). Must be column in data_dt
#' @param y_v name of y-axis variable (onnly used if stat_v = "identity")
#' @param stat_v argument to 'stat' parameter of geom_bar. 'count' is default and does not use y_v. "identity" requires y_v to be set
#' @param position_v argument to 'position' parameter of geom_bar. 'stack' (default) is counts; 'fill' turns to percentage by filling out of 1.
#' @param yLab_v name of y-axis label.
#' @param fill_v name of fill variable. Must be column in data_dt (only used if data_dt used)
#' @param fillColors_v optional named color vector. Values are colors, names are values of data_dt[[fill_v]]
#' @param keepXAxis_v logical indicating if x-axis labels should be included. Default is false. Should only be used if there aren't too many x-axis values
#' @param legendPos_v character vector indicating where to place legends. Either below plot ("bottom"), or to the right of the plot ("right")
#' @param title_v optional plot title
#' @param annot_lsv list of annotations to add below the plot. e.g. annot_lsv = list("outName1" = "colName1", "outName2" = "colName2")
#' @param annotColors_lsv list of colors for annotations. list element names are same as annot_lsv names; list elements are named color vectors, names are values of annot_lsv colNames in data_dt
#' @param backgroundCol_v sent to panel.background element of ggplot theme. Needs to be white for a specific use case where I want to plot a percentage using fill_v = "fill", but I only want a subset of the identities to be shown.
#' @details make a standard ggplot bar plot with extra annotations. Can provide a pre-existing plot instead of data.
#' If a plot is provided, the x_v variable is still required, along with annot_lsv and annotColors_lsv. 
#' Not required if plot is given: stat_v, position_v, fill_v, fillColors_v
#' Note that data_dt overrules plot_gg. If both are provided, plot_gg is ignored
#' @return list of gg objects output by ggarrange. Combo plot, data only, annotations only.
#' @export
annotatedBar <- function(data_dt = NULL, plot_gg = NULL, x_v, y_v = NULL, 
                         stat_v = "count", position_v = "stack", yLab_v = "count", 
                         fill_v = NULL, fillColors_v = NULL, keepXAxis_v = F, legendPos_v = "bottom",
                         title_v = NULL, annot_lsv = NULL, annotColors_lsv, backgroundCol_v = "white") {
  
  if (!is.null(data_dt)) {
    
    if (is.null(fill_v)) stop("Need to provide fill column.")
    
    ### Base plot depends on stat argument
    if (stat_v == "identity") {
      if (is.null(y_v)) stop("Must provide y_v argument if stat is identity.\n")
      plot_gg <- ggplot(data = data_dt, aes(x = !!sym(x_v), y = !!sym(y_v), fill = !!sym(fill_v)))
    } else if (stat_v == "count") {
      plot_gg <-  ggplot(data = data_dt, aes(x = !!sym(x_v), fill = !!sym(fill_v)))
    } else {
      stop("Can only handle 'identity' or 'count' for stat_v argument")
    } # fi
    
    ### Add the rest to the base plot
    plot_gg <- plot_gg +
      geom_bar(stat = stat_v, position = position_v) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + # make 0 the bottom
      theme(legend.position = "bottom", 
            plot.title = element_text(size = 26, hjust = 0.5),
            axis.text = element_text(size = 16),
            panel.background = element_rect(fill = backgroundCol_v)) +
      ylab(yLab_v)
    
    ### Add title
    if (!is.null(title_v)) plot_gg <- plot_gg + ggtitle(title_v)
    
    ### Add colors
    if (!is.null(fillColors_v)) plot_gg <- plot_gg + scale_fill_manual(values = fillColors_v, breaks = names(fillColors_v))
    
  } else if (!is.null(plot_gg)) {
    
    ### Extract data from plot
    data_dt <- as.data.table(plot_gg$data)
    
    ### Factorize x-axis to maintain order
    if (!is.factor(data_dt[[x_v]])) {
      data_dt[[x_v]] <- factor(data_dt[[x_v]], levels = unique(data_dt[[x_v]]))
    }
    
  } else {
    
    stop("Must provide either data_dt or plot_gg.\n")
    
  } # fi !is.null(data_dt)
  
  ### Begin annotations
  if (is.null(annot_lsv)) {
    warning("annot_lsv not provided...just make a regular barplot if you don't have annotations.\n")
    return(plot_gg)
  } else {
    
    ### Extract legend
    mainLegend_gg <- g_legend(plot_gg + theme(legend.text = element_text(size = 18), legend.title = element_text(size = 20),
                                              plot.margin = unit(c(t = 0, r = 1, b = 1, l = 1), "cm")))
    plot_gg <- plot_gg + theme(legend.position = "none")

    ### Extract axis
    mainXAxis_gg <- ggpubr::as_ggplot(cowplot::get_x_axis(plot_gg + theme(plot.margin = unit(c(t = -0.1, r = 1, b = 1, l = 1), "cm"))))
    plot_gg <- plot_gg + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

    ### Set margin
    plot_gg <- plot_gg + theme(plot.margin = unit(c(t = 1, r = 1, b = 0, l = 1), "cm"))

    ### Build bar charts
    bar_lsgg <- barLegend_lsgg <- list()
    for (i in 1:length(annot_lsv)) {

      ### Get info
      currName_v <- names(annot_lsv)[i]
      currColName_v <- annot_lsv[[currName_v]]

      ### Make bar
      currBar_gg <- suppressWarnings(ggplot(data_dt, aes(x = !!sym(x_v), y = 1, fill = !!sym(currColName_v)))) +
        geom_bar(stat = "identity", width = 1) +
        scale_y_continuous(limits = c(0,1)) +
        theme(panel.background = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "bottom",
              panel.grid = element_blank()) +
        scale_fill_manual(values = annotColors_lsv[[currName_v]], breaks = names(annotColors_lsv[[currName_v]])) +
        labs(fill = currName_v)

      ### Split legend
      currBarLegend_gg <- suppressWarnings(g_legend(currBar_gg))
      currBar_gg <- currBar_gg + theme(legend.position = "none")

      ### Set margin
      currBar_gg <- currBar_gg + theme(plot.margin = unit(c(t = -0.1, r = 1, b = -0.1, l = 1), "cm"))

      ### Save
      bar_lsgg[[currName_v]] <- currBar_gg
      barLegend_lsgg[[currName_v]] <- currBarLegend_gg

    } # for i

    ### New legend version - this is messed up. Doesn't work with the scale fill manual...
    ### As it is now, annot_lsv has to be list("Name" = "Name") and annotColors_lsv must be list("Name" = colorVector_v)
    tmp <- melt(data_dt[,mget(c(x_v, names(annot_lsv)))], measure.vars = names(annot_lsv))
    tmpColor <- unlist(annotColors_lsv)
    names(tmpColor) <- gsub("^.*\\.", "", names(tmpColor))
    ### This is new...need to check with clinical3x plots and tumorBurden to make sure works for both situations.
    if (is.logical(all.equal(tmpColor[names(annotColors_lsv)], NA_character_))) {
      tmpColor <- tmpColor[names(annotColors_lsv)]
      comboBarLegend_gg <- ggplot(tmp, aes(x = !!sym(x_v), y = 1, fill = variable)) + geom_boxplot() +
        scale_fill_manual(values = tmpColor, breaks = names(tmpColor)) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              plot.subtitle = element_text(hjust = 0.5, size = 14),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14))
    } else {
      comboBarLegend_gg <- ggplot(tmp, aes(x = !!sym(x_v), y = 1, fill = value)) + geom_boxplot() +
        scale_fill_manual(values = tmpColor, breaks = names(tmpColor)) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              plot.subtitle = element_text(hjust = 0.5, size = 14),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14))
    }

    comboBarLegend_gg <- comboBarLegend_gg + theme(legend.margin = margin(t = 0, r = 0.5, b = 0, l = 0.5, "cm"))
    comboBarLegend_gg <- get_legend(comboBarLegend_gg)

    ### Add axis
    if (keepXAxis_v) {
      toPlot_lsgg <- c(list("plot" = plot_gg), bar_lsgg, list("axis" = mainXAxis_gg))
      heights_v <- c(2, rep(0.1, length(bar_lsgg)), 1)
    } else {
      toPlot_lsgg <- c(list("plot" = plot_gg), bar_lsgg)
      heights_v <- c(2, rep(0.1, length(bar_lsgg)))
    }

    # toPlotLegend_lsgg <- c(list("plot" = mainLegend_gg), barLegend_lsgg)
    toPlotLegend_lsgg <- c(list("plot" = mainLegend_gg, "annot" = comboBarLegend_gg))

    if (length(bar_lsgg) > 3) heights_v[1] <- 1.9

    ### Combine plots
    combo_gg <- suppressWarnings(ggpubr::ggarrange(plotlist = toPlot_lsgg, ncol = 1, nrow = length(toPlot_lsgg), align = "v", heights = heights_v))

    ### Combine legends
    if (length(toPlotLegend_lsgg) > 3) {nrow_v <- 2; ncol_v <- ceiling(length(toPlotLegend_lsgg)/2)} else {nrow_v <- 1; ncol_v <- length(toPlotLegend_lsgg)}
    comboLegend_gg <- ggpubr::ggarrange(plotlist = toPlotLegend_lsgg, nrow = nrow_v, ncol = ncol_v)

    ### Combine both
    if (legendPos_v == "bottom") {
      finalCombo_gg <- ggpubr::ggarrange(plotlist = list(combo_gg, comboLegend_gg), ncol = 1, nrow = 2, heights = c(2, 1))
    } else if (legendPos_v == "right") {
      finalCombo_gg <- ggpubr::ggarrange(plotlist = list(combo_gg, comboLegend_gg), ncol = 2, nrow = 1, widths = c(2, 1))
    } else {
      stop("Only 'bottom' and 'right' are available for legend position.")
    }

    ### Make output
    out_ls <- list("annotPlot" = finalCombo_gg, "plot" = combo_gg, "legend" = comboLegend_gg)

    return(out_ls)
    
  } # fi is.null(annot_lsv)
} # annotatedBar

###################################################################################################

#' Get Tx Label Positions
#' @description
#' Get x,y positions for labels
#' @param data_dt data.table
#' @param treatCol_v variable to annotate treatments. Should be 'treatment'
#' @param treatOrder_v treatment order (should be made already in plot function)
#' @param xVar_v variable to put on x-axis. Should be 'mouse'
#' @param yVar_v column that has y-variable. Should be in format Dstart_to_Dend_pctChange
#' @param groupVar_v column that has y-var groupings. If NULL, is assumed to be gsub("pctChange", "growthRate", yVar_v)
#' @param maxLabel_v determines y-axis of labels. either 'global' to reference the global maximum value, or "individual" to reference each treatment
#' @return table of labels and x,y coordinates
#' @export
getTxLabelPositions <- function(data_dt, treatCol_v, treatOrder_v, 
                                xVar_v, yVar_v, groupVar_v, maxLabel_v = "global") {
  
  treats_dt <- as.data.table(table(data_dt[[treatCol_v]]))
  
  ### Global max
  max_v <- max(data_dt[[yVar_v]]) * .9
  
  out_lsdt <- list()
  for (i in 1:length(treatOrder_v)) {
    
    ### Get a treatment
    currTreat_v <- treatOrder_v[i]
    
    ### Subset for just those entries
    treat_dt <- data_dt[get(treatCol_v) == currTreat_v,]
    
    ### Get the middle entry
    middle_dt <- treat_dt[ceiling(nrow(treat_dt)/2),]
    
    ### Increase its y value 1.5x (to make the text above the data)
    ### Or use the global max
    if (maxLabel_v == "global") {
      middle_dt[[yVar_v]] <- max_v
    } else if (maxLabel_v == "individual") {
      middle_dt[[yVar_v]] <- middle_dt[[yVar_v]] * 1.5
    } else {
      stop(sprintf("max_v argument can be 'global' or 'individual'. %s provided.\n", maxLabel_v))
    } # fi
    
    ### Subset
    middle_dt <- middle_dt[,mget(c(xVar_v, treatCol_v, groupVar_v, yVar_v))]
    
    ### Add count
    middle_dt$N <- treats_dt[V1 == currTreat_v, N]
    
    ### Add to list
    out_lsdt[[currTreat_v]] <- middle_dt
  } # for i
  
  out_dt <- do.call(rbind, out_lsdt)
  return(out_dt)
  
} # getTxLabelPositions

###################################################################################################

#' Add Gaps
#' @description
#' Add a spacer between groups in the chart to make it easier to read
#' @param data_dt data to plot for waterfall (already wrangled)
#' @param treats_v should be treatOrder_v arg from tumorWaterfall()
#' @param treatCol_v should be treatCol_v arg from tumorWaterfall() (treatment)
#' @param xVar_v should be xVar_v arg (mouse)
#' @param yVar_v yVar that needs 0 padded. (yVar_v for tumorWaterfall() and 'pctChange' for multiDate())
#' @param groupVar_v should be tumorWaterfall()'s groupVar_v or 'growthRate' for multiDate()
#' @param annotCol_v optional additional annotation, passed from main function
#' @return input data_dt with one empty-value padder row after each treatment
#' @export
addGaps <- function(data_dt, treats_v, treatCol_v, xVar_v, yVar_v, groupVar_v, annotCol_v = NULL) {
  
  ### Empty list
  tmp_lsdt <- list()
  
  ### Loop for each treat
  for (t_v in treats_v) {
    
    ### Get treatment only
    curr_dt <- data_dt[get(treatCol_v) == t_v,]
    
    ### Add first row back to the bottom
    curr_dt <- do.call(rbind, list(curr_dt, curr_dt[1,]))
    
    ### Change x-var to be unique
    curr_dt[nrow(curr_dt), (xVar_v) := paste0("gap_", t_v)]
    
    ### Make y-var 0 so it doesn't show
    curr_dt[nrow(curr_dt), (yVar_v) := 0]
    
    ### Replace group var value. (Don't think this does anything)
    curr_dt[nrow(curr_dt), (groupVar_v) := curr_dt[[groupVar_v]][1]]
    
    ### Replace annot column
    if (!is.null(annotCol_v)) curr_dt[nrow(curr_dt), (annotCol_v) := "gap"]
    
    ### Add back to list
    tmp_lsdt[[t_v]] <- curr_dt
    
  } # for t_v
  
  ### Re-bind
  data_dt <- do.call(rbind, tmp_lsdt)
  
  ### Remove last gap
  data_dt <- data_dt[1:(nrow(data_dt)-1),]
  
  ### Output
  return(data_dt)
  
} # addGaps

###################################################################################################

#' Simple Stats
#' @description
#' Run simple one-variable stats to compare group means
#' @param data_dt data.table
#' @param depVar_v dependent variable - numeric
#' @param indVar_v independent variable. If not factored, use indOrder_v
#' @param indOrder_v optional order for independent variable
#' @param compType_v character indicating which statistical test to use
#' @details
#' Can do either ANOVA with Tukey's post-hoc ('anova'), 
#' Kruskal-Wallis with Dunn's test (KW), or Wilcoxon Rank-Sum (wilcox)
#' with multiple testing correction
#' @return data.table with results
#' @export
simpleStats <- function(data_dt, depVar_v, indVar_v, indOrder_v = NULL, compType_v) {
  
  ### Handle independent variable order
  if (is.null(indOrder_v)) {
    if (!is.factor(data_dt[[indVar_v]])) {
      data_dt[[indVar_v]] <- factor(data_dt[[indVar_v]], 
                                    levels = unique(data_dt[[indVar_v]]))
    } # fi factor
  } else {
    data_dt[[indVar_v]] <- factor(as.character(data_dt[[indVar_v]]), levels = indOrder_v)
  } # fi null order
  
  ### Build formula
  formula_v <- as.formula(paste(depVar_v, indVar_v, sep = " ~ "))
  
  ### Run tests
  if (compType_v %in% c("a", "anova")) {
    
    ### Run test
    aovRes <- aov(formula_v, data_dt)
    aovSummary <- summary(aovRes)
    tukeyRes <- TukeyHSD(aovRes)
    
    ### Reformat
    tukeyRes_dt <- convertDFT(tukeyRes$Short, newName_v = "Comparison")
    aovSummary_dt <- data.table("Comparison" = "ANOVA", "diff" = aovSummary[[1]][1,4],
                                "lwr" = NA, "upr" = NA, "p adj" = aovSummary[[1]][1,5])
    out_dt <- rbind(aovSummary_dt, tukeyRes_dt)
    colnames(out_dt)[colnames(out_dt) == "p adj"] <- "pAdj"
    
  } else if (compType_v %in% c("k", "K", "KW")) {
    
    ### Run Test
    KW <- kruskal.test(formula_v, data_dt)
    dunn <- FSA::dunnTest(formula_v, data = data_dt, method = 'bh')
    
    ### Reformat
    D_df <- dunn$res
    D_df <- rbind(c("Kruskal-Wallis: all", KW$statistic, KW$p.value, NA), D_df)
    out_dt <- as.data.table(D_df)
    numCols_v <- c("Z", "P.unadj", "P.adj")
    out_dt[, (numCols_v) := lapply(.SD, as.numeric), .SDcols = numCols_v]
    colnames(out_dt)[colnames(out_dt) == "P.adj"] <- "pAdj"
    
  } else if (compType_v %in% c("w", "W", "wilcox")) {
    
    ### Get pairwise combinations of treatments
    txCombo_lsv <- combn(levels(data_dt[[indVar_v]]), 2, simplify = F)
    wilcoxOut_mat <- do.call(rbind, txCombo_lsv)
    
    ### Run wilcox test for each pairwise comp
    wilcox_lsw <- lapply(txCombo_lsv, function(x) {
      wilcox.test(formula_v, data = data_dt[get(indVar_v) %in% x,])})
    
    ### Extract p-value and adjust it
    WP_v <- unlist(lapply(wilcox_lsw, function(x) x$p.value))
    WPadj_v <- p.adjust(WP_v, method = 'BH')
    #WPadj_v <- p.adjust(WP_v, method = 'hochberg')
    
    ### Combine with treatment names and convert
    wilcoxOut_mat <- cbind(wilcoxOut_mat, WP_v, WPadj_v)
    out_dt <- as.data.table(wilcoxOut_mat)
    colnames(out_dt) <- c("Treat1", "Treat2", "P-value", "pAdj")
    
    ### Convert numeric columns
    numCols_v <- c("P-value", "pAdj")
    out_dt[, (numCols_v) := lapply(.SD, as.numeric), .SDcols = numCols_v]
    
    ### Add comparison column
    out_dt[, Comparison := paste(Treat2, Treat1, sep = " - ")]
    
  } else {
    stop(sprintf("compType_v can only be 'anova', 'KW', or 'wilcox', %s provided.", compType_v))
  } # run tests
  
  ### Output
  return(out_dt)
  
} # simpleStats