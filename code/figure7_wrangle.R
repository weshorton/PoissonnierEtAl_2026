##################################################################
######################## FIGURE 7 WRANGLE ########################
##################################################################

### Set to "/" for CodeOcean and "./" for github
baseDir_v <- "./"

### Read in data used for main figure 7 and wrangle for plots

source(file.path(baseDir_v, "code/figure7_sourceRef.R"))
library(data.table)
library(ggplot2)
library(immunarch)
dir.create(file.path(baseDir_v, "data/fig7/rds"))

###
### Read in ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Tumor
tumorDiv_dt <- fread(file.path(baseDir_v, "data/fig7/diversity/tumor_divAnalysis.txt"))
tumorMeta_dt <- fread(file.path(baseDir_v, "data/fig7/meta/Tumor_meta.txt"))
geoTumorMeta_dt <- fread(file.path(baseDir_v, "data/fig7/meta/Tumor_samples.csv"))
tumorClones_lsdt <- readDir(file.path(baseDir_v, "data/fig7/clones/tumor/"))
tumorFreqGroups_dt <- fread(file.path(baseDir_v, "data/fig7/cloneDiv/Tumor_cumClonalFreqs.txt"))
tumorTopClones_dt <- fread(file.path(baseDir_v, "data/fig7/topClones/Tumor_topCloneFreq.txt"))

### Lung
lungDiv_dt <- fread(file.path(baseDir_v, "data/fig7/diversity/lung_divAnalysis.txt"))
lungMeta_dt <- fread(file.path(baseDir_v, "data/fig7/meta/Lung_meta.txt"))
geoLungMeta_dt <- fread(file.path(baseDir_v, "data/fig7/meta/Lung_samples.csv"))
lungClones_lsdt <- readDir(file.path(baseDir_v, "data/fig7/clones/lung/"))
lungFreqGroups_dt <- fread(file.path(baseDir_v, "data/fig7/cloneDiv/Lung_cumClonalFreqs.txt"))
lungTopClones_dt <- fread(file.path(baseDir_v, "data/fig7/topClones/Lung_topCloneFreq.txt"))

### Blood
bloodDiv_dt <- fread(file.path(baseDir_v, "data/fig7/diversity/blood_divAnalysis.txt"))
bloodMeta_dt <- fread(file.path(baseDir_v, "data/fig7/meta/Blood_meta.txt"))
geoBloodMeta_dt <- fread(file.path(baseDir_v, "data/fig7/meta/Blood_samples.csv"))
bloodClones_lsdt <- readDir(file.path(baseDir_v, "data/fig7/clones/blood/"))
bloodFreqGroups_dt <- fread(file.path(baseDir_v, "data/fig7/cloneDiv/Blood_cumClonalFreqs.txt"))
bloodTopClones_dt <- fread(file.path(baseDir_v, "data/fig7/topClones/Blood_topCloneFreq.txt"))

###
### Plot info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Tumor
tumorTitleLevels_v <- c("d80 WTNT", "d80 Untreated", 
                        "d105 WTNT", "d105 Untreated", 
                        "d105 PTX", "d105 PLX",
                        "d105 PLX+PTX", "d105 PTX+aPD-1", 
                        "d105 PLX+PTX+aPD-1_R", "d105 PLX+PTX+aPD-1_NR", 
                        "d105 aPD-L1", "d105 PTX+aPD-L1", "d105 PLX+PTX+aPD-L1",
                        "d95 PTX", "d95 PLX+PTX+aPD-1_R", "d95 PLX+PTX+aPD-1_NR")

### Lung
lungTitleLevels_v <- c("d80 WTNT", "d80 Untreated", 
                       "d105 WTNT", "d105 Untreated", 
                       "d105 PTX", "d105 PLX",
                       "d105 PLX+PTX", "d105 PTX+aPD-1", 
                       "d105 PLX+PTX+aPD-1_R", "d105 PLX+PTX+aPD-1_NR", 
                       "d105 aPD-L1", "d105 PTX+aPD-L1", "d105 PLX+PTX+aPD-L1",
                       "d95 PTX", "d95 PLX+PTX", "d95 PLX+PTX+aPD-1_R", "d95 PLX+PTX+aPD-1_NR")

### Blood
bloodTitleLevels_v <- c("d105 WTNT", "d95 Untreated",
                        "d105 Untreated", "d105 PTX",
                        "d105 PLX+PTX", "d105 PTX+aPD-1",
                        "d105 PLX+PTX+aPD-1_R", "d105 PLX+PTX+aPD-1_NR",
                        "d95 PTX", "d95 PLX+PTX",
                        "d95 PLX+PTX+aPD-1_R", "d95 PLX+PTX+aPD-1_NR")

### Jaccard levels
treatCompLevels_v <- c("PLX+PTX+aPD-1_R_PLX+PTX+aPD-1_NR", 
                       "PLX+PTX+aPD-1_R_PLX+PTX+aPD-L1",
                       "PLX+PTX+aPD-1_NR_PLX+PTX+aPD-L1")

### Colors
colors_dt <- setDT(readxl::read_xlsx(file.path(baseDir_v, "data/fig7/colorCodes.xlsx")))
colors_v <- colors_dt$Hex
names(colors_v) <- colors_dt$Treatment

### Clone colors
cloneColors_dt <- setDT(readxl::read_xlsx(file.path(baseDir_v, "data/fig7/cloneColors.xlsx")))
cloneColors_v <- cloneColors_dt$Hex
names(cloneColors_v) <- cloneColors_dt$Clonotype

### Top Clone Colors
topCloneColors_v <- c("#1B4B7E", "#7ECCE7", "#F5EA31", "#DA322A", "#21A849")
names(topCloneColors_v) <- c("5", "6-25", "26-50", "51-100", "rest")

### Freq Group Colors
freqGroupColors_v <- c("#292E6E", "#7ECCE7", "#F9EF7D", "#E57570", "#158643")
names(freqGroupColors_v) <- c("Rare", "Small", "Medium", "Large", "Hyperexpanded")

### Plot Levels
plotInfo_lsv <- list("diversity" = c("Shannon Entropy", "Clonality"),
                     "freqGroups" = names(freqGroupColors_v),
                     "freqGroupColors" = freqGroupColors_v,
                     "topClones" = names(topCloneColors_v),
                     "topCloneColors" = topCloneColors_v,
                     "treatColors" = colors_v,
                     "trackCloneColors" = cloneColors_v)
saveRDS(plotInfo_lsv, file = file.path(baseDir_v, "data/fig7/rds/plotInfo_lsv.rds"))

###
### Diversity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

divOut_lsdt <- list()

###
### Tumor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove file extension and subset columns
tumorDiv_dt$File <- gsub("_[a-z]+.*$|\\..*$", "", tumorDiv_dt$File)
tumorDiv_dt <- tumorDiv_dt[,mget(c("File", "Shannon Entropy", "Clonality"))]

### Update geo sample and subset columns
geoTumorMeta_dt[, Sample := paste(Batch, Sample, sep = "_")]
idCols_v <- c("Sample", "geoSample", "Treatment", "Experiment", "Title")
geoTumorMeta_dt <- geoTumorMeta_dt[,mget(idCols_v)]

### Merge
tumorDiv_dt <- merge(geoTumorMeta_dt, tumorDiv_dt, by.x = "Sample", by.y = "File", sort = F)

### Factorize
tumorDiv_dt$Title <- factor(tumorDiv_dt$Title, levels = tumorTitleLevels_v)

### Sort
setorder(tumorDiv_dt, Title)

### Add to list
divOut_lsdt[["Tumor"]] <- tumorDiv_dt

###
### Lung ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove file extension and subset columns
lungDiv_dt$File <- gsub("_[a-z]+.*$|\\..*$", "", lungDiv_dt$File)
lungDiv_dt <- lungDiv_dt[,mget(c("File", "Shannon Entropy", "Clonality"))]

### Update geo sample and subset columns
geoLungMeta_dt[, Sample := paste(Batch, Sample, sep = "_")]
idCols_v <- c("Sample", "geoSample", "Treatment", "Experiment", "Title")
geoLungMeta_dt <- geoLungMeta_dt[,mget(idCols_v)]

### Merge
lungDiv_dt <- merge(geoLungMeta_dt, lungDiv_dt, by.x = "Sample", by.y = "File", sort = F)

### Factorize
lungDiv_dt$Title <- factor(lungDiv_dt$Title, levels = lungTitleLevels_v)

### Sort
setorder(lungDiv_dt, Title)

### Add to list
divOut_lsdt[["Lung"]] <- lungDiv_dt

###
### Blood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove file extension and subset columns
bloodDiv_dt$File <- gsub("_[a-z]+.*$|\\..*$", "", bloodDiv_dt$File)
bloodDiv_dt <- bloodDiv_dt[,mget(c("File", "Shannon Entropy", "Clonality"))]

### Update geo sample and subset columns
geoBloodMeta_dt[, Sample := paste(Batch, Sample, sep = "_")]
idCols_v <- c("Sample", "geoSample", "Treatment", "Experiment", "Title")
geoBloodMeta_dt <- geoBloodMeta_dt[,mget(idCols_v)]

### Merge
bloodDiv_dt <- merge(geoBloodMeta_dt, bloodDiv_dt, by.x = "Sample", by.y = "File", sort = F)

### Factorize
bloodDiv_dt$Title <- factor(bloodDiv_dt$Title, levels = bloodTitleLevels_v)

### Sort
setorder(bloodDiv_dt, Title)

### Add to list
divOut_lsdt[["Blood"]] <- bloodDiv_dt

###
### Freq Groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

freqGroupOut_lsdt <- list()

###
### Tumor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Merge with meta
tumorFreqGroups_dt <- merge(geoTumorMeta_dt, tumorFreqGroups_dt, by = "Sample", sort = F)

### Factorize
tumorFreqGroups_dt$Title <- factor(tumorFreqGroups_dt$Title, levels = tumorTitleLevels_v)

### Make numeric
for (col in plotInfo_lsv$freqGroups) set(tumorFreqGroups_dt, j = col, value = as.numeric(tumorFreqGroups_dt[[col]]))

### Set order
setorder(tumorFreqGroups_dt, Title)

### Add to list
freqGroupOut_lsdt[["Tumor"]] <- tumorFreqGroups_dt

###
### Lung ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Merge with meta
lungFreqGroups_dt <- merge(geoLungMeta_dt, lungFreqGroups_dt, by = "Sample", sort = F)

### Factorize
lungFreqGroups_dt$Title <- factor(lungFreqGroups_dt$Title, levels = lungTitleLevels_v)

### Make numeric
for (col in plotInfo_lsv$freqGroups) set(lungFreqGroups_dt, j = col, value = as.numeric(lungFreqGroups_dt[[col]]))

### Set order
setorder(lungFreqGroups_dt, Title)

### Add to list
freqGroupOut_lsdt[["Lung"]] <- lungFreqGroups_dt

###
### Blood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Merge with meta
bloodFreqGroups_dt <- merge(geoBloodMeta_dt, bloodFreqGroups_dt, by = "Sample", sort = F)

### Factorize
bloodFreqGroups_dt$Title <- factor(bloodFreqGroups_dt$Title, levels = bloodTitleLevels_v)

### Make numeric
for (col in plotInfo_lsv$freqGroups) set(bloodFreqGroups_dt, j = col, value = as.numeric(bloodFreqGroups_dt[[col]]))

### Set order
setorder(bloodFreqGroups_dt, Title)

### Add to list
freqGroupOut_lsdt[["Blood"]] <- bloodFreqGroups_dt

###
### Top Clone Freq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

topClones_lsdt <- list()

###
### Tumor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove file extension
tumorTopClones_dt$Sample <- gsub("_[a-z]+.*$|\\..*$", "", tumorTopClones_dt$Sample)

### Merge
tumorTopClones_dt <- merge(geoTumorMeta_dt, tumorTopClones_dt, by = "Sample", sort = F)

### Factorize
tumorTopClones_dt$Title <- factor(tumorTopClones_dt$Title, levels = tumorTitleLevels_v)

### Set order
setorder(tumorTopClones_dt, Title)

### Add to list
topClones_lsdt[["Tumor"]] <- tumorTopClones_dt

###
### Lung ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove file extension
lungTopClones_dt$Sample <- gsub("_[a-z]+.*$|\\..*$", "", lungTopClones_dt$Sample)

### Merge
lungTopClones_dt <- merge(geoLungMeta_dt, lungTopClones_dt, by = "Sample", sort = F)

### Factorize
lungTopClones_dt$Title <- factor(lungTopClones_dt$Title, levels = lungTitleLevels_v)

### Set order
setorder(lungTopClones_dt, Title)

### Add to list
topClones_lsdt[["Lung"]] <- lungTopClones_dt

###
### Blood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove file extension
bloodTopClones_dt$Sample <- gsub("_[a-z]+.*$|\\..*$", "", bloodTopClones_dt$Sample)

### Merge
bloodTopClones_dt <- merge(geoBloodMeta_dt, bloodTopClones_dt, by = "Sample", sort = F)

### Factorize
bloodTopClones_dt$Title <- factor(bloodTopClones_dt$Title, levels = bloodTitleLevels_v)

### Set order
setorder(bloodTopClones_dt, Title)

### Add to list
topClones_lsdt[["Blood"]] <- bloodTopClones_dt

###
### Jaccard ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

jaccard_lsdt <- list()

###
### Tumor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove extra file stuff
names(tumorClones_lsdt) <- gsub("_[a-z].*$", "", names(tumorClones_lsdt))

### Clean sample column, change column names
tumorClones_lsdt <- sapply(names(tumorClones_lsdt), function(x) {
  y <- tumorClones_lsdt[[x]]
  y$sample <- x
  colnames(y)[colnames(y) == "normC"] <- "Clones"
  colnames(y)[colnames(y) == "normFreq"] <- "Proportion"
  colnames(y)[colnames(y) == "Seq"] <- "CDR3.nt"
  colnames(y)[colnames(y) == "aaSeq"] <- "CDR3.aa"
  colnames(y)[colnames(y) == "V"] <- "V.name"
  colnames(y)[colnames(y) == "J"] <- "J.name"
  return(y)}, simplify = F, USE.NAMES = T)

### Subset meta
geoTumorSubMeta_dt <- geoTumorMeta_dt[Experiment == "d105" &
                                        Treatment %in% c("PLX+PTX+aPD-1_R", "PLX+PTX+aPD-1_NR", "PLX+PTX+aPD-L1"),]

### Subset data
subTumorClones_lsdt <- tumorClones_lsdt[geoTumorSubMeta_dt$Sample]

### Calculate overlap
tumorOverlap <- repOverlap(subTumorClones_lsdt, .method = "jaccard", .col = "v+aa+j")

### Melt and rename
tumorOverlapMelt_dt <- melt(convertDFT(as.data.frame(tumorOverlap)), id.vars = "V1")
colnames(tumorOverlapMelt_dt) <- c("S1", "S2", "Jaccard")

### Remove NA's (self-comparisons)
tumorOverlapMelt_dt <- tumorOverlapMelt_dt[!is.na(Jaccard)]

### Add treatment
tumorOverlapMelt_dt <- merge(tumorOverlapMelt_dt, geoTumorSubMeta_dt[,mget(c("Sample", "Treatment"))],
                             by.x = "S1", by.y = "Sample", sort = F)
colnames(tumorOverlapMelt_dt)[colnames(tumorOverlapMelt_dt) == "Treatment"] <- "T1"
tumorOverlapMelt_dt <- merge(tumorOverlapMelt_dt, geoTumorSubMeta_dt[,mget(c("Sample", "Treatment"))],
                             by.x = "S2", by.y = "Sample", sort = F)
colnames(tumorOverlapMelt_dt)[colnames(tumorOverlapMelt_dt) == "Treatment"] <- "T2"

### Combine treatments
tumorOverlapMelt_dt[,treatmentComparison := paste(T1, T2, sep = "_")]

### Subset
tumorOverlapMelt_dt <- tumorOverlapMelt_dt[treatmentComparison %in% treatCompLevels_v]

### Change Sample names
tumorOverlapMelt_dt <- merge(geoTumorSubMeta_dt[,mget(c("geoSample", "Sample"))], tumorOverlapMelt_dt,
                             by.x = "Sample", by.y = "S1")
tumorOverlapMelt_dt$Sample <- NULL
colnames(tumorOverlapMelt_dt)[colnames(tumorOverlapMelt_dt) == "geoSample"] <- "S1"

tumorOverlapMelt_dt <- merge(geoTumorSubMeta_dt[,mget(c("geoSample", "Sample"))], tumorOverlapMelt_dt,
                             by.x = "Sample", by.y = "S2")
tumorOverlapMelt_dt$Sample <- NULL
colnames(tumorOverlapMelt_dt)[colnames(tumorOverlapMelt_dt) == "geoSample"] <- "S2"

### Factorize and sort
tumorOverlapMelt_dt$treatmentComparison <- factor(tumorOverlapMelt_dt$treatmentComparison, levels = treatCompLevels_v)
setorder(tumorOverlapMelt_dt, treatmentComparison)

### Add to list
jaccard_lsdt[["Tumor"]] <- tumorOverlapMelt_dt

###
### Tumor Top 100 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset data
subTumorClonesTop100_lsdt <- sapply(subTumorClones_lsdt, function(x) {
  setorder(x, cols = -Proportion)
  x <- x[1:100,]
  return(x)}, simplify = F, USE.NAMES = T)

### Calculate overlap
tumorOverlapTop100 <- repOverlap(subTumorClonesTop100_lsdt, .method = "jaccard", .col = "v+aa+j")

### Melt and rename
tumorOverlapTop100Melt_dt <- melt(convertDFT(as.data.frame(tumorOverlapTop100)), id.vars = "V1")
colnames(tumorOverlapTop100Melt_dt) <- c("S1", "S2", "Jaccard")

### Remove NA's (self-comparisons)
tumorOverlapTop100Melt_dt <- tumorOverlapTop100Melt_dt[!is.na(Jaccard)]

### Add treatment
tumorOverlapTop100Melt_dt <- merge(tumorOverlapTop100Melt_dt, geoTumorSubMeta_dt[,mget(c("Sample", "Treatment"))],
                                   by.x = "S1", by.y = "Sample", sort = F)
colnames(tumorOverlapTop100Melt_dt)[colnames(tumorOverlapTop100Melt_dt) == "Treatment"] <- "T1"
tumorOverlapTop100Melt_dt <- merge(tumorOverlapTop100Melt_dt, geoTumorSubMeta_dt[,mget(c("Sample", "Treatment"))],
                                   by.x = "S2", by.y = "Sample", sort = F)
colnames(tumorOverlapTop100Melt_dt)[colnames(tumorOverlapTop100Melt_dt) == "Treatment"] <- "T2"

### Combine treatments
tumorOverlapTop100Melt_dt[,treatmentComparison := paste(T1, T2, sep = "_")]

### Subset
tumorOverlapTop100Melt_dt <- tumorOverlapTop100Melt_dt[treatmentComparison %in% treatCompLevels_v]

### Change Sample names
tumorOverlapTop100Melt_dt <- merge(geoTumorSubMeta_dt[,mget(c("geoSample", "Sample"))], tumorOverlapTop100Melt_dt,
                                   by.x = "Sample", by.y = "S1")
tumorOverlapTop100Melt_dt$Sample <- NULL
colnames(tumorOverlapTop100Melt_dt)[colnames(tumorOverlapTop100Melt_dt) == "geoSample"] <- "S1"

tumorOverlapTop100Melt_dt <- merge(geoTumorSubMeta_dt[,mget(c("geoSample", "Sample"))], tumorOverlapTop100Melt_dt,
                                   by.x = "Sample", by.y = "S2")
tumorOverlapTop100Melt_dt$Sample <- NULL
colnames(tumorOverlapTop100Melt_dt)[colnames(tumorOverlapTop100Melt_dt) == "geoSample"] <- "S2"

### Factorize and sort
tumorOverlapTop100Melt_dt$treatmentComparison <- factor(tumorOverlapTop100Melt_dt$treatmentComparison, levels = treatCompLevels_v)
setorder(tumorOverlapTop100Melt_dt, treatmentComparison)

### Add to list
jaccard_lsdt[["TumorTop100"]] <- tumorOverlapTop100Melt_dt

###
### Lung ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove extra file stuff
names(lungClones_lsdt) <- gsub("_[a-z].*$", "", names(lungClones_lsdt))

### Clean sample column, change column names
lungClones_lsdt <- sapply(names(lungClones_lsdt), function(x) {
  y <- lungClones_lsdt[[x]]
  y$sample <- x
  colnames(y)[colnames(y) == "normC"] <- "Clones"
  colnames(y)[colnames(y) == "normFreq"] <- "Proportion"
  colnames(y)[colnames(y) == "Seq"] <- "CDR3.nt"
  colnames(y)[colnames(y) == "aaSeq"] <- "CDR3.aa"
  colnames(y)[colnames(y) == "V"] <- "V.name"
  colnames(y)[colnames(y) == "J"] <- "J.name"
  return(y)}, simplify = F, USE.NAMES = T)

### Subset meta
geoLungSubMeta_dt <- geoLungMeta_dt[Experiment == "d105" &
                                      Treatment %in% c("PLX+PTX+aPD-1_R", "PLX+PTX+aPD-1_NR", "PLX+PTX+aPD-L1"),]

### Subset data
subLungClones_lsdt <- lungClones_lsdt[geoLungSubMeta_dt$Sample]

### Calculate overlap
lungOverlap <- repOverlap(subLungClones_lsdt, .method = "jaccard", .col = "v+aa+j")

### Melt and rename
lungOverlapMelt_dt <- melt(convertDFT(as.data.frame(lungOverlap)), id.vars = "V1")
colnames(lungOverlapMelt_dt) <- c("S1", "S2", "Jaccard")

### Remove NA's (self-comparisons)
lungOverlapMelt_dt <- lungOverlapMelt_dt[!is.na(Jaccard)]

### Add treatment
lungOverlapMelt_dt <- merge(lungOverlapMelt_dt, geoLungSubMeta_dt[,mget(c("Sample", "Treatment"))],
                            by.x = "S1", by.y = "Sample", sort = F)
colnames(lungOverlapMelt_dt)[colnames(lungOverlapMelt_dt) == "Treatment"] <- "T1"
lungOverlapMelt_dt <- merge(lungOverlapMelt_dt, geoLungSubMeta_dt[,mget(c("Sample", "Treatment"))],
                            by.x = "S2", by.y = "Sample", sort = F)
colnames(lungOverlapMelt_dt)[colnames(lungOverlapMelt_dt) == "Treatment"] <- "T2"

### Combine treatments
lungOverlapMelt_dt[,treatmentComparison := paste(T1, T2, sep = "_")]

### Subset
lungOverlapMelt_dt <- lungOverlapMelt_dt[treatmentComparison %in% treatCompLevels_v]

### Change Sample names
lungOverlapMelt_dt <- merge(geoLungSubMeta_dt[,mget(c("geoSample", "Sample"))], lungOverlapMelt_dt,
                            by.x = "Sample", by.y = "S1")
lungOverlapMelt_dt$Sample <- NULL
colnames(lungOverlapMelt_dt)[colnames(lungOverlapMelt_dt) == "geoSample"] <- "S1"

lungOverlapMelt_dt <- merge(geoLungSubMeta_dt[,mget(c("geoSample", "Sample"))], lungOverlapMelt_dt,
                            by.x = "Sample", by.y = "S2")
lungOverlapMelt_dt$Sample <- NULL
colnames(lungOverlapMelt_dt)[colnames(lungOverlapMelt_dt) == "geoSample"] <- "S2"

### Factorize and sort
lungOverlapMelt_dt$treatmentComparison <- factor(lungOverlapMelt_dt$treatmentComparison, levels = treatCompLevels_v)
setorder(lungOverlapMelt_dt, treatmentComparison)

### Add to list
jaccard_lsdt[["Lung"]] <- lungOverlapMelt_dt

###
### Clonotype Tracking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

trackCloneQuery_lsdt <- list()
trackCloneCorpus_lsdt <- list()

###
### Tumor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset data (this is for the query)
subTumorClonesTop50_lsdt <- sapply(subTumorClones_lsdt, function(x) {
  setorder(x, cols = -Proportion)
  if (nrow(x) >= 50) x <- x[1:50,]
  return(x)}, simplify = F, USE.NAMES = T)

### Shared 3x aPD-1 R (this is for the query)
tumor3xRSamples_v <- geoTumorSubMeta_dt[Treatment == "PLX+PTX+aPD-1_R", Sample]
tumorTop50_3xR_dt <- do.call(rbind, subTumorClonesTop50_lsdt[names(subTumorClonesTop50_lsdt) %in% tumor3xRSamples_v])
tumorSharedTop50_3xR_dt <- getShared(tumorTop50_3xR_dt)

### Shared 3x aPD-L1 (this is for the query)
tumor3xaPDL1Samples_v <- geoTumorSubMeta_dt[Treatment == "PLX+PTX+aPD-L1", Sample]
tumorTop50_3xaPDL1_dt <- do.call(rbind, subTumorClonesTop50_lsdt[names(subTumorClonesTop50_lsdt) %in% tumor3xaPDL1Samples_v])
tumorSharedTop50_3xaPDL1_dt <- getShared(tumorTop50_3xaPDL1_dt)

### Add to output
trackCloneQuery_lsdt[["Tumor3xR"]] <- tumorSharedTop50_3xR_dt
trackCloneQuery_lsdt[["Tumor3xaPDL1"]] <- tumorSharedTop50_3xaPDL1_dt
trackCloneCorpus_lsdt[["TumorTop50"]] <- subTumorClonesTop50_lsdt
trackCloneCorpus_lsdt[["TumorAll"]] <- subTumorClones_lsdt

###
### Lung ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Subset data
subLungClonesTop50_lsdt <- sapply(subLungClones_lsdt, function(x) {
  setorder(x, cols = -Proportion)
  if (nrow(x) >= 50) x <- x[1:50,]
  return(x)}, simplify = F, USE.NAMES = T)

### Shared 3x aPD-1 R (this is for the query)
lung3xRSamples_v <- geoLungSubMeta_dt[Treatment == "PLX+PTX+aPD-1_R", Sample]
lungTop50_3xR_dt <- do.call(rbind, subLungClonesTop50_lsdt[names(subLungClonesTop50_lsdt) %in% lung3xRSamples_v])
lungSharedTop50_3xR_dt <- getShared(lungTop50_3xR_dt)

### Shared 3x aPD-L1 (this is for the query)
lung3xaPDL1Samples_v <- geoLungSubMeta_dt[Treatment == "PLX+PTX+aPD-L1", Sample]
lungTop50_3xaPDL1_dt <- do.call(rbind, subLungClonesTop50_lsdt[names(subLungClonesTop50_lsdt) %in% lung3xaPDL1Samples_v])
lungSharedTop50_3xaPDL1_dt <- getShared(lungTop50_3xaPDL1_dt)

### Add to output
trackCloneQuery_lsdt[["Lung3xR"]] <- lungSharedTop50_3xR_dt
trackCloneQuery_lsdt[["Lung3xaPDL1"]] <- lungSharedTop50_3xaPDL1_dt
trackCloneCorpus_lsdt[["LungTop50"]] <- subLungClonesTop50_lsdt

###
### Blood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Remove extra file stuff
names(bloodClones_lsdt) <- gsub("_[a-z].*$", "", names(bloodClones_lsdt))

### Clean sample column, change column names
bloodClones_lsdt <- sapply(names(bloodClones_lsdt), function(x) {
  y <- bloodClones_lsdt[[x]]
  y$sample <- x
  colnames(y)[colnames(y) == "nb.clone.count"] <- "Clones"
  colnames(y)[colnames(y) == "nb.clone.fraction"] <- "Proportion"
  colnames(y)[colnames(y) == "nSeqCDR3"] <- "CDR3.nt"
  colnames(y)[colnames(y) == "aaSeqCDR3"] <- "CDR3.aa"
  colnames(y)[colnames(y) == "V segments"] <- "V.name"
  colnames(y)[colnames(y) == "J segments"] <- "J.name"
  return(y)}, simplify = F, USE.NAMES = T)

### Subset meta
geoBloodSubMeta_dt <- geoBloodMeta_dt[Experiment == "d105" &
                                        Treatment %in% c("PLX+PTX+aPD-1_R", "PLX+PTX+aPD-1_NR", "PLX+PTX+aPD-L1"),]

### Subset data (this is for the query)
subBloodClones_lsdt <- bloodClones_lsdt[geoBloodSubMeta_dt$Sample]

### Subset data for top 50 (this is for the query)
subBloodClonesTop50_lsdt <- sapply(subBloodClones_lsdt, function(x) {
  setorder(x, cols = -Proportion)
  if (nrow(x) >= 50) x <- x[1:50,]
  return(x)}, simplify = F, USE.NAMES = T)

### Add to output
trackCloneCorpus_lsdt[["BloodAll"]] <- subBloodClones_lsdt
trackCloneCorpus_lsdt[["BloodTop50"]] <- subBloodClonesTop50_lsdt

###
### Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

saveRDS(divOut_lsdt, file = file.path(baseDir_v, "data/fig7/rds/diversityMetrics_lsdt.rds"))
saveRDS(freqGroupOut_lsdt, file = file.path(baseDir_v, "data/fig7/rds/freqGroups_lsdt.rds"))
saveRDS(topClones_lsdt, file = file.path(baseDir_v, "data/fig7/rds/topClones_lsdt.rds"))
saveRDS(jaccard_lsdt, file = file.path(baseDir_v, "data/fig7/rds/jaccard_lsdt.rds"))
saveRDS(trackCloneQuery_lsdt, file = file.path(baseDir_v, "data/fig7/rds/trackCloneQuery_lsdt.rds"))
saveRDS(trackCloneCorpus_lsdt, file = file.path(baseDir_v, "data/fig7/rds/trackCloneCorpus_lsdt.rds"))

