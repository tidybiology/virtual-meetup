#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE29797", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL11180", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] 
# "GSExxxxx" = "Series accession number". If the Series is associated with multiple Platforms (e.g. different chip types),
# you will be asked to select the Platform of interest
# Also, this changes the class from "list" to 
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase"
# and typeof() changes from "list" to "S4" object

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset)) # This step slightly changes the column names to match the toptable format

# group names for all samples
gsms <- "X0X0X0X1X1X1XXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) # This is a way of checking if the data have ALREADY been log transformed.
# For example, if `qx[5]` (the 99th quantile) really is greater than 100, then the anti-log value will be more than 2^100,
# which is crazy! If this `qx[5]` really is greater than 100, then the data have likely not yet been log-transformed
# This is GEO2R's "auto-detect feature" for determining if values have been log-transformed
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# Boxplot for selected GEO samples
palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE29797", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2)
# You do this to assess whether the distribution of values across Samples are median-centered, which is generally 
# indicative that the data are normalized and cross-comparable

# Finally, while most data have been pre-processed, if you want confirmation of this, check the "data processing"
# field - 
gset@phenoData@data[["data_processing"]]
# "The data were analyzed with Partek Genomics Suite 6.5 using the RMA algorithmn" 
# So, in this case, you've confirmed that background correction and quantile normalization have been performed!

# If all this looks good, proceed!

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl # before this step, `gset$description` has names like `V3`, `V5`, and so on, which then change to 
# `G0`, `G0`, and so on after you run `gset$description <- fl`
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01) # 0.01 refers to the "assumed proportion of genes which are differentially expressed"
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# alternatively (preferably? :))
design_alt <- model.matrix(~ description, gset)
fit_alt <- lmFit(gset, design_alt)
cont.matrix.alt <- makeContrasts(descriptionG1, levels=design_alt)
fit2_alt <- contrasts.fit(fit_alt, cont.matrix.alt)
fit2_alt <- eBayes(fit2_alt, 0.01)
tT_alt <- topTable(fit2_alt, adjust="fdr", sort.by="B", number=250)

identical(tT$ID, tT_alt$ID) # They're the same!

# final point about `ex` and `gset``- they're effectively interchangeable! I find it easier to use `ex` to design more
# complex design matrices
# replace `gset` with `ex`
design_alt <- model.matrix(~ description, gset)
fit_ex <- lmFit(ex, design_alt)
cont.matrix.ex <- makeContrasts(descriptionG1, levels=design_alt)
fit2_ex <- contrasts.fit(fit_ex, cont.matrix.ex)
fit2_ex <- eBayes(fit2_ex, 0.01)
tT_ex <- topTable(fit2_ex, adjust="fdr", sort.by="B", number=250)
# While the logFC, adjusted p-values and so on are the same, notice that `tT_ex` only has 6 columns, as opposed to 27, like
# the other two

# To save the table
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Simple!
# Some questions -
# Do we assume that the data have already been quantile normalized because that is one of the pre-processing steps?
# From GEO - for the most part, yes, but USERS SHOULD CONFIRM THIS!
