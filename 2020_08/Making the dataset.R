library(Biobase)
library(GEOquery)
library(limma)

gset <- getGEO("GSE29797", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL11180", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] 

fvarLabels(gset) <- make.names(fvarLabels(gset)) 
# group names for all samples
gsms <- "X0X0X0X1X1X1XXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

ex <- exprs(gset)
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
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(gset))

tT %>% 
  rownames_to_column("probe") -> tT
ex %>% 
  data.frame(.) %>% 
  rownames_to_column("probe") -> ex

ex %>% 
  inner_join(tT, by = "probe") %>% 
  select(Gene.symbol, contains("GSM")) -> ex

ex %>% 
  slice(1:250) -> ex

write_csv(ex, "df.csv")  
  
