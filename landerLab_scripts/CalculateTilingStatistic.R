##################################################################
## Jesse Engreitz
## April 28, 2016
##
## Example command:
## use .r-3.1.1-bioconductor-3.0; R_LIBS_SITE=; reuse Python-2.7
## Rscript CalculateTilingStatistic.R --input /seq/lincRNA/cfulco/MYCTiling/160303.EssentialityTiling/160328.HiSeq.CRISPRiScreensMoreMaterial/160329.ScreenSummaries/160501.CiTile.database.GATA1test.bed --output /seq/lincRNA/RAP/Paper/MYCTiling/Figures/WindowSize/test.bed --scoreColumn log2FC --T0columns RAW_CiTile.1.T0,RAW_CiTile.8.T0

## Russell Ryan mods 2019 - Option of normalizing log2FoldChange values of filtered sgRNAs to median or mean of negative controls

suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option(c("-i", "--input"), type="character", help="Input database file containing at least the following columns: chr, start, end, GuideSequence, an arbitrarily named score column, and a 'target' column that includes negative_control labeled guides"),
  make_option(c("-o", "--output"), type="character", help="Output prefix for tab-delimited text file"),
  make_option(c("--T0columns"), type="character", help="Comma-separated list of columns to apply the minT0count filter"),
  make_option(c("-c", "--scoreColumn"), type="character", default="log2FC", help="Column name calculating scores and statistics"),
  make_option(c("-w", "--grnaWindowSize"), type="numeric", default=20, help="Number of guides to tile in windows"),
  make_option(c("-m", "--maxSpan"), type="numeric", default=500, help="Max span of the guides in a given window"),
  make_option(c("-s", "--minOffTargetScore"), type="numeric", default=50, help="Minimum off-target score"),
  make_option(c("-t", "--minT0count"), type="numeric", default=50, help="Minimum T0 counts"),
  make_option(c("--minGuideLength"), type="numeric", default=20, help="Minimum length of guideRNAs to use"),
  make_option(c("--maxGuideLength"), type="numeric", default=21, help="Maximum length of guideRNAs to use"),
  make_option(c("--maxGsInGuideSequence"), type="numeric", default=10, help="Maximum number of G's allows in guide sequence"),
  make_option(c("-n", "--normVals"), type="character", default="median", help="Control sgRNA normalization: 'median', 'mean', or 'none'"), ## RR added for normalization
  make_option("--codeDir", type="character", default="/seq/lincRNA/RAP/Promoters/LanderLab-EP-Prediction/", help="Absolute path to code repository")
)
opt <- parse_args(OptionParser(option_list=option.list))

suppressPackageStartupMessages(source(paste0(opt$codeDir, "/src/libs/JuicerUtilities.R")))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "/src/CRISPRScreen/JuicerCRISPRi.R")))



## PARAMETERS
t0columns <- strsplit(gsub("-",".",opt$T0columns), ",")[[1]]
window.size <- opt$grnaWindowSize
max.span <- opt$maxSpan
data.col <- gsub("-",".",opt$scoreColumn)
normVals <- opt$normVals #RR added for normalization

x <- read.delim(gzfile(opt$input))
stopifnot(all(c("chr","start","end",data.col,"target","GuideSequence") %in% colnames(x)))

x$nG <- sapply(as.character(as.matrix(x$GuideSequence)), function(x) sum(strsplit(x, "")[[1]] == "G"))
x$GuideLength <- sapply(as.character(as.matrix(x$GuideSequence)), nchar)

to.use <- subset(x, OffTargetScore >= opt$minOffTargetScore & 
                    nG <= opt$maxGsInGuideSequence & 
                    GuideLength >= opt$minGuideLength & 
                    GuideLength <= opt$maxGuideLength)

#TWS Fix: must remove NA scores in order to do calculations.
to.use <- to.use[!is.na(to.use[,data.col]),]

for (col.name in t0columns) {
  to.use <- subset(to.use, get(col.name) >= opt$minT0count)
}
cat(nrow(to.use),"of",nrow(x),"total guideRNAs passed the filters\n")
stopifnot(nrow(to.use) > 0)

control.vals <- subset(to.use, target == "negative_control")[,data.col]
cat("Using",length(control.vals),"as negative control guides\n")
stopifnot(length(control.vals) > 0)

## RR added normlization of log2FoldChange values
if (normVals == "none") {
  cat("No normalization of log2FoldChange to negative control values")
  to.use$mean_log2FC_norm <- to.use[,data.col]
  x <- x[,data.col]
  
} else if (normVals == "mean") {
  cat("Normalizing log2FoldChange to mean of negative control values")
  to.use$mean_log2FC_norm <- to.use[,data.col] - mean(subset(to.use, target == "negative_control")[,data.col])
  x$mean_log2FC_norm <- x[,data.col] - mean(subset(to.use, target == "negative_control")[,data.col])
  
} else if (normVals == "median") {
  cat("Normalizing log2FoldChange to median of negative control values")
  to.use$mean_log2FC_norm <- to.use[,data.col] - median(subset(to.use, target == "negative_control")[,data.col])
  x$mean_log2FC_norm <- x[,data.col] - median(subset(to.use, target == "negative_control")[,data.col])
  
} else {
  cat("Normalizing log2FoldChange to median of negative control values by default")
  to.use$mean_log2FC_norm <- to.use[,data.col] - median(subset(to.use, target == "negative_control")[,data.col])
  x$mean_log2FC_norm <- x[,data.col] - median(subset(to.use, target == "negative_control")[,data.col])
}

#data.col.norm <- "mean_log2FC_norm"

write.table(x, file=paste0(sub("txt","",opt$input),"norm.txt"), sep='\t', quote=F, row.names=F, col.names=T, na = "")
write.table(to.use, file=paste0(sub("txt","",opt$input),"filt.norm.txt"), sep='\t', quote=F, row.names=F, col.names=T, na = "")
cat("Before calculateSlidingWindowStatistic")
result <- calculateSlidingWindowStatistic(subset(to.use, target != "negative_control"), control.vals, window.size, max.span, data.col=data.col) #data.col=data.col.norm #values are normalized within function
cat("After calculateSlidingWindowStatistic")
write.table(result, file=paste0(opt$output, "_tiling_stat.txt"), sep='\t', quote=F, row.names=F, col.names=T)
cat("Table written")

quit()
#TWS - No need for these additional files to be created. Not valid bedgraph in any case.
#to.write.mean <- with(result, data.frame(chr=chr, start=start, end=end, score=-mean))
#write.table(to.write.mean, file=paste0(opt$output,".mean.bedgraph"), sep='\t', quote=F, row.names=F, col.names=F)

#to.write.u <- with(result, data.frame(chr=chr,start=start,end=end,score=log10(fdr.utest)*sign(mean)))
#write.table(to.write.u, file=paste0(opt$output,".UFDR.bedgraph"), sep='\t', quote=F, row.names=F, col.names=F)

#to.write.gRNA <- with(subset(to.use,!is.na(start)), data.frame(chr=chr, start=start, end=end, score=-get(data.col)))
#write.table(to.write.gRNA, file=paste0(opt$output,".IndividualGuides.bedgraph", sep='\t', quote=F, row.names=F, col.names=F)


