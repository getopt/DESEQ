# Usage:
# R --vanilla --args --table=<file.tab>
#                    --header=<YES|NO>
#                    --rownamesCol=<N>
#                    --fitType=<paramteric|local>
#                    --method=<pooled|blind>
#                    --sharingMode=<maximum|fit-only>
#                    --treatCountCols=<1-based column numbers, csv>
#                    --contrCountCols=<1-based column numbers, csv>
#                    --outfiBaseTAB=<path/nameBase>
#                    --outfiBasePNG=<path/nameBase> 
#                    --filterNoPolya

# Description:
#   --table           <   file with first column as row names and two columns
#                         providing counts (values will be rounded to integers)
#   --header          <   "YES" or "NO" (default "NO"; header isn't important
#                         as such, only important for correct reading of the table)
#   --rownamesCol     <   column number for the column that provides rownames


#  --treatCountCols   <   <1-based column numbers, csv>
#  --contrCountCols   <   <1-based column numbers, csv>


#   --method          <   "pooled" or "blind" method of fitting dispersion model. 
#                         Use "pooled" if you have replicates, "blind" if you 
#                         do not have replicates. Default: "pooled"

#   --sharingMode     <   "maximum" or "fit-only" mode of fitting dispersion 
#                         model. Use "maximum" if you have replicates, 
#                         "fit-only" if you do not have replicates. Default: "maximum"

#   --fitType         <   "parametric" or "local" mode of fitting dispersion model.
#                         "parametric" (default) is recommended. But it doesn't 
#                          work sometimes for unknown reasons, suggested alternative 
#                         "local". Default: "parametric"


#   --outfiBaseTAB    <   path/nameBase of for the output files: 
#                                             ${base}.deseq.Results_all.tab
#                                             ${base}.deseq.Results_all.tab.SKIPPED_IDS 
#   --outfiBasePNG    <   path/nameBase of for output PDFs. Two will be written: 
#                                                        ${base}.deseq.dispersion.png 
#                                                        ${base}.deseq.diffexpres.png 


#   --filterNoPolya   <  set to "YES" if to remove certain gene names from the
#                        table (those that are matched by "toFilter") by reg exp 
#                        matching to gene names: 
#                                                   "tRNA" "rRNA" "\\)n"
#                                                   "4.5S" "5S" "7S" "_rich",
#                                                   "\\|U\\d+" "snoRNA" "mir-"

# NOTE: more filtering commented out inside: 
#   rows that have 0-value in all columns are not included into the analysis

library(R.utils)
library(DESeq)

args <- R.utils::commandArgs(asValues=TRUE)

table.file       <-    NULL

treatCountCols   <-    NULL
contrCountCols   <-    NULL
filterNoPolya    <-    "NO"
skipped_ids      <-    NULL

rownamesCol      <-    NULL 
outfiBaseTAB     <-    "outfi"
outfiBasePNG     <-    "outfi"
padjCutoff       <-    0.2
pvalCutoff       <-    0.05
headerPresent    <-    "NO"

### Params for 'estimateDispersions'
# with     replicates use "pooled" "maximum"   "parametric"
# without  replicates use "blind"  "fit-only"  "parametric"|"local" 
method          <-    "pooled"
sharingMode     <-    "maximum"
fitType         <-    "parametric"

if(!is.null(args[["table"]]))
    table.file        <-   args$table
if(!is.null(args[["treatCountCols"]]))
    treatCountCols    <-   as.numeric(unlist(strsplit(args$treatCountCols, ",")))
if(!is.null(args[["contrCountCols"]]))
    contrCountCols    <-   as.numeric(unlist(strsplit(args$contrCountCols, ",")))
if(!is.null(args[["rownamesCol"]]))
    rownamesCol       <-   as.numeric(args$rownamesCol)
if(!is.null(args[["outfiBaseTAB"]]))
    outfiBaseTAB      <-   args$outfiBaseTAB
if(!is.null(args[["outfiBasePNG"]]))
    outfiBasePNG      <-   args$outfiBasePNG
if(!is.null(args[["header"]]))
    headerPresent     <-   args$header
if(!is.null(args[["fitType"]]))
    fitType           <-   args$fitType
if(!is.null(args[["method"]]))
    method            <-   args$method
if(!is.null(args[["sharingMode"]]))
    sharingMode       <-   args$sharingMode
if(!is.null(args[["filterNoPolya"]]))
    filterNoPolya   <-   args$filterNoPolya


if (is.null(table.file) | is.null(treatCountCols) |  is.null(contrCountCols))
    stop("must have --table=<table.file> and 
                    --treatCountCols=<1-based column numbers, csv> options and
                    --contrCountCols=<1-based column numbers, csv> options" )

if(headerPresent == "NO"){
    table           <-   read.table(table.file, sep = "\t", header = FALSE, as.is = T)
    rownames(table) <-   table[,rownamesCol]
}
if(headerPresent == "YES"){
    table         <-   read.table(table.file, sep = "\t", header = TRUE, as.is = T, skip = 1)
    rownames(table) <-   table[,rownamesCol]
}
if(headerPresent != "YES" & headerPresent != "YES"){
    stop("cannot interpret --header <YES|NO> option!")
}

if(filterNoPolya == "YES"){
    toFilter      <-   c("tRNA","rRNA","\\)n","4.5S","5S","7S","_rich","\\|U\\d+","snoRNA","mir-")
    indToFilter   <-   unlist(lapply(toFilter, function(x){grep(x,rownames(table), perl = T)}))
    skipped_ids   <-   c(skipped_ids, rownames(table[indToFilter,]))
    if(length(indToFilter) > 0){
        table         <-   table[-(indToFilter),]
    }
}

treatCountTable    <-   round(table[,treatCountCols, drop = F])
contrCountTable    <-   round(table[,contrCountCols, drop = F])
colnames(treatCountTable) <- paste("treatment_", c(1:length(treatCountTable[1,])), sep = "")
colnames(contrCountTable) <- paste("control_", c(1:length(contrCountTable[1,])), sep = "")

countTable <- cbind(treatCountTable, contrCountTable)

#
#  f i l t e r i n g 
# 

# # remove raws if   A N Y   of the two samples have 0-count
# zeroCheck <- apply(countTable, 1, 
#                    function(x){
#                        zeroCheckCurrent <- FALSE
#                        for(i in c(1:length(x))){
#                            if(x[i] == 0){
#                                zeroCheckCurrent <- TRUE
#                            }
#                        }
#                        return(zeroCheckCurrent)
#                    })
# skipped_ids <- c(skipped_ids, rownames(countTable[zeroCheck,]))
# countTable  <- countTable[!zeroCheck,]            
 
# # remove raws if   B O T H   of the two samples have 0-count
# # WARNING: you get infinities
# zeroCheck <- apply(countTable, 1, 
#                    function(x){
#                        zeroCheckCurrent <- FALSE
#                        zeroCheck1 <- FALSE
#                        zeroCheck2 <- FALSE
#                        for(i in c(1:length(x))){
#                            if(x[i] == 0 & i == 1){
#                                zeroCheck1 <- TRUE
#                            }
#                            if(x[i] == 0 & i == 2){
#                                zeroCheck2 <- TRUE
#                            }
#                        }
#                        if(zeroCheck1 == TRUE & zeroCheck2 == TRUE){
#                            zeroCheckCurrent <- TRUE
#                        }
#                        return(zeroCheckCurrent)
#                    })
# skipped_ids <- c(skipped_ids, rownames(countTable[zeroCheck,]))
# countTable <- countTable[!zeroCheck,]            

# 
# MAKE COUNT DATASET
# 
designTable     <-    data.frame(row.names = colnames(countTable), 
                                 condition = c(
                                   rep("treated",    times = length(treatCountTable[1,])),
                                   rep("untreated",  times = length(contrCountTable[1,]))
                                               ))

conds           <-    designTable$condition
cds             <-    newCountDataSet(countTable, conds)

# #
# # m o r e  f i l t e r i n g 
# # 
# exclude lower counts 
# rs  <- rowSums(counts(cds))
# use <- (rs > quantile(rs, .75))
# cds <- cds[use,]


#
# S I Z E   F A C T O R S  &  D I S P E R S I O N
#
cds   <- estimateSizeFactors(cds)
cds   <- estimateDispersions(cds, method = method, sharingMode= sharingMode , fitType = fitType)
# cds   <- estimateDispersions(cds, method="per-condition", sharingMode="maximum" ,
#                                                   fitType = fitType) ## this has to be tested

outfile  <-  paste(outfiBasePNG, ".deseq.dispersion.png", sep = "")
mainT    <-  sub("..*\\/","", outfile)
 png(outfile, height = 800, width = 800)
    plot(rowMeans(counts(cds, normalized=TRUE)), fitInfo(cds)$perGeneDispEsts,
                                             pch = '.', log="xy", main = mainT)
    xg <- 10^seq( -.5, 5, length.out=300)
    lines(xg, fitInfo(cds)$dispFun(xg), col="red")
 dev.off()


#
# N B I N O M  T E S T
#
res    <- nbinomTest(cds, "untreated", "treated") # THIS IS TESTED CORRECT DIRECTIONALITY 
                                                  # WHEN DATA ORIGINALLY
                                                  # SUPPLIED IN
                                                  # "TREATMENT,CONTROL"
                                                  # DIRECTION THE RESULTING
                                                  # POSITIVE FOLD CHANGE MEANS
                                                  # UPREGULATION IN TREATMENT
                                                  # NEGATVIE
                                                  # DOWNREGULATION IN TREATMENT

outfile  <-  paste(outfiBasePNG, ".deseq.diffexpres.png", sep = "")
mainT    <-  sub("..*\\/","", outfile)
 png(outfile, height = 800, width =  800)
 plot(res$baseMean, res$log2FoldChange,
     log="x",pch=20, 
     col = ifelse(res$pval < pvalCutoff, "blue", "black"),
     cex = ifelse(res$pval < pvalCutoff,  .6,     .3),
     cex.main = 0.7, main = mainT)
 points(res[res$padj < padjCutoff, "baseMean"], res[res$padj < padjCutoff, "log2FoldChange"],
                                                              col = "red", pch = 19, cex = 1)
 dev.off()

#
# W R I T E   O U T P U T
#
res      <-  res[order(res$foldChange, decreasing = T),]
outfile  <-  paste(outfiBaseTAB, ".deseq.Results_all.tab", sep = "")
write.table(res, file = outfile, sep = "\t", quote = F, row.names = F)
write(skipped_ids, file = paste(outfile, ".SKIPPED_IDS", sep = ""), sep = "\n")

# #
# # subset for padjCutoff and write as a separate table
# #
# resSig   <-  res[res$padj < padjCutoff,]
# padjForName <- padjCutoff * 10
# outfile  <-  paste(outfiBaseTAB, ".deseq.Results_padj", padjForName, ".tab", sep = "")
# write.table(resSig, file = outfile, sep = "\t", quote = F, row.names = F)

