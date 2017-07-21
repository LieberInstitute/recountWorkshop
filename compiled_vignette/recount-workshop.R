## ----Figure1, out.width="100%", fig.align="center", fig.cap = "Overview of the data available in recount2. Reads (pink boxes) aligned to the reference genome can be used to compute a base-pair coverage curve and identify exon-exon junctions (split reads). Gene and exon count matrices are generated using annotation information providing the gene (green boxes) and exon (blue boxes) coordinates together with the base-level coverage curve. The reads spanning exon-exon junctions (jx) are used to compute a third count matrix that might include un-annotated junctions (jx 3 and 4). Without using annotation information, expressed regions (orange box) can be determined from the base-level coverage curve to then construct data-driven count matrices.", echo = FALSE----
knitr::include_graphics("Figure1.png")

## ----Figure2, out.width="100%", fig.align="center", fig.cap = "recount2 provides coverage count matrices in RangedSummarizedExperiment (rse) objects. Once the rse object has been downloaded and loaded into R, the feature information is accessed with rowRanges(rse) (blue box), the counts with assays(rse)\\$counts (pink box) and the sample metadata with colData(rse) (green box). The sample metadata can be expanded using add\\_predictions(rse) (orange box) or with custom code (brown box) matching by a unique sample identifier such as the SRA Run id. The rse object is inside the purple box and matching data is highlighted in each box.", echo = FALSE----
knitr::include_graphics("Figure2.png")

## ----"install", eval = FALSE-----------------------------------------------
#  ## Install packages from Bioconductor
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(c("recount", "GenomicRanges", "DESeq2", "ideal", "regionReport",
#      "clusterProfiler", "org.Hs.eg.db", "gplots", "derfinder",
#      "rtracklayer", "GenomicFeatures", "bumphunter", "derfinderPlot",
#      "devtools"))

## ----"load libraries", message = FALSE, warning = FALSE--------------------
library("recount")
library("GenomicRanges")
library("DESeq2")
library("ideal")
library("regionReport")
library("clusterProfiler")
library("org.Hs.eg.db")
library("gplots")
library("derfinder")
library("rtracklayer")
library("GenomicFeatures")
library("bumphunter")
library("derfinderPlot")
library("devtools")

## ----Figure3, out.width="100%", fig.align="center", fig.cap = "RNA-seq starting data. 16 RNA-seq un-aligned RNA-seq reads 3 base-pairs long are shown (pink boxes) along a reference genome 16 base-pairs long (white box).", echo = FALSE----
knitr::include_graphics("Figure3.png")

## ----Figure4, out.width="100%", fig.align="center", fig.cap = "Aligned RNA-seq reads. Spice-aware RNA-seq aligners such as Rail-RNA are able to find the coordinates to which the reads map, even if they span exon-exon junctions (connected boxes). Rail-RNA soft clips some reads (purple boxes with rough edges) such that a portion of these reads align to the reference genome.", echo = FALSE----
knitr::include_graphics("Figure4.png")

## ----Figure5, out.width="100%", fig.align="center", fig.cap = "Gene annotation. A single gene with two isoforms composed by three distinct exons (blue boxes) is illustrated. Exons 1 and 3 share the first five base-pairs while exon 2 is common to both isoforms.", echo = FALSE----
knitr::include_graphics("Figure5.png")

## ----"disjoin", message = FALSE--------------------------------------------
library("GenomicRanges")
exons <- GRanges("seq", IRanges(start = c(1, 1, 13), end = c(5, 8, 15)))
exons
disjoin(exons)

## ----Figure6, out.width="100%", fig.align="center", fig.cap = "Disjoint exons. Windows of distinct exonic sequence for the example gene. Disjoint exons 1 and 2 form exon 1.", echo = FALSE----
knitr::include_graphics("Figure6.png")

## ----Figure7, out.width="100%", fig.align="center", fig.cap = "Base-pair coverage counting for exonic base-pairs. At each exonic base-pair we compute the number of reads overlapping that given base-pair. The first base (orange arrow) has 3 reads overlapping that base-pair. Base-pair 11 has a coverage of 3 but does not overlap known exonic sequence, so that information is not used for the gene and exon count matrices (grey arrow). If a read partially overlaps exonic sequence, only the portion that overlaps is used in the computation (see right most read).", echo = FALSE----
knitr::include_graphics("Figure7.png")

## ----Figure8, out.width="100%", fig.align="center", fig.cap = "Exon and gene coverage counts. The coverage counts for each disjoint exon are the sum of the base-pair coverage. The gene coverage count is the sum of the disjoint exons coverage counts.", echo = FALSE----
knitr::include_graphics("Figure8.png")

## ----"coverage"------------------------------------------------------------
## Take the example and translate it to R code
library("GenomicRanges")
reads <- GRanges("seq", IRanges(
    start = rep(
        c(1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 14), 
        c(3, 1, 2, 1, 2, 1, 2, 1, 2, 4, 1)
    ), width = rep(
        c(1, 3, 2, 3, 1, 2, 1, 3, 2, 3, 2, 1, 3),
        c(1, 4, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1)
    )
))
## Get the base-level genome coverage curve
cov <- as.integer(coverage(reads)$seq)

## AUC
sum(cov)

## ----"coverage-reproduce", eval = FALSE------------------------------------
#  ## Code for reproducing the bottom portion of Figure 8.
#  pdf("base_pair_coverage.pdf", width = 20)
#  par(mar = c(5, 6, 4, 2) + 0.1)
#  plot(cov, type = "o", col = "violetred1", lwd = 10, ylim = c(0, 5),
#       xlab = "Genome", ylab = "Coverage", cex.axis = 2, cex.lab = 3,
#       bty = "n")
#  polygon(c(1, seq_len(length(cov)), length(cov)), c(0, cov, 0),
#          border = NA, density = -1, col = "light blue")
#  points(seq_len(length(cov)), cov, col = "violetred1", type = "o",
#         lwd = 10)
#  dev.off()

## ----Figure9, out.width="100%", fig.align="center", fig.cap = "Area under coverage (AUC). The area under coverage is the sum of the base-pair coverage for all positions in the genome regardless of the annotation. It is the area under the base-level coverage curve shown as the light blue area under the pink curve.", echo = FALSE----
knitr::include_graphics("Figure9.png")

## ----"example_scaled", message = FALSE-------------------------------------
## Check that the number of reads is less than or equal to 40 million
## after scaling.
library("recount")
rse_scaled <- scale_counts(rse_gene_SRP009615, round = FALSE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

## ----Figure10, out.width="100%", fig.align="center", fig.cap = "Exon-exon junctions go beyond the annotation. Reads spanning exon-exon junctions are highlighted and compared against the annotation. Three of them match the annotated junctions, but one (blue and orange read) spans an un-annotated exon-exon junction with the left end matching the annotation and the right end hinting at a possible new isoform for this gene (blue and orange isoform).", echo = FALSE----
knitr::include_graphics("Figure10.png")

## ----Figure11, out.width="100%", fig.align="center", fig.cap = "Intron retention events. Some reads might align with known intronic segments of the genome and provide information for exploring intron retention events (pink read). Some might support an intron retention event or a new isoform when coupled with exon-exon junction data (orange read).", echo = FALSE----
knitr::include_graphics("Figure11.png")

## ----Figure12, out.width="100%", fig.align="center", fig.cap = "Exon boundaries. Reads that go beyond the known exon boundaries can inform us of whether the annotated boundaries are correct or if there was a run-off transcription event.", echo = FALSE----
knitr::include_graphics("Figure12.png")

## ----"localPath"-----------------------------------------------------------
## Locate path with the data
library("recountWorkshop")
local_path <- system.file("extdata", "SRP056604", package = "recountWorkshop")
dir(local_path)

## ----"testing locally", echo = FALSE, eval = FALSE-------------------------
#  ## For testing
#  if(FALSE) {
#      local_path <- "/Users/lcollado/Dropbox/Code/recountWorkshop/inst/extdata/SRP056604"
#  }

## ----"download gene"-------------------------------------------------------
library("recount")

## Find the project id by searching abstracts of studies
abstract_search("hippocampi of 4")

## Download the data if it is not there
if(!file.exists(file.path(local_path, "rse_gene.Rdata"))) {
    ## In case you decide to download the data instead of using the
    ## pre-installed data
    local_path <- "SRP056604"
    download_study("SRP056604", type = "rse-gene")
}

## Check that the file was downloaded
file.exists(file.path(local_path, "rse_gene.Rdata"))

## Load the data
load(file.path(local_path, "rse_gene.Rdata"))

## ----colData---------------------------------------------------------------
## One row per sample, one column per phenotype variable
dim(colData(rse_gene))

## Mostly technical variables are included
colnames(colData(rse_gene))

## ----"explore colData"-----------------------------------------------------
## Input reads: number reported by SRA might be larger than number
## of reads Rail-RNA downloaded
colData(rse_gene)[, c("read_count_as_reported_by_sra", "reads_downloaded")]
summary(colData(rse_gene)$proportion_of_reads_reported_by_sra_downloaded)

## AUC information used by scale_counts() by default
head(colData(rse_gene)$auc)

## Alternatively, scale_scounts() can use the number of mapped reads
## and other information
colData(rse_gene)[, c("mapped_read_count", "paired_end", "avg_read_length")]

## ----sharq-----------------------------------------------------------------
## SHARQ tissue predictions: not present for all studies
head(colData(rse_gene)$sharq_beta_tissue)
head(colData(rse_gene_SRP009615)$sharq_beta_tissue)

## ----characteristics-------------------------------------------------------
## GEO information was absent for the SRP056604 data set
colData(rse_gene)[, c("geo_accession", "title", "characteristics")]

## GEO information for the SRP009615 data set
head(colData(rse_gene_SRP009615)$geo_accession)
head(colData(rse_gene_SRP009615)$title, 2)
head(colData(rse_gene_SRP009615)$characteristics, 2)

## Similar but not exactly the same wording used for two different samples
colData(rse_gene_SRP009615)$characteristics[[1]]
colData(rse_gene_SRP009615)$characteristics[[11]]

## For study SRP056604 we have characteristics information
## Note the > symbol in the age for the second sample
colData(rse_gene)$characteristics[[1]]
colData(rse_gene)$characteristics[[2]]

## ----characteristics_process-----------------------------------------------
## Get the case status: either LOAD or control
colData(rse_gene)$case <- factor(
    sapply(colData(rse_gene)$characteristics, function(x)
        ifelse(any(grepl("LOAD", x)), "control", "LOAD"))
)

## Get the sex
colData(rse_gene)$sex <- factor(
    sapply(colData(rse_gene)$characteristics, function(x)
        ifelse(any(grepl("female", x)),
    "female", "male"))
)

## Extract the age. Note that one of them is a bit more complicated.
colData(rse_gene)$age <- sapply(colData(rse_gene)$characteristics,
    function(x) {
        y <- x[grep("age.*(yrs)", x)]
        as.integer(gsub("age.*\\(yrs\\): |>", "", y))
    }
)

## ----"explore case sex age"------------------------------------------------
table(colData(rse_gene)$case, colData(rse_gene)$sex)
table(colData(rse_gene)$case, colData(rse_gene)$age)

## ----"add_predictions"-----------------------------------------------------
## Before adding predictions
dim(colData(rse_gene))

## Add the predictions
rse_gene <- add_predictions(rse_gene)

## After adding the predictions
dim(colData(rse_gene))

## Explore the variables
colData(rse_gene)[, 25:ncol(colData(rse_gene))]

## ----"sex predictions"-----------------------------------------------------
table("Observed" = colData(rse_gene)$sex,
    "Predicted" = colData(rse_gene)$predicted_sex)

## ----"matching sex", out.width="100%", fig.align="center"------------------
## Is the sex matching? TRUE for yes.
colData(rse_gene)$matching_sex <- as.character(colData(rse_gene)$sex) ==
    as.character(colData(rse_gene)$predicted_sex)

## Matching sex vs other variables
table(colData(rse_gene)$matching_sex, colData(rse_gene)$case)
table(colData(rse_gene)$matching_sex, colData(rse_gene)$age)
boxplot(colData(rse_gene)$mapped_read_count ~ colData(rse_gene)$matching_sex,
    ylab = "Mapped read count", xlab = "Matching sex")

## ----"sra_run_table"-------------------------------------------------------
## Save the information from 
## https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP056604
## to a table. We saved the file as SRP056604/SraRunTable.txt.
file.exists(file.path(local_path, "SraRunTable.txt"))

## Read the table
sra <- read.table(file.path(local_path, "SraRunTable.txt"),
    header = TRUE, sep = "\t")

## Explore it
head(sra)

## We will remove some trailing '_s' from the variable names
colnames(sra) <- gsub("_s$", "", colnames(sra))

## Choose some variables we want to add
sra_vars <- c("braak_stage", "gender", "tissue")

## Re-organize the SRA table based on the SRA Run ids we have
sra <- sra[match(colData(rse_gene)$run, sra$Run), ]

## Double check the order
identical(colData(rse_gene)$run, as.character(sra$Run))

## Append the variables of interest
colData(rse_gene) <- cbind(colData(rse_gene), sra[, sra_vars])

## Final dimensions
dim(colData(rse_gene))

## Explore result
colData(rse_gene)[, 38:ncol(colData(rse_gene))]

## The sex information from the 'characteristics' and the
## SRA Run Selector match
table(colData(rse_gene)$gender, colData(rse_gene)$sex)

## ----"scale_counts"--------------------------------------------------------
## Scale counts
rse_gene_scaled <- scale_counts(rse_gene)

## To highlight that we scaled the counts
rm(rse_gene)

## ----"deseq2_de", out.width="100%", fig.align="center"---------------------
library("DESeq2")

## Specify design and switch to DESeq2 format
dds <- DESeqDataSet(rse_gene_scaled, ~ sex + age + case)

## Perform DE analysis
dds <- DESeq(dds, test = "LRT", reduced = ~ sex + age, fitType = "local")
res <- results(dds, alpha = 0.01)

## Explore results
plotMA(res, main="DESeq2 results for SRP056604")

## ----"ideal volcano"-------------------------------------------------------
## Make a volcano plot
library("ideal")
plot_volcano(res, FDR = 0.01)

## ----"create_report", message = FALSE, warning = FALSE, results = "hide", eval = FALSE----
#  ## Make a report with the results
#  library("regionReport")
#  DESeq2Report(dds, res = res, project = "SRP056604",
#      intgroup = c("sex", "case"), outdir = ".",
#      output = "SRP056604_main-results")

## ----"browse_report", eval = FALSE-----------------------------------------
#  browseURL("SRP056604_main-results.html")

## ----"go_analysis", out.width="100%", fig.align="center"-------------------
library("clusterProfiler")
library("org.Hs.eg.db")

## Remember that dds had ENSEMBL ids for the genes
ensembl <- gsub("\\..*", "", rownames(dds))
head(ensembl)

## Not all genes have a p-value
table(!is.na(res$padj))

## Perform enrichment analysis for Biological Process (BP)
## Note that the argument is keytype instead of keyType in Bioconductor 3.5
enrich_go <- enrichGO(gene = ensembl[which(res$padj < 0.05)],
    OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
    universe = ensembl[!is.na(res$padj)])

## Visualize enrichment results
dotplot(enrich_go)

## ----"secondary analysis", out.width="100%", fig.align="center"------------
## DE analysis checking what genes are different by matching sex
dds2 <- DESeqDataSet(rse_gene_scaled, ~ matching_sex)
dds2 <- DESeq(dds2, test = "LRT", reduced = ~ 1, fitType = "local")
res2 <- results(dds2, alpha = 0.01)

## Visually inspect results
plotMA(res2, main="DESeq2 results for SRP056604 - sex predictions")

## Lets add gene symbols to the volcano plot, they are stored in
## the rowRanges slot of the rse object
rowRanges(rse_gene_scaled)
res2$symbol <- sapply(rowRanges(rse_gene_scaled)$symbol, "[[", 1)

## Select some DE genes
intgenes <- res2$symbol[which(res2$padj < 0.0005)]
intgenes <- intgenes[!is.na(intgenes)]

## Make the volcano plot
plot_volcano(res2, FDR = 0.01, intgenes = intgenes)

## ----"create_report2", message = FALSE, warning = FALSE, results = "hide", eval = FALSE----
#  DESeq2Report(dds2, res = res2, project = "SRP056604 - matching sex",
#      intgroup = c("sex", "predicted_sex"), outdir = ".",
#      output = "SRP056604_sex-predictions")

## ----"matching sex de genes chr"-------------------------------------------
sort(table(seqnames(rowRanges(rse_gene_scaled))[which(res2$padj < 0.01)]))

## ----"exon_de_analysis", out.width="100%", fig.align="center"--------------
## Download the data if it is not there
if(!file.exists(file.path(local_path, "rse_exon.Rdata"))) {
    ## In case you decide to download the data instead of using the
    ## pre-installed data
    local_path <- "SRP056604"
    download_study("SRP056604", type = "rse-exon")
}

## Load the data
load(file.path(local_path, "rse_exon.Rdata"))

## Scale and add the metadata (it is in the same order)
identical(colData(rse_exon)$run, colData(rse_gene_scaled)$run)
colData(rse_exon) <- colData(rse_gene_scaled)
rse_exon_scaled <- scale_counts(rse_exon)
## To highlight that we scaled the counts
rm(rse_exon)

## Filter lowly expressed exons: reduces the object size
## and run time
filter_exon <- rowMeans(assays(rse_exon_scaled)$counts) > 5
round(table(filter_exon) / length(filter_exon) * 100, 2)

## Perform the filtering and change default names
rse_e <- rse_exon_scaled[filter_exon, ]
rowRanges(rse_e)$gene_id <- rownames(rse_e)
rownames(rse_e) <- paste0("exon_", seq_len(nrow(rse_e)))

## Create DESeq2 object for the exon data
dds_exon <- DESeqDataSet(rse_e, ~ sex + age + case)

## Perform DE analysis
dds_exon <- DESeq(dds_exon, test = "LRT", reduced = ~ sex + age,
    fitType = "local")
res_exon <- results(dds_exon, alpha = 0.01)

## Explore results
plotMA(res_exon, main="DESeq2 results for SRP056604 - exon level")
plot_volcano(res_exon, FDR = 0.01)

## ----"gene_exon", out.width="100%", fig.align="center"---------------------
## Get the gene ids for genes that are DE at the gene level or that have at
## least one exon with DE signal.
genes_w_de_exon <- unique(rowRanges(rse_e)$gene_id[which(res_exon$padj < 0.01)])
genes_de <- rownames(rse_gene_scaled)[which(res$padj < 0.01)]

## Make a venn diagram
library("gplots")
vinfo <- venn(list("genes" = genes_de, "exons" = genes_w_de_exon),
    names = c("genes", "exons"), show.plot = FALSE) 
plot(vinfo) +
    title("Genes with DE signal: at the gene and exon levels")

## ----"gene_exon_match", out.width="100%", fig.align="center"---------------
## Keep only the DE exons that are from a gene that is also DE
top_exon_de <- res_exon[intersect(which(res_exon$padj < 0.01),
    which(rowRanges(rse_e)$gene_id %in% 
    attr(vinfo, "intersections")[["genes:exons"]])), ]
## Add the gene id
top_exon_de$gene_id <- rowRanges(rse_e)$gene_id[match(rownames(top_exon_de),
    rownames(rse_e))]
    
## Find the fold change that is the most extreme among the DE exons of a gene
exon_max_fc <- tapply(top_exon_de$log2FoldChange, top_exon_de$gene_id,
    function(x) { x[which.max(abs(x))] })

## Keep only the DE genes that match the previous selection
top_gene_de <- res[match(names(exon_max_fc), rownames(res)), ]

## Make the plot
plot(top_gene_de$log2FoldChange, exon_max_fc, pch = 20,
    col = adjustcolor("black", 1/2),
    ylab = "Most extreme log FC at the exon level among DE exons",
    xlab = "Log fold change (FC) at the gene level",
    main = "DE genes with at least one DE exon")
abline(a = 0, b = 1, col = "red")
abline(h = 0, col = "grey80")
abline(v = 0, col = "grey80")

## ----"de_genes_by_chr"-----------------------------------------------------
sort(table(seqnames(rowRanges(rse_gene_scaled)[which(res$padj < 0.01)])),
    decreasing = TRUE)

## ----"identify regions", eval = .Platform$OS.type != "windows"-------------
## Define expressed regions for study SRP056604, only for chromosome 12
regions <- expressed_regions("SRP056604", "chr12", cutoff = 5L,
    maxClusterGap = 3000L, outdir = local_path)
    
## Explore the resulting expressed regions
regions
summary(width(regions))
table(width(regions) >= 60)

## Keep only the ones that are at least 60 bp long
regions <- regions[width(regions) >= 60]
length(regions)

## ----"build_rse_ER", eval = .Platform$OS.type != "windows"-----------------
## Compute coverage matrix for study SRP056604, only for chromosome 12
## Takes about 45 seconds with local data
## and about 70 seconds with data from the web
system.time(rse_er <- coverage_matrix("SRP056604", "chr12", regions,
    chunksize = length(regions), outdir = local_path))

## Use the expanded metadata we built for the gene model
colData(rse_er) <- colData(rse_gene_scaled)

## Normally, we would scale the counts but that leads
## to very small numbers in this example.
colSums(assays(scale_counts(rse_er))$counts)
summary(rowMeans(assays(scale_counts(rse_er))$counts))

## So instead of scaling, we will transform to integers
## keeping numbers around the 1 million counts
colSums(assays(rse_er)$counts)
assays(rse_er)$counts <- round(assays(rse_er)$counts, 0)
colSums(assays(rse_er)$counts)

## ----"er_de_analysis", eval = .Platform$OS.type != "windows", out.width="100%", fig.align="center"----
## Define DESeq2 object
dds_er <- DESeqDataSet(rse_er, ~ sex + age + case)
dds_er <- DESeq(dds_er, test = "LRT", reduced = ~ sex + age, fitType = "local")
res_er <- results(dds_er, alpha = 0.05)

## Visually inspect results
plotMA(res_er, main="DESeq2 results for SRP056604 - DERs")
plot_volcano(res_er, FDR = 0.05)

## ----"create_report3", message = FALSE, warning = FALSE, results = "hide", eval = FALSE----
#  DESeq2Report(dds_er, res = res_er, project = "SRP056604 - DERs",
#      intgroup = c("sex", "case"), outdir = ".",
#      output = "SRP056604_DERs")

## ----"sort_qvalue", eval = .Platform$OS.type != "windows"------------------
## Sort regions by q-value
regions_by_padj <- regions[order(res_er$padj, decreasing = FALSE)]

## Look at the top 10
regions_by_padj[1:10]
width(regions_by_padj[1:10])

## ----"find_bws", eval = .Platform$OS.type != "windows"---------------------
## Construct the list of BigWig URLs
## They have the following form:
## http://duffel.rail.bio/recount/
## project id
## /bw/
## sample run id
## .bw
bws_web <- paste0("http://duffel.rail.bio/recount/SRP056604/bw/",
    colData(rse_er)$bigwig_file)

## Note that they are also present in the recount_url data.frame
bws_url <- recount_url$url[match(colData(rse_er)$bigwig_file,
    recount_url$file_name)]
identical(bws_web, bws_url)

## Local bigwigs
bws <- file.path(local_path, "bw", colData(rse_er)$bigwig_file)
all(file.exists(bws))

## Use the sample run ids as the sample names
names(bws) <- colData(rse_er)$run

## ----"add_padding", eval = .Platform$OS.type != "windows"------------------
## Add 100 bp padding on each side
regions_resized <- resize(regions_by_padj[1:10],
    width(regions_by_padj[1:10]) + 200, fix = "center")

## ----"regionCov", eval = .Platform$OS.type != "windows"--------------------
## Get the bp coverage data for the plots
library("derfinder")
regionCov <- getRegionCoverage(regions = regions_resized, files = bws,
    targetSize = 40 * 1e6 * 100, totalMapped = colData(rse_er)$auc,
    verbose = FALSE)

## ----"gencode_txdb", eval = .Platform$OS.type != "windows"-----------------
## Import the Gencode v25 hg38 gene annotation
library("rtracklayer")

## If accessing from the web
if(FALSE) {
    gencode_v25_hg38 <- import(paste0(
        "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/",
        "gencode.v25.annotation.gtf.gz"))  
}

## Using the file locally, included in the workshop package
gencode_v25_hg38 <- import(system.file("extdata",
    "gencode.v25.annotation.gtf.gz", package = "recountWorkshop"))
            
## Keep only the chr12 info
gencode_v25_hg38 <- keepSeqlevels(gencode_v25_hg38, "chr12",
    pruning.mode = "coarse")

## Get the chromosome information for hg38
library("GenomicFeatures")
chrInfo <- getChromInfoFromUCSC("hg38")
chrInfo$chrom <- as.character(chrInfo$chrom)
chrInfo <- chrInfo[chrInfo$chrom %in% seqlevels(regions), ]
chrInfo$isCircular <- FALSE

## Assign the chromosome information to the object we will use to
## create the txdb object
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
    genome = "hg38"))
seqinfo(gencode_v25_hg38) <- si

## Switch from Gencode gene ids to Ensembl gene ids
gencode_v25_hg38$gene_id <- gsub("\\..*", "", gencode_v25_hg38$gene_id)

## Create the TxDb object
gencode_v25_hg38_txdb <- makeTxDbFromGRanges(gencode_v25_hg38)

## Explore the TxDb object
gencode_v25_hg38_txdb

## ----"bump_ann", eval = .Platform$OS.type != "windows"---------------------
library("bumphunter")
## Annotate all transcripts for gencode v25 based on the TxDb object
## we built previously.
ann_gencode_v25_hg38 <- annotateTranscripts(gencode_v25_hg38_txdb,
    annotationPackage = "org.Hs.eg.db",
    mappingInfo = list("column" = "ENTREZID", "keytype" = "ENSEMBL",
    "multiVals" = "first"))
    
## Annotate the regions of interest
## Note that we are using the original regions, not the resized ones
nearest_ann <- matchGenes(regions_by_padj[1:10], ann_gencode_v25_hg38)

## ----"make_gs", eval = .Platform$OS.type != "windows"----------------------
## Create the genomic state object using the gencode TxDb object
gs_gencode_v25_hg38 <- makeGenomicState(gencode_v25_hg38_txdb,
    chrs = seqlevels(regions))
    
## Annotate the original regions
regions_ann <- annotateRegions(regions_resized,
    gs_gencode_v25_hg38$fullGenome)

## ----"region_plots", eval = .Platform$OS.type != "windows", out.width="100%", fig.align="center"----
library("derfinderPlot")
plotRegionCoverage(regions = regions_resized, regionCoverage = regionCov, 
   groupInfo = colData(rse_er)$case,
   nearestAnnotation = nearest_ann, 
   annotatedRegions = regions_ann,
   txdb = gencode_v25_hg38_txdb,
   scalefac = 1, ylab = "Coverage (RP40M, 100bp)",
   ask = FALSE, verbose = FALSE)

## ----sessionInfo----------------------------------------------------------------------------------
## Pandoc information
rmarkdown::pandoc_version()

## Time for reproducing this workflow, in minutes
round(proc.time()[3] / 60, 1)

options(width = 100)
library("devtools")
session_info()

