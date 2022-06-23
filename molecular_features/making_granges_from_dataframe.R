library(PharmacoGx)
library(SummarizedExperiment)
library(GenomicRanges)
library(data.table)

# -- Get our example data
data(GDSCsmall)
SE <- molecularProfilesSlot(GDSCsmall)$rna
features <- as.data.frame(rowData(SE))
setDT(features)  # to data.table by reference

# -- Get our annotations
gencodev33_url <- "https://github.com/BHKLAB-Pachyderm/Annotations/raw/master/Gencode.v33.annotation.RData"
genecode_file <- file.path(tempdir(), "gencodev33.rda")
download.file(
    gencodev33_url,
    destfile=genecode_file
)
(object_names <- load(genecode_file))

# get features_gene table
gencodev33 <- as.data.table(get(object_names[2]))
# remove versions to match the
gencodev33[, gene_id_no_ver := gsub("\\..*$", "", gene_id)]

# -- Attach to our feature data
ranged_features <- merge.data.table(
    features,
    gencodev33,
    by.x="EnsemblGeneId",
    by.y="gene_id_no_ver",
    sort=FALSE,
    all.x=TRUE
)

# handle unmapped features
ranged_features[
    is.na(start),
    c("start", "end", "length", "strand") := list(-1, -1, 0, "*")
]
# drop Y chromosome version if duplicated genes
# (lexically X is before Y so first always select X)
ranged_features <- ranged_features[, lapply(.SD, first), by=gene_id]

# make sure duplicated were only due to chromosomes
stopifnot(all(ranged_features$rownames == features$features))

# build the GRanges object
row_ranges <- makeGRangesFromDataFrame(
    ranged_features,
    keep.extra.columns=TRUE  # retain metadata
)
names(row_ranges) <- row_ranges$rownames

RSE <- SummarizedExperiment(
    assays=assays(SE),
    rowRanges=row_ranges,
    colData=colData(SE),
    metadata=metadata(SE)
)