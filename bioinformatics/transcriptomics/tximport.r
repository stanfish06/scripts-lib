#!/usr/bin/env Rscript --vanilla

# ========================= Helper functions =========================
# Read and parse tx2gene.tsv
read_transcript_info <- function(tinfo_path, tximport_obj) {
  info <- file.info(tinfo_path)
  if (info$size == 0) {
    stop("tx2gene file is empty")
  }

  transcript_info <- read.csv(
    tinfo_path,
    sep = "\t",
    header = TRUE,
    col.names = c("tx", "gene_id", "gene_name")
  )

  extra <- setdiff(
    rownames(tximport_obj[[1]]),
    as.character(transcript_info[["tx"]])
  )
  transcript_info <- rbind(
    transcript_info,
    data.frame(tx = extra, gene_id = extra, gene_name = extra)
  )
  transcript_info <- transcript_info[
    match(rownames(tximport_obj[[1]]), transcript_info[["tx"]]),
  ]
  rownames(transcript_info) <- transcript_info[["tx"]]

  list(
    transcript = transcript_info,
    gene = unique(transcript_info[, 2:3]),
    tx2gene = transcript_info[, 1:2]
  )
}

# Output and save matrices to csvs
build_table <- function(data, row_meta) {
  data.frame(
    cbind(row_meta, data),
    check.names = FALSE
  )
}
write_table <- function(out_file, file_prefix) {
  file_name <- paste0(file_prefix, ".", out_file$out_file_suffix)
  message("Writing ", file_name)
  write.table(
    build_table(out_file$data, out_file$row_meta),
    file_name,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
# ====================================================================

# =========================== Import data ===========================
# Identify all quant.sf files in current directory
# Assume the directory right above the quant.sf is the condition name
fns <- list.files(
  ".",
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)
if (length(fns) == 0) {
  stop("No files found")
}
message("Found ", length(fns), " quant.sf file(s)")
sample_names <- basename(dirname(fns))
names(fns) <- sample_names
message("Samples: ", paste(sample_names, collapse = ", "))

# tximport is the recommended io library for RNAseq analysis
# Check https://github.com/thelovelab/tximport.git for source code
library(tximport)

# Load in data and transcript-gene mapping
message("Running tximport (type = salmon, txOut = TRUE)...")
txi <- tximport(
  fns,
  type = "salmon",
  txOut = TRUE,
  dropInfReps = FALSE
)
message(
  "tximport complete: ", nrow(txi[["counts"]]), " transcripts x ",
  ncol(txi[["counts"]]), " samples"
)
# --------------------------------------------------
# tx2gene.tsv:
#     transcript_id	gene_id	gene_name
#     ENST00000832824.1	ENSG00000290825.2	DDX11L16
#     ENST00000832825.1	ENSG00000290825.2	DDX11L16
#     ENST00000832826.1	ENSG00000290825.2	DDX11L16
# --------------------------------------------------
message("Reading tx2gene.tsv...")
transcript_info <- read_transcript_info(
  tinfo_path = "tx2gene.tsv",
  tximport_obj = txi
)
message(
  "tx2gene mapping loaded: ", nrow(transcript_info$tx2gene), " transcripts, ",
  nrow(transcript_info$gene), " genes"
)
# ===================================================================

# =========================== Transcript expresssion ===========================
# Setting parameters for writing tables
out_file_prefix <- "salmon.merged"
# Transcript outputs
out_files <- list(
  list(
    data = txi[["abundance"]],
    row_meta = transcript_info$transcript,
    out_file_suffix = "transcript_tpm.tsv"
  ),
  list(
    data = txi[["counts"]],
    row_meta = transcript_info$transcript,
    out_file_suffix = "transcript_counts.tsv"
  ),
  list(
    data = txi[["length"]],
    row_meta = transcript_info$transcript,
    out_file_suffix = "transcript_lengths.tsv"
  )
)
# Note:
# - <txi[[counts]]> is the raw count matrix for transcripts
# - <txi[[abundance]]> is the tpm matrix for transcripts
# - <txi[[length]]> is the transcript lengths in bp (not kb)
# ==============================================================================

# =================================== Gene expresssion ===================================
# Compute gene-level expression data
# Read tximport repo's summarizeToGene.R for more details
tx2gene <- transcript_info$tx2gene
message("Summarizing to gene level (countsFromAbundance = 'no')...")
gi <- summarizeToGene(txi, tx2gene = tx2gene)
message("Summarizing to gene level (countsFromAbundance = 'lengthScaledTPM')...")
gi.ls <- summarizeToGene(
  txi,
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)
message("Summarizing to gene level (countsFromAbundance = 'scaledTPM')...")
gi.s <- summarizeToGene(
  txi,
  tx2gene = tx2gene,
  countsFromAbundance = "scaledTPM"
)
message("Gene-level summarization complete: ", nrow(gi[["counts"]]), " genes")
gene_info <- transcript_info$gene[
  match(rownames(gi[[1]]), transcript_info$gene[["gene_id"]]),
]
rownames(gene_info) <- NULL
out_files <- c(out_files, list(
  list(
    data = gi[["length"]],
    row_meta = gene_info,
    out_file_suffix = "gene_lengths.tsv"
  ),
  list(
    data = gi[["abundance"]],
    row_meta = gene_info,
    out_file_suffix = "gene_tpm.tsv"
  ),
  list(
    data = gi[["counts"]],
    row_meta = gene_info,
    out_file_suffix = "gene_counts.tsv"
  ),
  list(
    data = gi.ls[["counts"]],
    row_meta = gene_info,
    out_file_suffix = "gene_counts_length_scaled.tsv"
  ),
  list(
    data = gi.s[["counts"]],
    row_meta = gene_info,
    out_file_suffix = "gene_counts_scaled.tsv"
  )
))
# Note (see tximport source code):
# - <gi[["length"]]> is the estimated lengths for genes in bp (not kb)
#   + it is computed by doing weighted (by transcript TPM) average of transcript lengths
#     > weightedLength <- rowsum(abundanceMatTx * lengthMatTx, geneId)
#     > lengthMat <- weightedLength / abundanceMat
# - <gi[["abundance"]]> is the tpm matrix for genes
#   + it is just abundanceMat <- rowsum(abundanceMatTx, geneId), sum of transcript TPMs
# - <gi[["counts"]]> is the count matrix for genes
#   + it is just countsMat <- rowsum(countsMatTx, geneId), sum of transcript counts
# - <gi.ls[["counts"]]> is the count matrix for genes derived from the gene abundanceMat
# with gene length correction
#   + it uses the estimated average gene length across samples, so you remove sample-wise
#   gene length bias
#     > countsSum <- colsumfun(countsMat)
#     > newCounts <- abundanceMat * rowMeans(lengthMat)
#     > newSum <- colsumfun(newCounts)
#     > countsMat <- t(t(newCounts) * (countsSum/newSum))
#   + it is the inverse of the TPM formula
#   + it is the count matrix used by deseq2
# - <gi.s[["counts"]]> is the count matrix for genes derived from the gene abundanceMat
# without gene length correction
#   + no multiplication of the average gene length
#     > countsSum <- colsumfun(countsMat)
#     > newCounts <- abundanceMat
#     > newSum <- colsumfun(newCounts)
#     > countsMat <- t(t(newCounts) * (countsSum/newSum))
# - To get TPM-like quantity, use gi[["counts"]] and gi.ls[["counts"]], average
# gi[["length"]] across samples to obtain gene lengths, then follow standard TPM formula
# ========================================================================================


# ===================== Output data =====================
message("Writing ", length(out_files), " output table(s)...")
lapply(out_files, write_table, out_file_prefix)
message("Done.")
# =======================================================
