library(optparse)
suppressMessages(library(dplyr))
library(readr)
library(stringr)

# ---------------------------------------------------------------------------
# Command-line options
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("-s", "--sample_rgx"), type = "character", default = "^([^-]+)",
              help = "regex for sample [default= %default]", metavar = "character"),
  make_option(c("-r", "--group_rgx"), type = "character", default = "^([^_]+)",
              help = "regex for group [default= %default]", metavar = "character"),
  make_option(c("-g", "--geno_rgx"), type = "character", default = "^([^_]+)",
              help = "regex for genotype [default= %default]", metavar = "character"),
  make_option(c("-c", "--cond_rgx"), type = "character", default = "^([^_]+)",
              help = "regex for condition [default= %default]", metavar = "character"),
  make_option(c("-u", "--unit_rgx"), type = "character", default = "^([^_]+)",
              help = "regex for unit [default= %default]", metavar = "character"),
  make_option(c("-i", "--fq_dir"), type = "character",
              help = "directory containing fastq.gz files", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------
# Discover all fastq.gz / fq.gz files in the directory (non-recursive)
# ---------------------------------------------------------------------------
fq_files <- file.path(
  opt$fq_dir,
  list.files(opt$fq_dir, pattern = "\\.(fastq|fq)\\.gz$", recursive = FALSE)
) %>%
  normalizePath() %>%
  unlist()

# ---------------------------------------------------------------------------
# Keep only:
#   - R1 files from paired-end data  (contain -R1_ or -R1.)
#   - files with no R1/R2 marker     (single-end)
# R2 files are intentionally excluded here; they are derived below.
# ---------------------------------------------------------------------------
fq1_files <- fq_files[
  str_detect(basename(fq_files), "-R1[_\\.]") |     # PE: R1 read
  !str_detect(basename(fq_files), "-R[12][_\\.]")   # SE: no R1/R2 marker
]

# ---------------------------------------------------------------------------
# Build the samplesheet data frame
# ---------------------------------------------------------------------------
df <- data.frame(fq1 = fq1_files) %>%
  mutate(
    # Extract metadata from file names using user-supplied regex patterns
    sample    = str_extract(basename(fq1), opt$sample_rgx),
    group     = str_extract(basename(fq1), opt$group_rgx),
    genotype  = str_extract(basename(fq1), opt$geno_rgx),
    condition = str_extract(basename(fq1), opt$cond_rgx),
    unit      = str_extract(basename(fq1), opt$unit_rgx),

    # Derive fq2 path for PE files by replacing -R1 with -R2.
    # Crucially, also check the R2 file exists on disk:
    # files named *-R1.fastq.gz that have no R2 partner are treated as SE.
    fq2 = ifelse(
      str_detect(basename(fq1), "-R1[_\\.]") &
        file.exists(str_replace(fq1, "-R1([_\\.])", "-R2\\1")),
      str_replace(fq1, "-R1([_\\.])", "-R2\\1"),
      NA_character_
    ),

    strandedness = "reverse"
  ) %>%
  arrange(sample)

# ---------------------------------------------------------------------------
# Validations
# ---------------------------------------------------------------------------

# fq1 paths must be unique (no duplicate R1 / SE files)
stopifnot(length(df$fq1) == length(unique(df$fq1)))

# fq2 paths must be unique among paired rows
paired_fq2 <- df$fq2[!is.na(df$fq2)]
stopifnot(length(paired_fq2) == length(unique(paired_fq2)))

# All fq1 files must exist on disk
stopifnot(all(file.exists(df$fq1)))

# All derived fq2 files must exist on disk (skip NA rows = SE samples)
stopifnot(all(is.na(df$fq2) | file.exists(df$fq2)))

# fq1 column must never contain R2 files
stopifnot(all(!str_detect(basename(df$fq1), "-R2[_\\.]")))

# fq2 column must contain R2 files where not NA
stopifnot(all(is.na(df$fq2) | str_detect(basename(df$fq2), "-R2[_\\.]")))

# Every fastq file found in the directory must appear exactly once
# in either fq1 or fq2 — catches missing pairs or untracked files
listed_files <- c(df$fq1, df$fq2[!is.na(df$fq2)])
stopifnot(setequal(normalizePath(listed_files), normalizePath(fq_files)))

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
# Always write all three columns (sample, fq1, fq2).
# Missing fq2 (SE samples) are written as empty string, not NA.
write_csv(df %>% select(sample, fq1, fq2), "samplesheet.csv", na = "")

message("Done. Samplesheet written to samplesheet.csv")
message("  Total samples : ", nrow(df))
message("  Paired-end    : ", sum(!is.na(df$fq2)))
message("  Single-end    : ", sum(is.na(df$fq2)))
