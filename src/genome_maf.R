##### Generate a genome maf file #####

source("src/libs.R")

trios_md <- read_tsv("data/metadata/metadata_sequencing.tsv")

trios_maf <- get_ssm_by_samples(
    trios_md$DNAseq_sample_id,
    seq_type = c("genome", "capture"),
    subset_from_merge = FALSE,
    min_read_support = 0
)

write_tsv(trios_maf, "data/maf/genomes_augmented.grch37.maf")
