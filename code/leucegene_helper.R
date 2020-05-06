collate_vartypes <- function(results) {
    # take results from MINTIE and collate
    # vartypes into the following classes:
    # fusion, NSV, TSV and unknown

    # categories
    fusions <- c("FUS")
    splicevars <- c("AS", "EE", "NE", "NEJ", "PNJ", "RI")
    tsvs <- c("DEL", "INS")

    results$class <- "Unknown"
    results$class[results$variant_type %in% fusions] <- "Fusion"
    results$class[results$variant_type %in% tsvs] <- "Transcribed structural variant"
    results$class[results$variant_type %in% splicevars] <- "Novel splice variant"

    # fix variants marked as fusions as TSVs if they
    # are fusions (i.e. hard-clips) within the same gene
    is_tsv <- results$overlapping_genes %>% str_split(":|\\|") %>%
                lapply(., duplicated) %>%
                lapply(., any) %>%
                unlist()
    results$class[is_tsv] <- "Transcribed structural variant"

    return(results)
}

get_samples_with_kmt2a_sv <- function(results, controls) {
    # return unique samples with KMT2A hard/soft
    # clipped variant transcripts
    samples <- results %>%
                filter(controls == controls,
                       variant_type %in% c("FUS", "UN"),
                       overlapping_genes %like% "KMT2A") %>%
                pull(sample) %>%
                unique()
    return(samples)
}
