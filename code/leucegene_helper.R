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
    con_results <- filter(results, controls == controls)
    get_samples_with_variant(con_results, "KMT2A", "PTD") %>%
        return()
}

get_varfilter <- function(variant_type) {
    # return what MINTIE might classify
    # the input variant type as
    varfilter <- NULL
    if (variant_type == "fusion") {
        varfilter <- c("FUS", "UN")
    } else  if (variant_type == "ITD") {
        varfilter <- c("INS", "FUS", "UN")
    } else if (variant_type == "PTD") {
        varfilter <- c("FUS", "UN")
    }
    return(varfilter)
}

get_gene_hits <- function(results, gene) {
    # return rows from MINTIE results
    # containing variant gene
    if (gene == "") {
        return(data.frame())
    }
    vargenes <-  results$overlapping_genes %>%
        str_split("\\||:")
    gene_hit <- lapply(vargenes, function(x) any(x %in% gene)) %>%
        unlist()
    return(results[gene_hit,])
}

get_samples_with_variant <- function(results, gene, variant_type) {
    # return unique samples with variant from MINTIE's results
    varfilter <- get_varfilter(variant_type)
    samples <- get_gene_hits(results, gene) %>%
                filter(variant_type %in% varfilter) %>%
                pull(sample) %>%
                unique()
    return(samples)
}

is_variant_in_sample <- function(srx_id, gene1, gene2, vartype, results) {
    # checks if variant is in MINTIE output for the given
    # sample, genes and variant type
    varfilter <- get_varfilter(vartype)
    hits <- filter(results, sample == srx_id) %>%
                filter(variant_type %in% varfilter)
    hits_gene1 <- get_gene_hits(hits, gene1)
    hits_gene2 <- get_gene_hits(hits, gene2)

    return(nrow(hits_gene1) > 0 | nrow(hits_gene2) > 0)
}
