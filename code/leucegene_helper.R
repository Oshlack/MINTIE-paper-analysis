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

load_controls_comparison <- function(results_dir, logFC_cutoff = 5) {
    # load MINTIE results from analysis
    # containing sets of controls
    results <- NULL
    controls <- list.files(results_dir)
    for (control in controls) {
        tmp <- str_c(results_dir, control, sep = "/") %>%
            list.files(., full.names = TRUE) %>%
            lapply(., read.delim) %>%
            rbindlist() %>%
            filter(logFC > logFC_cutoff) %>%
            mutate(controls = control)
        results <- rbind(results, tmp)
    }
    return(results)
}

get_results_by_gene <- function(results) {
    # split MINTIE's results so that each unique
    # gene in the overlapping genes field appears
    # on a separate record row
    var_genes <- results$overlapping_genes %>%
        str_split("\\||:")

    repeat_rows <- rep(1:nrow(results), sapply(var_genes, length))
    results_by_gene <- data.table(results[repeat_rows,])
    results_by_gene$gene <- unlist(var_genes)

    return(results_by_gene)
}

get_results_summary <- function(results, group_var_name = "controls") {
    # extract variant genes and make summary
    results_by_gene <- get_results_by_gene(results)
    results_summary <- results_by_gene[, length(unique(gene)), by = c("sample", "group_var")]
    results_summary <- results_summary %>% arrange(group_var, V1)

    # reorder samples by total variants
    sample_totals <- data.table(results_summary)[, sum(V1), by = c("sample")] %>%
                        arrange(V1)
    results_summary$sample <- factor(results_summary$sample,
                                      levels = sample_totals$sample)

    # reorder by totals across different group_var
    control_totals <- data.table(results_summary)[, sum(V1), by = c("group_var")] %>%
                        arrange(V1)
    results_summary$group_var <- factor(results_summary$group_var,
                                       levels = control_totals$group_var)
    colnames(results_summary)[2] <- group_var_name

    return(results_summary)
}

get_samples_with_kmt2a_sv <- function(results, control_group) {
    # return unique samples with KMT2A hard/soft
    # clipped variant transcripts
    con_results <- filter(results, controls == control_group)
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

get_gene_stats <- function(bygene, results) {
    # print how many variants across how many
    # samples were found for the given gene
    # expects a gene summary table
    gene_summary <- filter(results, gene == bygene)
    paste("We found",
          gene_summary$var_count %>% sum(),
          "variants across",
          gene_summary$sample %>% unique %>% length(),
          "samples in",
          bygene) %>%
        return()
}

get_chess_genes <- function(file) {
    # takes CHESS reference gene file path,
    # loads it and splits gene names
    chess_genes <- read.delim(file)

    genes <- chess_genes$Gene_Name %>%
        str_split(",")
    repeat_rows <- rep(1:nrow(chess_genes),
                       sapply(genes, length))
    chess_genes <- chess_genes[repeat_rows,]
    chess_genes$gene <- unlist(genes)

    return(chess_genes)
}
