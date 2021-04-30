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
            list.files(., full.names=TRUE) %>%
            lapply(., read.delim) %>%
            rbindlist(., fill=TRUE) %>%
            # filter(logFC > logFC_cutoff) %>%
            mutate(controls = control)
        results <- rbind(results, tmp, fill=TRUE)
    }
    results <- results[results$logFC > logFC_cutoff | is.na(results$logFC),]
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

get_results <- function(vdir) {
    # get results from directory (e.g. arriba, squid)
    all_results <- NULL
    rout <- list.files(vdir, full.names = TRUE)
    for(results_file in rout) {
        results <- read.delim(results_file, sep='\t')
        results$sample <- str_split(results_file, '/')[[1]] %>%
            tail(., 1) %>%
            str_split('_') %>%
            unlist() %>%
            head(., 1)
        all_results <- rbind(all_results, results)
    }
    return(all_results)
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
        varfilter <- c("INS", "UN", "IGR")
    } else if (variant_type == "PTD") {
        varfilter <- c("UN", "IGR")
    }
    return(varfilter)
}

get_varfilter_arriba <- function(variant_type) {
    # return what Arriba might classify
    # the input variant type as
    if (variant_type == "fusion") {
        varfilter <- "translocation|inversion"
    } else  if (variant_type == "ITD") {
        varfilter <- "ITD|duplication"
    } else if (variant_type == "PTD") {
        varfilter <- "duplication"
    }
    return(varfilter)
}

get_varfilter_tap <- function(variant_type) {
    # return what TAP might classify
    # the input variant type as
    if (variant_type == "fusion") {
        varfilter <- "fusion"
    } else  if (variant_type == "ITD") {
        varfilter <- "ITD"
    } else if (variant_type == "PTD") {
        varfilter <- "duplication"
    }
    return(varfilter)
}

get_gene_hits <- function(results, gene, return_bool=FALSE) {
    # return rows from MINTIE results
    # containing variant gene
    if (gene == "") {
        return(data.frame())
    }
    vargenes <-  results$overlapping_genes %>%
        str_split("\\||:")
    gene_hit <- lapply(vargenes, function(x) any(x %in% gene)) %>%
        unlist()
    if(return_bool) {
        return(gene_hit)
    } else {
        return(results[gene_hit,])
    }
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

is_variant_in_mintie_results <- function(srx_id, gene1, gene2, vartype, results) {
    # checks if variant is in MINTIE output for the given
    # sample, genes and variant type
    varfilter <- get_varfilter(vartype)
    hits <- filter(results, sample == srx_id) %>%
                filter(variant_type %in% varfilter)
    hits_gene1 <- get_gene_hits(hits, gene1, return_bool = TRUE)
    if (gene2 != "") {
        hits_gene2 <- get_gene_hits(hits, gene2, return_bool = TRUE)
        sc <- (hits_gene1 | hits_gene2) & hits$variant_type == "UN"
        return(any(hits_gene1 & hits_gene2) | any(sc))
    } else {
        return(any(hits_gene1))
    }
}

is_variant_in_squid_results <- function(sample, gene1, gene2, vgx, squid_results) {
    # checks if variant in sample in gene1 and gene2 is present in squid results
    squid <- squid_results[squid_results$sample == sample,]
    if(nrow(squid) == 0) {
        return(FALSE)
    }
    chr1 <- sapply(squid$X..chrom1, convert_chrom, add_chr=FALSE)
    loc1_grx <- GRanges(seqnames = chr1,
                        ranges = IRanges(start = squid$start1, end = squid$end1))

    chr2 <- sapply(squid$chrom2, convert_chrom, add_chr=FALSE)
    loc2_grx <- GRanges(seqnames = chr2,
                        ranges = IRanges(start = squid$start2, end = squid$end2))

    v <- vgx[vgx$genes %in% c(gene1, gene2)]
    olaps1 <- findOverlaps(loc1_grx, v) %>% suppressWarnings()
    olaps2 <- findOverlaps(loc2_grx, v) %>% suppressWarnings()

    if(gene2 == '') {
        return(1 %in% subjectHits(olaps1) & 1 %in% subjectHits(olaps2))
    } else {
        poss1 <- 1 %in% subjectHits(olaps1) & 2 %in% subjectHits(olaps2)
        poss2 <- 1 %in% subjectHits(olaps2) & 2 %in% subjectHits(olaps1)
        return(poss1 | poss2)
    }
}

is_variant_in_results <- function(sample, gene1, gene2, variant_type, caller, results) {
    # checks if variant in sample in gene1 and
    # gene2 is present for the given caller results
    if (caller == "arriba") {
        varfilter <- get_varfilter_arriba(variant_type)
        g1_field <- "X.gene1"; g2_field <- "gene2"
    } else if (caller == "tap") {
        varfilter <- get_varfilter_tap(variant_type)
        g1_field <- "X5.gene"; g2_field <- "X3.gene"
    }
    res <- results[results$sample %in% sample,]
    if (gene2 == '') {
        found <- any(res[, g1_field] %in% gene1 & res[, g2_field] %in% gene1)
    } else {
        found_genes <- any(res[, g1_field] %in% gene1 & res[, g2_field] %in% gene2)
        found <- found_genes | any(res[, g1_field] %in% gene2 & res[, g2_field] %in% gene1)
    }
    return(sum(found) > 0)
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
