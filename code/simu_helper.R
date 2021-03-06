split_coords <- function(x) {
    # split locations in the format
    # chr:start-end and return chrom,
    # start and ends as separate values
    loc <- strsplit(x, ':')[[1]]
    chrom <- loc[1]
    pos <- strsplit(loc[2], '-')[[1]]
    start <- pos[1]; end <- strsplit(pos[2], '\\(')[[1]][1]
    return(c(chrom, start, end))
}

get_loc_data <- function(locs) {
    # apply coordinate splitting to whole
    # data frame and return results
    loc_data <- data.frame(t(sapply(locs, split_coords)))
    loc_data[,2:3] <- apply(loc_data[,2:3], 2, as.numeric)
    colnames(loc_data) <- c('chr', 'start', 'end')
    return(loc_data)
}

convert_chrom <- function(chrom, add_chr = TRUE) {
    # add or remove 'chr' prefix from chrom
    if(add_chr) {
        new_chrom <- ifelse(chrom == 'MT', 'chrM', paste0('chr', chrom))
    } else {
        new_chrom <- ifelse(chrom == 'chrM', 'MT', gsub('chr', '', chrom))
    }
    return(new_chrom)
}

get_vartype_from_gene <- function(gene, fus_truth, tsv_nsv_truth) {
    # takes a gene and the truth data from the simulations,
    # and returns  the variant type of the simulated gene
    tsv_hit <- tsv_nsv_truth[tsv_nsv_truth$gene %in% gene,]
    fus_hit <- fus_truth[fus_truth$gene1 %in% gene | fus_truth$gene2 %in% gene ,]

    if (nrow(tsv_hit) > 0) {
        return(tsv_hit$vartype[1])
    } else if (nrow(fus_hit) > 0) {
        return(fus_hit$fusion_type[1])
    } else {
        return("FP")
    }
}

is_classification_correct <- function(row) {
    # returns TRUE if MINTIE's classification
    # is plausible given a truth dataframe
    # containing the variant type and the
    # matching truth vartype of the simulated gene
    vtruth <- row['truth_vartype']
    vtype <- row['variant_type']

    if (vtype == 'UN') {
        return(NA)
    } else if (vtruth == vtype) {
        return(TRUE)
    } else if (vtruth == 'canonical_fusion' & vtype == c('FUS')) {
        return(TRUE)
    } else if (vtruth == 'unpartnered_fusion' & vtype %in% c('FUS', 'NE')) {
        return(TRUE)
    } else if (vtruth == 'INS_fusion' & vtype %in% c('FUS', 'INS')) {
        return(TRUE)
    } else if (vtruth == 'EE_fusion' & vtype %in% c('FUS', 'EE')) {
        return(TRUE)
    } else if (vtruth == 'NE_fusion' & vtype %in% c('FUS', 'NE')) {
        return(TRUE)
    } else if (vtruth == 'DEL' & vtype == 'NEJ') {
        return(TRUE)
    } else if (vtruth == 'ITD' & vtype == 'INS') {
        return(TRUE)
    } else if (vtruth %in% c('INV', 'PTD') & vtype == 'IGR') {
        return(TRUE)
    } else if (vtruth == 'US' & vtype == 'AS') {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

get_barnacle_loc <- function(row, side_A = TRUE) {
    # return separated chrom, start and end
    # from barnacle data for the given side
    side <- ifelse(side_A, 'ALIGN_A:', 'ALIGN_B')
    tmp <- strsplit(row, side)[[1]][2]
    tmp <- strsplit(tmp, '=')[[1]][2]
    tmp <- strsplit(tmp, '\\(')[[1]][1]
    tmp <- strsplit(tmp, ':')[[1]]
    chrom <- tmp[1]
    pos1 <- strsplit(tmp[2], '-')[[1]][1]
    pos2 <- strsplit(tmp[2], '-')[[1]][2]
    return(c(chrom, pos1, pos2))
}

get_barnacle_grx <- function(barnacle_results, side_A = TRUE) {
    # return Genomic Ranges object from barnacle
    # results data for the given sidee.
    bn_locs <- data.frame(t(sapply(barnacle_results, get_barnacle_loc, side_A=side_A)))
    rownames(bn_locs) <- 1:nrow(bn_locs)

    # make sure end and start are the right way around
    loc1 <- apply(bn_locs[,2:3], 1, function(x){min(as.numeric(x))})
    loc2 <- apply(bn_locs[,2:3], 1, function(x){max(as.numeric(x))})
    barnacle_grx <- GRanges(seqnames = bn_locs$X1, ranges = IRanges(start = loc1, end = loc2))

    return(barnacle_grx)
}

get_granges <- function(truth_df, convert_chrom = FALSE, add_chr = TRUE) {
    # get Genomic Ranges from truth data
    # the convert_chrom option will add or subtract 'chr'
    # prefix from  the locations (specified by 'add_chr')
    if('loc1' %in% colnames(truth_df)) {
        loc1 <- get_loc_data(truth_df$loc1)
        loc2 <- get_loc_data(truth_df$loc2)
        if (convert_chrom) {
            loc1$chr <- sapply(loc1$chr, convert_chrom, add_chr = add_chr)
            loc2$chr <- sapply(loc2$chr, convert_chrom, add_chr = add_chr)
        }
        grx1 <- GRanges(seqnames = loc1$chr,
                        ranges = IRanges(start = loc1$start,
                                         end = loc1$end))
        grx2 <- GRanges(seqnames = loc2$chr,
                        ranges = IRanges(start = loc2$start,
                                         end = loc2$end))
        return(list(grx1, grx2))
    } else if('loc' %in% colnames(truth_df)) {
        loc <- get_loc_data(truth_df$loc)
        if (convert_chrom) {
            loc$chr <- sapply(loc$chr, convert_chrom, add_chr = add_chr)
        }
        grx <- GRanges(seqnames = loc$chr,
                       ranges = IRanges(start = loc$start,
                                        end = loc$end))
        return(grx)
    } else if('chrom' %in% colnames(truth_df)) {
        if (convert_chrom) {
            truth_df$chrom <- sapply(truth_df$chrom, convert_chrom, add_chr = add_chr)
        }
        grx <- GRanges(seqnames = truth_df$chrom,
                       ranges = IRanges(start = truth_df$start,
                                        end = truth_df$end),
                       genes = truth_df$gene)
        return(grx)
    } else {
        stop('No coordinate info detected.')
    }
}

get_arriba_genes <- function(genes) {
    # split and return all genes found in Arriba's results
    genes <- strsplit(genes, ',')[[1]]
    genes <- gsub('\\([0-9]+\\)', '', genes)
    return(genes)
}

get_kissplice_hits <- function(bam_loc, grx) {
    # get hits from Genomic Ranges object from the KisSplice
    # bam file
    param <- ScanBamParam(which = grx, what = c("pos"))
    sbam <- scanBam(bam_loc, param = param)
    hits <- sapply(sbam, function(x){length(x$pos) > 0})
    return(hits)
}

get_hits <- function(loc1_grx, loc2_grx, fus1_grx_truth, fus2_grx_truth, tsv_nsv_grx_truth, bgenes_grx) {
    # Takes genomic ranges from the caller, as well as the
    # fusion, TSV and NSV results and returns the hits.

    # check all combinations of overlaps for fusions
    olaps1 <- suppressWarnings(findOverlaps(loc1_grx, fus1_grx_truth))
    olaps2 <- suppressWarnings(findOverlaps(loc2_grx, fus2_grx_truth))
    olaps3 <- suppressWarnings(findOverlaps(loc2_grx, fus1_grx_truth))
    olaps4 <- suppressWarnings(findOverlaps(loc1_grx, fus2_grx_truth))

    fus_truth_hits <- Reduce(union, Map(subjectHits, list(olaps1, olaps2, olaps3, olaps4)))
    caller_hits <- Reduce(union, Map(queryHits, list(olaps1, olaps2, olaps3, olaps4)))

    # check overlaps for TSVs and NSVs
    olaps1 <- suppressWarnings(findOverlaps(loc1_grx, tsv_nsv_grx_truth))
    olaps2 <- suppressWarnings(findOverlaps(loc2_grx, tsv_nsv_grx_truth))

    tsv_nsv_truth_hits <- union(subjectHits(olaps1), subjectHits(olaps2))
    caller_hits <- union(caller_hits, union(queryHits(olaps1), queryHits(olaps2)))

    # check for overlaps in background genes
    olaps1 <- suppressWarnings(findOverlaps(loc1_grx, bgenes_grx))
    olaps2 <- suppressWarnings(findOverlaps(loc2_grx, bgenes_grx))
    bg_hits <- union(queryHits(olaps1), queryHits(olaps2))

    results <- list(fus_truth_hits, tsv_nsv_truth_hits, caller_hits, bg_hits)
    names(results) <- c('fus_truth_hits', 'tsv_nsv_truth_hits', 'caller_hits', 'bg_hits')
    return(results)
}

get_hits_oneloc <- function(loc_grx, fus1_grx_truth, fus2_grx_truth, tsv_nsv_grx_truth, bgenes_grx) {
    # check all combinations of overlaps for fusions
    olaps1 <- suppressWarnings(findOverlaps(loc_grx, fus1_grx_truth))
    olaps2 <- suppressWarnings(findOverlaps(loc_grx, fus2_grx_truth))

    fus_truth_hits <- union(subjectHits(olaps1), subjectHits(olaps2))
    caller_hits <- union(queryHits(olaps1), queryHits(olaps2))

    # check overlaps for TSVs and NSVs
    olaps <- suppressWarnings(findOverlaps(loc_grx, tsv_nsv_grx_truth))

    tsv_nsv_truth_hits <- subjectHits(olaps)
    caller_hits <- union(caller_hits, queryHits(olaps))

    # check for overlaps in background genes
    olaps <- suppressWarnings(findOverlaps(loc_grx, bgenes_grx))
    bg_hits <- queryHits(olaps)

    results <- list(fus_truth_hits, tsv_nsv_truth_hits, caller_hits, bg_hits)
    names(results) <- c('fus_truth_hits', 'tsv_nsv_truth_hits', 'caller_hits', 'bg_hits')
    return(results)
}

append_results <- function(results, method, found_fus, found_tsv, n_fp) {
    # append results of the given method to the collated
    # results summary table
    results[, method] <- 0
    results[names(found_fus), method] <- as.numeric(found_fus)
    results[names(found_tsv), method] <- as.numeric(found_tsv)
    results["FP", method] <- n_fp
    return(results)
}

get_experiment_results <- function(results, exp_dir) {
    # get results from coverage and dispersion experiments
    files <- list.files(exp_dir, full.names = TRUE)

    for (results_file in files) {
        # load results file
        mintie_results <- read.delim(results_file, sep="\t")

        # extract genes
        mintie_vargenes <- sapply(mintie_results$overlapping_genes, function(x){strsplit(x, "\\||:")[[1]]})
        mvg <- unlist(mintie_vargenes)
        mvg <- mvg[mvg != ""]

        # count found fusions, TSVs and NSVs
        found_fus <- table(fus_truth$fusion_type[fus_truth$gene1 %in% mvg | fus_truth$gene2 %in% mvg])
        found_tsv <- table(tsv_nsv_truth$vartype[tsv_nsv_truth$gene %in% mvg])

        # count false positives genes
        fp <- sapply(mintie_vargenes, function(x){!any(x %in% var_genes_truth)})
        fp_genes <- unlist(sapply(mintie_results$overlapping_genes[fp], function(x){strsplit(x, "\\||:")}))
        fp_genes <- fp_genes[fp_genes != ""]
        n_fp <- length(unique(fp_genes))

        cname <- str_split(results_file, "/")[[1]] %>%
                    tail(1) %>%
                    gsub("allvars-case-", "", .) %>%
                    gsub("_results.tsv.gz", "", .)
        results <- append_results(results, cname, found_fus, found_tsv, n_fp)
    }
    return(results)
}
