plot_simu_benchmarking <- function(results, method_name) {
    # plot bar chart results of variants found
    # for a single method
    ggplot(results, aes_string("vartype", method_name, fill = "class")) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() +
        theme(axis.text.x = element_text()) +
        ylab("variants found") +
        xlab("") +
        scale_fill_manual(values = cols) +
        theme(legend.position = "right") +
        coord_flip()
}

plot_all_simu_benchmarks <- function(results) {
    # plot bar chart results of variants found
    # for all methods
    ggplot(results, aes(vartype, value, group = method, fill = class)) +
        geom_bar(stat="identity", position="dodge") +
        theme_bw() +
        theme(axis.text.x = element_text()) +
        facet_wrap(~method, nrow = 3) +
        coord_flip() +
        ylab("variants found") +
        xlab("") +
        scale_fill_manual(values = cols) +
        theme(legend.position = "bottom")
}
