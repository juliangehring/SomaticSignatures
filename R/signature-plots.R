
plotSamplesObserved <- function(s, group = "study") {

    df = as(mcols(s), "data.frame")
    facet_string = paste0(group, " ~ alteration")
    p = ggplot(df) + geom_bar(aes_string(x = "context", fill = "alteration")) +
        facet_grid(facet_string, scales = "free_y") +
            .theme_ss + theme(legend.position = "none") +
                xlab("Motif") + ylab("Frequency")

    return(p)
}


plotSignatureMap <- function(s) {

    w_df = .meltSignatures(s$w)
    p_w = ggplot(w_df) + geom_tile(aes_string(y = "motif", x = "signature",
        fill = "value")) + scale_fill_gradient2(name = "") +
            .theme_ss + xlab("Signature") + ylab("Motif")
  
    return(p_w)
}


plotSignatures <- function(s) {

    w_df = .meltSignatures(s$w)
    p = ggplot(w_df) + geom_bar(aes_string(x = "context", y = "value",
        fill = "alteration"), stat = "identity", position = "identity") +
            facet_grid(signature ~ alteration, scales = "free_y") + .theme_ss +
                theme(legend.position = "none") + xlab("Motif") + ylab("Contribution")
    
    return(p)
}


plotSampleMap <- function(s) {

    h_df = melt(s$h, varnames = c("sample", "signature"))
    h_df$signature = factor(h_df$signature, mixedsort(levels(h_df$signature)))
    p_h = ggplot(h_df) + geom_tile(aes_string(y = "sample", x = "signature",
        fill = "value")) + scale_fill_gradient2(name = "") + .theme_ss +
            xlab("Signature") + ylab("Sample")
  
    return(p_h)
}


plotSamples <- function(s) {

    w_df = melt(s$h, varnames = c("sample", "signature"))
    context_string = str_split_fixed(w_df$sample, pattern = " ", n = 2)
    w_df$alteration = context_string[ ,1]
    w_df$context = context_string[ ,2]
    p = ggplot(w_df) + geom_bar(aes_string(x = "context", y = "value",
        fill = "alteration"), stat = "identity", position = "identity") +
            facet_grid(signature ~ alteration, scales = "free_y")  + .theme_ss +
                theme(legend.position = "none") +
                    xlab("Sample") + ylab("Signature Contribution")

    return(p)
}


.theme_ss <- theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5),
          axis.text = element_text(size = 6, family = "mono"))


.meltSignatures <- function(x, vars = c("motif", "signature")) {
    
    w_df = melt(x, varnames = vars)
    context_string = str_split_fixed(w_df$motif, pattern = " ", n = 2)
    w_df$alteration = context_string[ ,1]
    w_df$context = context_string[ ,2]
    w_df$signature = factor(w_df$signature, mixedsort(levels(w_df$signature)))

    return(w_df)
}
