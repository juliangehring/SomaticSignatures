plotSpectrum <- function(x, colorby = c("sample", "alteration")) {

    colorby = match.arg(colorby)

    ## reused part of 'meltSignatures'
    w_df = melt(x, varnames = c("motif", "sample"))
    w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
    w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)

    p = ggplot(w_df)
    p = p + geom_bar(aes_string(x = "context", y = "value", fill = colorby),
        stat = "identity", position = "identity")
    p = p + facet_grid(sample ~ alteration)
    p = p + theme_ss() + theme_small_axis()
    p = p + theme(legend.position = "none")
    p = p + scale_fill_brewer(palette = "Set3")
    p = p + xlab("Motif") + ylab("Contribution")

    return(p)
}


plotMutationSpectrum <- function(vr, group, colorby = c("sample", "alteration"), normalize = TRUE) {

    m = motifMatrix(vr, group, normalize = normalize)
    
    p = plotSpectrum(m, colorby)

    return(p)
}


plotObservedSpectrum <- function(s, colorby = c("sample", "alteration")) {

    return(plotSpectrum(observed(s), colorby))
}


plotFittedSpectrum <- function(s, colorby = c("sample", "alteration")) {

    return(plotSpectrum(fitted(s), colorby))
}


plotSignatureMap <- function(s) {

    w_df = meltSignatures(signatures(s))

    p = ggplot(w_df)
    p = p + geom_tile(aes_string(y = "motif", x = "signature", fill = "value"))
    p = p + scale_fill_gradient2(name = "", space = "Lab")
    p = p + theme_ss() + theme_small_axis(x = FALSE)
    p = p + xlab("Signature") + ylab("Motif")
  
    return(p)
}


plotSignatures <- function(s, normalize = FALSE, percent = FALSE) {

    h = signatures(s)
    if(normalize) {
        h = t(t(h) / colSums(h))
        if(percent) {
            h = h * 100
        }
    }
    w_df = meltSignatures(h)

    p = ggplot(w_df)
    p = p + geom_bar(aes_string(x = "context", y = "value", fill = "alteration"),
        stat = "identity", position = "identity")
    p = p + facet_grid(signature ~ alteration)
    p = p + theme_ss() + theme_small_axis()
    p = p + theme(legend.position = "none")
    p = p + scale_fill_brewer(palette = "Set3")
    p = p + xlab("Motif") + ylab("Contribution")

    return(p)
}


plotSampleMap <- function(s) {

    h_df = melt(samples(s), varnames = c("sample", "signature"))
    h_df$signature = factor(h_df$signature)

    ## lower zlim depending on data:
    ## - 0 for NMF-like methods
    ## - min(value) for others
    zmin = min(h_df$value, na.rm = TRUE)
    zmin = ifelse(zmin >= 0, 0, zmin)

    p = ggplot(h_df)
    p = p + geom_tile(aes_string(y = "sample", x = "signature", fill = "value"), color = "black")
    p = p + scale_fill_gradient2(name = "Contribution", limits = c(zmin, NA), space = "Lab")
    p = p + theme_ss()
    p = p + xlab("Signature") + ylab("Sample")
  
    return(p)
}


plotSamples <- function(s, normalize = FALSE, percent = FALSE) {

    h = samples(s)
    if(normalize) {
        h = h / rowSums(h)
        if(percent) {
            h = h * 100
        }
    }
    w_df = melt(h, varnames = c("sample", "signature"))
    w_df$signature = factor(w_df$signature)

    p = ggplot(w_df)
    p = p + geom_bar(aes_string(x = "sample", y = "value", fill = "signature"),
        color = "black", size = 0.3, stat = "identity", position = "stack")
    p = p + theme_ss()
    p = p + scale_fill_brewer(palette = "Set3")
    p = p + xlab("Sample") + ylab("Signature Contribution")

    return(p)
}


theme_ss <- function() {

    t = theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              axis.text.y = element_text(hjust = 0.5))

    return(t)
}

theme_small_axis <- function(x = TRUE, y = TRUE, size = 6, family = "mono") {
    ## decrease the x/y-axis label size
    template = element_text(size = size, family = family)
    t = theme_ss()
    if(x)
        t = t + theme(axis.text.x = template)
    if(y)
        t = t + theme(axis.text.y = template)
        
    return(t)
}


meltSignatures <- function(x, vars = c("motif", "signature")) {
    
    w_df = melt(x, varnames = vars)
    w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
    w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)
    levels = unique(w_df$signature)
    labels = signatureLabels(length(levels))
    if(all(levels %in% labels))
        levels = labels
    w_df$signature = factor(w_df$signature, levels)

    return(w_df)
}


signatureLabels <- function(n) {
    labels = sprintf("S%d", 1:n)
    return(labels)
}
