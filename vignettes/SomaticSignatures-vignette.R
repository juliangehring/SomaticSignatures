set.seed(1)

options(width = 70)

library(knitr)

inlineCode <- function(file, format = c) {
    file_exist = sapply(file, file.exists)
    file = file[file_exist]
    if(length(file) == 0)
        return("")
    style = sapply(file,
        function(file) {
            paste(readLines(file), collapse = "\n")},
        USE.NAMES = FALSE)
    style = sapply(style, format)
    style = paste(style, "\n", collapse = "\n\n")
    return(style)
}

knitrHeader <- function(css, js) {
    header = opts_knit$get("header")
    if(!missing(css) && !identical(css, character())) {
        header["highlight"] = inlineCode(css)
    }
    if(!missing(js) && !identical(js, character())) {
        header["js"] = inlineCode(js, formatInlineJS)
    }
    return(header)
}

base_dir = system.file(package = "SomaticSignatures")
css_path = file.path(base_dir, "css", "bioc.css")

opts_knit$set(self.contained = TRUE,
              upload.fun = image_uri,
              header = knitrHeader(css = css_path))

opts_chunk$set(comment = "  ",
               fig.path = "",
               fig.align = "center",
               out.width = "50%",
               dpi = 300,
               indent = 10,
               cache = FALSE,
               cache.path = "../cache")

knit_hooks$set(fig.cap = function(before, options, envir) {
    if(!before) {
        paste('<p class="caption">',options$fig.cap,"</p>",sep="")
    }
})
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg19)
sca_metadata = scaMetadata()

sca_metadata
sca_data = unlist(scaLoadDatasets())
sca_data = unlist(sca_data)

short_names = gsub("(.*)_(.*)", "\\1", rownames(sca_metadata))
names(sca_data) = sca_data$study = factor(rep(short_names, 
                       times = elementLengths(sca_data)))
sca_data = sca_data[sca_data$Variant_Type %in% "SNP"]
sca_data = keepSeqlevels(sca_data, hsAutosomes())

sca_vr = VRanges(seqnames(sca_data), ranges(sca_data), 
    ref = sca_data$Reference_Allele, alt = sca_data$Tumor_Seq_Allele2, 
    sampleNames = sca_data$Patient_ID, seqinfo = seqinfo(sca_data), 
    study = sca_data$study)
sca_vr = ucsc(sca_vr)

head(sca_vr, 3)
sort(table(sca_vr$study), decreasing = TRUE)
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19, unify = TRUE)
sca_mm = motifMatrix(sca_motifs, group = "study", normalize = TRUE)

head(round(sca_mm, 4))
svg(file="report/p_mutation_spectrum.svg")
plotMutationSpectrum(sca_motifs, "study")
dev.off()
n_sigs = 5

sigs_nmf = identifySignatures(sca_mm, n_sigs, nmfDecomposition)

sigs_pca = identifySignatures(sca_mm, n_sigs, pcaDecomposition)
sigs_nmf
sigs_pca
n_sigs = 2:8

gof_nmf= assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)

gof_pca = assessNumberSignatures(sca_mm, n_sigs, pcaDecomposition)
svg(file="p_gof_nmf.svg")
plotNumberSignatures(gof_nmf)
dev.off()
svg(file="p_gof_pca.svg")
plotNumberSignatures(gof_pca)
dev.off()
library(ggplot2)
svg(file="report/p_nmf_signatures_map.svg")
plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
dev.off()
svg(file="report/p_nmf_signatures.svg")
plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")
dev.off()
svg(file="report/p_nmf_observed.svg")
plotObservedSpectrum(sigs_nmf)
dev.off()
svg(file="report/p_nmf_fitted.svg")
plotFittedSpectrum(sigs_nmf)
dev.off()
svg(file="report/p_nmf_samples_map.svg")
plotSampleMap(sigs_nmf)
dev.off()
svg(file="report/p_nmf_samples.svg")
plotSamples(sigs_nmf)
dev.off()
svg(file="report/p_pca_signatures_map.svg")
plotSignatureMap(sigs_pca) + ggtitle("Somatic Signatures: PCA - Heatmap")
dev.off()
svg(file="report/p_pca_signatures.svg")
plotSignatures(sigs_pca) + ggtitle("Somatic Signatures: PCA - Barchart")
dev.off()
svg(file="report/p_pca_fitted.svg")
plotFittedSpectrum(sigs_pca)
dev.off()
svg(file="report/p_pca_observed.svg")
plotObservedSpectrum(sigs_pca)
dev.off()
clu_motif = clusterSpectrum(sca_mm, "motif")
svg(file="p_cluster_motifs.svg")
library(ggdendro)

p = ggdendrogram(clu_motif, rotate = TRUE)
p
dev.off()
library(sva)
sca_anno = as.data.frame(lapply(sca_metadata, unlist))

model_null = model.matrix(~ 1, sca_anno)

sca_mm_batch = ComBat(sca_mm, batch = sca_anno$Sequence_Source, mod = model_null)
k = 3
n = 1e4
       
hs_chrs = as(seqinfo(BSgenome.Hsapiens.UCSC.hg19), "GRanges")
hs_chrs = keepStandardChromosomes(hs_chrs)

k3_hs_chrs = kmerFrequency(BSgenome.Hsapiens.UCSC.hg19, n, k, hs_chrs)
k3_hs_chrs
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

k = 3
n = 1e4
    
hs_exons = reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
hs_exons = keepStandardChromosomes(hs_exons)

k3_exons = kmerFrequency(BSgenome.Hsapiens.UCSC.hg19, n, k, hs_exons)
data(kmers)
norms = k3wg / k3we
head(norms)

sca_mm_norm = normalizeMotifs(sca_mm, norms)
sca_gbm = sca_motifs[ names(sca_motifs) %in% "gbm"]
svg(file="report/p_rainfall_alteration.svg")
plotRainfall(sca_gbm, group = "alteration", size = 1)
dev.off()
citation("SomaticSignatures")
source("http://bioconductor.org/biocLite.R")
biocLite("SomaticSignatures")
library(VariantAnnotation)
vr = VRanges(
    seqnames = "chr1",
    ranges = IRanges(start = 1000, width = 1),
    ref = "A",
    alt = "C")

vr
sessionInfo()
