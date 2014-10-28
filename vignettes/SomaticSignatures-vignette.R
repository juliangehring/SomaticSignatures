set.seed(1)

options(width = 70)

library(knitr)

style_sheet = "bioc.css"
style = if(file.exists(style_sheet)) {
    paste(readLines(style_sheet), collapse = "\n")
}
    
opts_knit$set(self.contained = TRUE,
              upload.fun = image_uri,
              header = c(highlight = style))

opts_chunk$set(comment = "  ",
               fig.path = "",
               fig.align = "center",
               out.width = "80%",
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
library(ggplot2)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg19)
sca_metadata = scaMetadata()

print(sca_metadata)
sca_vr = scaSNVRanges()

head(sca_vr, 3)
sort(table(sca_vr$study), decreasing = TRUE)
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19, unify = TRUE)
sca_mm = motifMatrix(sca_motifs, group = "study", normalize = TRUE)

head(round(sca_mm, 4))
svg(file="report/p_mutation_spectrum.svg")
plotMutationSpectrum(sca_motifs, "study")
dev.off()
n_sigs = 5

sigs_nmf = identifySignatures(sca_mm, n_sigs, nmfSignatures)

sigs_pca = identifySignatures(sca_mm, n_sigs, pcaSignatures)
sigs_nmf
sigs_pca
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
svg(file="report/p_pca_observed.svg")
plotObservedSpectrum(sigs_pca)
dev.off()
svg(file="report/p_pca_fitted.svg")
plotFittedSpectrum(sigs_pca)
dev.off()
library(sva)

df = as(sca_metadata, "data.frame") ## sample x covariable
pheno = data.frame(s = unlist(df[ ,"Sequence_Source"]), c = unlist(df[ ,"Cancer_Type"]))
rownames(pheno) = gsub("(.*)_.*", "\\1", rownames(pheno))
mod = model.matrix(~ s + c, data = pheno)
mod0 = model.matrix(~ c, data = pheno)

sv = sva(sca_mm, mod, mod0, method = "irw")
k = 3
n = 1e5
chrs = "chr1"
    
chr1_ranges = as(seqinfo(BSgenome.Hsapiens.UCSC.hg19), "GRanges")
chr1_ranges = keepSeqlevels(chr1_ranges, chrs)

k3_chr1 = kmerFrequency(BSgenome.Hsapiens.UCSC.hg19, n, k, chr1_ranges)

k3_chr1
head(sca_mm)

data(kmers)
norms = k3wg / k3we
head(norms)

sca_norm = normalizeMotifs(sca_mm, norms)

head(sca_norm)
svg(file="report/p_samples_unnorm.svg")
plotSamplesObserved(sca_mm, group = "study")
dev.off()
svg(file="report/p_samples_norm.svg")
plotSamplesObserved(sca_norm, group = "study")
dev.off()
sca_gbm = sca_motifs[ names(sca_motifs) %in% "gbm"]
svg(file="report/p_rainfall_alteration.svg")
plotRainfall(sca_gbm, group = "alteration", size = 1)
dev.off()
library(SomaticSignatures)

library(COSMIC.67)
library(VariantAnnotation)
vcf_path = system.file("vcf", "cosmic_67.vcf.gz", package = "COSMIC.67", mustWork = TRUE)
genes = c("KRAS", "APC", "BRCA1", "BRCA2", "BRAF", "TP53")

vc = ucsc(readVcfAsVRanges(vcf_path, "ncbi37", ScanVcfParam(info = "GENE")))

vcs = vc[ vc$GENE %in% genes & isSNV(vc) ]
head(vcs)
table(vcs$GENE)
library(BSgenome.Hsapiens.UCSC.hg19)

vcs_motifs = mutationContext(ucsc(vcs), BSgenome.Hsapiens.UCSC.hg19, unify = TRUE)
vcs_mm = motifMatrix(vcs_motifs, group = "GENE", normalize = TRUE)

head(round(vcs_mm, 4))
svg(file="report/vcs_samples_observed.svg")
plotSamplesObserved(vcs_motifs, group = "GENE")
dev.off()
vcf_path = COSMIC.67:::cosmicVcfPath()

data(genesymbol, package = "biovizBase")

nice_genes = c("KRAS", "NRAS", "APC", "BRCA1", "BRCA2", "TP53")

roi = sort(ncbi(unstrand(genesymbol[nice_genes])))

param = ScanVcfParam(which = roi, info = "GENE")
vcf = readVcfAsVRanges(vcf_path, "ncbi37", param)
vcf$GENE = factor(sub("_.*", "", vcf$GENE))

vcf <- readVcfAsVRanges(vcf_path, "ncbi37", ScanVcfParam(info = "GENE"))

table(vcf$GENE)
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
