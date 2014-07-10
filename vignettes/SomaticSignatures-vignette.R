
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
               out.width = "100%",
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

library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(stringr)

library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg19)

sca_metadata = scaMetadata()

print(sca_metadata)

sca_all = scaLoadDatasets()

sca_merge = unlist(sca_all)
short_names = str_split_fixed(rownames(sca_metadata), "_", 2)[ ,1]
names(sca_merge) = sca_merge$study = factor(rep(short_names, times = elementLengths(sca_all)))
  
sca_merge = sca_merge[ sca_merge$Variant_Type %in% "SNP" ]
sca_merge = keepSeqlevels(sca_merge, hsAutosomes())

sort(table(sca_merge$study), decreasing = TRUE)

sort(table(sca_merge$Variant_Classification), decreasing = TRUE)

sca_vr = VRanges(
    seqnames(sca_merge),
    ranges(sca_merge),
    ref = sca_merge$Reference_Allele,
    alt = sca_merge$Tumor_Seq_Allele2,
    seqinfo = seqinfo(sca_merge))
sca_vr = ucsc(sca_vr)

head(sca_vr, 3)

sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19, unify = TRUE)

sca_motifs$study = sca_merge$study

head(sca_motifs, 3)

sca_occurrence = motifMatrix(sca_motifs, group = "study", normalize = TRUE)

head(round(sca_occurrence, 4))

plotSamplesObserved(sca_motifs, group = "study")

n_sigs = 5

sigs_nmf = nmfSignatures(sca_occurrence, r = n_sigs)

sigs_pca = pcaSignatures(sca_occurrence, r = n_sigs)

names(sigs_nmf)
  
sapply(sigs_nmf, dim)

head(sigs_nmf$w, 3)

head(sigs_nmf$h, 3)

plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")

plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")


ss = sigs_nmf$w
sn = ss / rowSums(ss)

sx = sigs_nmf
sx$w = sn

plotSampleMap(sigs_nmf)

plotSamples(sigs_nmf)

plotSignatureMap(sigs_pca) + ggtitle("Somatic Signatures: PCA")

plotSignatures(sigs_pca)

plotSampleMap(sigs_pca)

plotSamples(sigs_pca)

library(sva)
library(stringr)

df = as(sca_metadata, "data.frame") ## sample x covariable
pheno = data.frame(s = unlist(df[ ,"Sequence_Source"]), c = unlist(df[ ,"Cancer_Type"]))
rownames(pheno) = str_split_fixed(rownames(pheno), "_", 2)[ ,1]
mod = model.matrix(~ s + c, data = pheno)
mod0 = model.matrix(~ c, data = pheno)

sv = sva(sca_occurrence, mod, mod0, method = "irw")

k = 3
n = 1e5
chrs = "chr1"
    
chr1_ranges = as(seqinfo(BSgenome.Hsapiens.UCSC.hg19), "GRanges")
chr1_ranges = keepSeqlevels(chr1_ranges, chrs)

k3_chr1 = kmerFrequency(BSgenome.Hsapiens.UCSC.hg19, n, k, chr1_ranges)

k3_chr1

head(sca_occurrence)

data(kmers)
norms = k3wg / k3we
head(norms)

sca_norm = normalizeMotifs(sca_occurrence, norms)

head(sca_norm)

plotSamplesObserved(sca_occurrence, group = "study")

plotSamplesObserved(sca_norm, group = "study")

sca_gbm = sca_motifs[ names(sca_motifs) %in% "gbm"]

plotRainfall(sca_gbm, group = "alteration", size = 1)

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

vcs_occurrence = motifMatrix(vcs_motifs, group = "GENE", normalize = TRUE)

head(round(vcs_occurrence, 4))

plotSamplesObserved(vcs_motifs, group = "GENE")

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
