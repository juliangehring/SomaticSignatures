
## Compute the 3mer frequency for the genome and exome of humans
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

k = 3
n = 1e7
toplevel_chrs = paste0("chr", c(1:22, "X", "Y"))

wgs_ranges = as(seqinfo(BSgenome.Hsapiens.UCSC.hg19), "GRanges")
wgs_ranges = keepSeqlevels(wgs_ranges, toplevel_chrs)

wes_ranges = reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
wes_ranges = keepSeqlevels(wes_ranges, toplevel_chrs)
  
## wg
k3wg = kmerFrequency(BSgenome.Hsapiens.UCSC.hg19, n, k, wgs_ranges)
stopifnot(length(k3wg) == 64 + 1) ## all combinations + NNN
## drop the 'NNN'
k3wg = k3wg[names(k3wg) != "NNN"]
k3wg = k3wg / sum(k3wg)

## we
k3we = kmerFrequency(BSgenome.Hsapiens.UCSC.hg19, n, k, wes_ranges)
stopifnot(length(k3wg) == 64) ## all combinations

save(file = "package/SomaticSignatures/data/kmers.rda",
     list = c("k3wg", "k3we"),
     compress = "xz")

## figures
library(ggplot2)
library(reshape2)
dfg = melt(k3wg)
dfg$set = "k3wg"
dfe = melt(k3we)
dfe$set = "k3we"

df = rbind(dfg, dfe)

p = ggplot(df) + geom_bar(aes_string(x = "seq", y = "value", group = "set", fill = "set"), stat = "identity") + facet_grid(. ~ set) + theme_bw() + theme(legend.position = "none")
