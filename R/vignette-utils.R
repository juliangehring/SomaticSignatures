scaSNVRanges <- function(chrs = hsAutosomes()) {

    sca_all = SomaticCancerAlterations::scaLoadDatasets()
    sca_metadata = SomaticCancerAlterations::scaMetadata()

    sca_merge = unlist(sca_all)
    short_names = gsub("(.*)_(.*)", "\\1", rownames(sca_metadata))
    names(sca_merge) = sca_merge$study = factor(rep(short_names, times = elementLengths(sca_all)))
    
    sca_merge = sca_merge[ sca_merge$Variant_Type %in% "SNP" ]
    sca_merge = keepSeqlevels(sca_merge, chrs)
    
    sca_vr = VRanges(
        seqnames(sca_merge),
        ranges(sca_merge),
        ref = sca_merge$Reference_Allele,
        alt = sca_merge$Tumor_Seq_Allele2,
        sampleNames = sca_merge$Patient_ID,
        seqinfo = seqinfo(sca_merge),
        study = sca_merge$study)
    sca_vr = ucsc(sca_vr)
    
    return(sca_vr)
}
