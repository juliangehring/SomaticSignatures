library(testthat)
library(SomaticSignatures)

context("mutect")

mutect_path = system.file("examples", "mutect.tsv", package = "SomaticSignatures")

## random mutations
n = 100
ref = sample(1:4, n, TRUE)
alt = sapply(ref, function(r) sample(setdiff(1:4, r), 1))
vr = VRanges(sample(1:3, n, TRUE),
             IRanges(sample(1:10000, n), width = 1),
             ref = ref, alt = alt,
             totalDepth = rbinom(n, 50, 0.9),
             altDepth = rbinom(n, 40, 0.4), 
             sampleNames = sample(c("T1", "T2"), n, TRUE))

test_that("'readMutect' works", {

    expect_true( file.exists(mutect_path) )

    vr1 = readMutect(mutect_path)
    vr2 = readMutect(mutect_path, strip = TRUE)

    expect_is( vr1, "VRanges" )
    expect_is( vr2, "VRanges" )
})


context("granges-utils")

vr2 = readMutect(mutect_path, strip = TRUE)

test_that("'granges' works", {
    gr2 = granges(vr2)

    expect_is( gr2, "GRanges" )
})

test_that("'ucsc/ncbi' works", {
    gr2 = granges(vr2)
    gr2_ncbi = ncbi(gr2)
    gr2_ucsc = ucsc(gr2_ncbi)

    expect_identical(seqchar(gr2), seqchar(gr2_ucsc))
})


context("mutationContextMutect")

test_that("'mutationContextMutect' works", {

    vr1 = readMutect(mutect_path)

    ct1 = mutationContextMutect(vr1)

    expect_is( ct1, "VRanges" )

    expect_error(mutationContextMutect(vr1, "notthere"))

})


context("mutationContext")

test_that("'mutationContext' works", {

    if(require(BSgenome.Hsapiens.1000genomes.hs37d5)) {

        vr0 = readMutect(mutect_path)
        vr1 = vr0
        mcols(vr1) = NULL

        vr2 = mutationContextMutect(vr0)
        vr1 = mutationContext(ncbi(vr1), BSgenome.Hsapiens.1000genomes.hs37d5)
        expect_equal(as.character(vr1$alteration), as.character(vr2$alteration))
        expect_equal(as.character(vr1$context), as.character(vr2$context.1))

        vr9 = vr1[1]
        alt(vr9) = "CA"
        expect_error(mutationContext(vr9, BSgenome.Hsapiens.1000genomes.hs37d5))

        vr9 = vr1[1]
        ref(vr9) = "CAT"
        expect_error(mutationContext(vr9, BSgenome.Hsapiens.1000genomes.hs37d5))

    }
})

context("Variant plots")

test_that("'plotRainfall' work", {

    plotRainfall(granges(vr))

    plotRainfall(vr)
    plotRainfall(vr, "alt")

})

test_that("'plotVariantAbundance' work", {

    plotVariantAbundance(vr)

})


context("datasets")

test_that("datasets for processing are available", {

    data("sca_mm", package = "SomaticSignatures")
    expect_is( sca_mm, "matrix" )
    expect_equal( nrow(sca_mm), 96 )

    data("sca_sigs", package = "SomaticSignatures")
    expect_is( sigs_nmf, "MutationalSignatures" )
    expect_is( sigs_pca, "MutationalSignatures" )    

})

context("identifySignatures")

test_that("'identifySignatures' works", {

    data("sca_mm", package = "SomaticSignatures")

    sigs_nmf = identifySignatures(sca_mm, 3)
    expect_is( sigs_nmf, "MutationalSignatures" )

    ## as many signatures as samples
    sigs_nmf = identifySignatures(sca_mm, ncol(sca_mm))
    expect_is( sigs_nmf, "MutationalSignatures" )
})


test_that("'identifySignatures' handles bad inputs", {

    data("sca_mm", package = "SomaticSignatures")

    ## bad number of signatures
    expect_error( identifySignatures(sca_mm, 0),
                 ".*single, positive integer.*")

    ## 'decomposition' argument no function
    expect_error( identifySignatures(sca_mm, 3, "a"),
                 "*.must be a function.*")

    ## 'x' no matrix
    expect_error( identifySignatures(1:10, 3) )
})


context("assessNumberSignatures")

test_that("'assessNumberSignatures' works", {

    data("sca_mm", package = "SomaticSignatures")

    gof_nmf = assessNumberSignatures(sca_mm, 2:4)
    expect_is( gof_nmf, "data.frame" )
    plotNumberSignatures(gof_nmf)

    gof_pca = assessNumberSignatures(sca_mm, 2:4, pcaDecomposition)
    expect_is( gof_pca, "data.frame" )
    plotNumberSignatures(gof_pca)

})

test_that("'assessNumberSignatures' low-level functions work", {

    x = matrix(rnorm(100, 20))
    y = matrix(rnorm(100, 20))

    library(NMF)

    r1 = rss(x, y)
    e1 = evar(x, y)
    expect_true( r1 > 0 )
    expect_true( e1 <= 1 )

    r1 = rss(x, x)
    e1 = evar(x, x)
    expect_equal( r1, 0 )
    expect_equal( e1, 1 )
})


context("plotting")

test_that("plotting functions works", {

    data("sca_sigs", package = "SomaticSignatures")

    plotObservedSpectrum(sigs_nmf)
    plotObservedSpectrum(sigs_nmf, "sample")
    plotObservedSpectrum(sigs_nmf, "alteration")
    expect_error( plotObservedSpectrum(sigs_nmf, "no") )

    plotFittedSpectrum(sigs_nmf)
    plotFittedSpectrum(sigs_nmf, "sample")
    plotFittedSpectrum(sigs_nmf, "alteration")
    expect_error( plotFittedSpectrum(sigs_nmf, "no") )

    plotSignatures(sigs_nmf)
    plotSignatures(sigs_nmf, normalize = TRUE)
    plotSignatures(sigs_nmf, normalize = TRUE, percent = TRUE)
    plotSignatures(sigs_nmf, normalize = FALSE, percent = TRUE)

    plotSignatureMap(sigs_nmf)

    plotSamples(sigs_nmf)
    plotSamples(sigs_nmf, normalize = TRUE)
    plotSamples(sigs_nmf, normalize = TRUE, percent = TRUE)
    plotSamples(sigs_nmf, normalize = FALSE, percent = TRUE)

    plotSampleMap(sigs_nmf)

})

context("Human chromosomes")

test_that("'hs*' work", {

    expect_equal( length(hsToplevel()), 25 )

    expect_equal( length(hsAutosomes()), 22 )

    expect_equal( length(hsAllosomes()), 2 )

    expect_equal( length(hsLinear()), 24 )

})
