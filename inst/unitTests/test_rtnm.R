# Unit tests for MBR-class methods
test_mbr <- function()
{
    data("dt4rtn", package = "RTN")
    tfs1 <- dt4rtn$tfs[c("FOXM1", "E2F2")]
    tfs2 <- dt4rtn$tfs[c("PTTG1", "RARA")]
    ##mbr.preprocess
    rmbr <- mbr.preprocess(gexp=dt4rtn$gexp, regulatoryElements1=tfs1, regulatoryElements2=tfs2, gexpIDs=dt4rtn$gexpIDs)
    status <- mbr.get(rmbr, what="status")
    checkTrue(status["Preprocess"]=="[x]" && status[1]=="[x]" && status[1]=="[x]")
    ##mbr.permutation
    rmbr <- mbr.permutation(rmbr, nPermutations=10, estimator="pearson")
    status <- mbr.get(rmbr, what="status")
    checkTrue(status["Permutation"]=="[x]" && status[2]=="[x]" && status[2]=="[x]")
    ##mbr.bootstrap
    rmbr <- mbr.bootstrap(rmbr, estimator="pearson", nBootstrap=10)
    status <- mbr.get(rmbr, what="status")
    checkTrue(status["Bootstrap"]=="[x]" && status[3]=="[x]" && status[3]=="[x]")
    ##mbr.dpi.filter
    rmbr <- mbr.dpi.filter(rmbr)
    status <- mbr.get(rmbr, what="status")
    checkTrue(status["DPI.filter"]=="[x]" && status[4]=="[x]" && status[4]=="[x]")
    ##mbr.combine.TNIs
    tni1 <- mbr.get(rmbr, what="TNI1")
    tni2 <- mbr.get(rmbr, what="TNI2")
    rmbr <- tni2mbr.preprocess(tni1, tni2)
    status <- mbr.get(rmbr, what="status")
    checkTrue(status[1:4]=="[x]" && status[1:4]=="[x]" && status[1:4]=="[x]")
    ##mbr.association
    rmbr <- mbr.association(rmbr, prob=0, estimator="pearson")
    status <- mbr.get(rmbr, what="status")
    motifsInformation <- mbr.get(rmbr, what="motifsInformation")
    checkTrue(status["Association"]=="[x]" && is.data.frame(motifsInformation) && ncol(motifsInformation)==11)
    ##mbr.motifs
    rmbr <- mbr.duals(rmbr)
    motifsInformation <- mbr.get(rmbr, what="motifsInformation")
    checkTrue(is.data.frame(motifsInformation))
}
