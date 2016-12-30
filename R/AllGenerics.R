##generic functions
##------------------------------------------------------------------------------
setGeneric("mbr.preprocess",
           function(gexp, regulatoryElements1, regulatoryElements2, 
                    verbose=TRUE,...)
             standardGeneric("mbr.preprocess"), package="RTNmotifs")
setGeneric("mbr.permutation",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbr.permutation"), package="RTNmotifs")
setGeneric("mbr.bootstrap",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbr.bootstrap"), package="RTNmotifs")
setGeneric("mbr.dpi.filter",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbr.dpi.filter"), package="RTNmotifs")
setGeneric("mbr.association",
           function(object, regulatoryElements1=NULL, 
                    regulatoryElements2=NULL, minRegulonSize=50, 
                    prob=0.95, estimator='spearman', 
                    pAdjustMethod="BH", verbose=TRUE)
             standardGeneric("mbr.association"), package="RTNmotifs")
setGeneric("mbr.duals",
           function(object, supplementary.table = NULL,
                    evidenceColname=NULL, verbose = TRUE)
             standardGeneric("mbr.duals"), package="RTNmotifs")
setGeneric("tni2mbr.preprocess",
           function(tni1, tni2,  verbose = TRUE)
             standardGeneric("tni2mbr.preprocess"), package="RTNmotifs")
setGeneric("mbr_get",
           function(object, what="status")
             standardGeneric("mbr_get"), package="RTNmotifs")
