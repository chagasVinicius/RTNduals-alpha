##generic functions
##------------------------------------------------------------------------------
setGeneric("mbrPreprocess",
           function(gexp, regulatoryElements1, regulatoryElements2, 
                    verbose=TRUE,...)
             standardGeneric("mbrPreprocess"), package="RTNmotifs")
setGeneric("mbrPermutation",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrPermutation"), package="RTNmotifs")
setGeneric("mbrBootstrap",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrBootstrap"), package="RTNmotifs")
setGeneric("mbrDpiFilter",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrDpiFilter"), package="RTNmotifs")
setGeneric("mbrAssociation",
           function(object, regulatoryElements1=NULL, 
                    regulatoryElements2=NULL, minRegulonSize=50, 
                    prob=0.95, estimator='spearman', 
                    pAdjustMethod="BH", verbose=TRUE)
             standardGeneric("mbrAssociation"), package="RTNmotifs")
setGeneric("mbrDuals",
           function(object, supplementary.table = NULL,
                    evidenceColname=NULL, verbose = TRUE)
             standardGeneric("mbrDuals"), package="RTNmotifs")
setGeneric("tni2mbrPreprocess",
           function(tni1, tni2,  verbose = TRUE)
             standardGeneric("tni2mbrPreprocess"), package="RTNmotifs")
setGeneric("mbrGet",
           function(object, what="status")
             standardGeneric("mbrGet"), package="RTNmotifs")
