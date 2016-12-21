################################################################################
##########################         MBR-methods      ############################
################################################################################
##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "MBR",
          function(.Object, gexp, regulatoryElements1, regulatoryElements2)
           {
              #---checks
              if(missing(gexp)) stop("NOTE: 'gexp' is missing ", call.=FALSE)
              if(missing(regulatoryElements1)) stop("NOTE: 'regulatoryElements1' is missing", call.=FALSE)
              if(missing(regulatoryElements2)) stop("NOTE: 'regulatoryElements2' is missing", call.=FALSE)
              mbr.checks(name='gexp', gexp)
              mbr.checks(name='regulatoryElements1', regulatoryElements1)
              mbr.checks(name='regulatoryElements2', regulatoryElements2)

              #---creating TNIs
              regulonsTNI1 <- new("TNI", gexp=gexp, transcriptionFactors=regulatoryElements1)
              regulonsTNI2 <- new("TNI", gexp=gexp, transcriptionFactors=regulatoryElements2)
              .Object@TNI1 <- regulonsTNI1
              .Object@TNI2 <- regulonsTNI2

              #---status
              .Object@status <- rep('[ ]', 1, 5)
              names(.Object@status) <- c('Preprocess', 'Permutation', 'Bootstrap', 'DPI.filter', 'Association')

              # #---summary info
              # ##TNIs
              # sum.info <- list ()
              # sum.info$TNIs <- matrix (NA,2,2)
              # rownames (sum.info$TNIs) <- c('TNI1', 'TNI2')
              # colnames(sum.info$TNIs)<-c("targets", "edges")
              # .Object@summary <- sum.info

              ##parameters
              #association
              sum.info.para <- list()
              sum.info.para$TNIs$perm <- NA
              sum.info.para$TNIs$boot <- NA
              sum.info.para$TNIs$dpi <- NA
              sum.info.para$MBR$association <- matrix(NA, 1, 3)
              colnames(sum.info.para$MBR$association) <- c('minRegulonSize', 'prob', 'estimator')
              rownames(sum.info.para$MBR$association) <- 'Parameter'
              # #motifs
              # sum.info.para$motifs <-
              .Object@para <- sum.info.para
              ##summary
              #motifsInformation
              sum.info.summary <- list()
              sum.info.summary$MBR$Duals <- matrix(NA, 1, 3)
              colnames(sum.info.summary$MBR$Duals) <- c('numberDuals', 'topRValue', 'bottomRValue')
              rownames(sum.info.summary$MBR$Duals) <- 'duals'
              .Object@summary <- sum.info.summary

              return(.Object)
           }
          )

##------------------------------------------------------------------------------
#' A preprocessing function for objects of class MBR.
#'
#'Constructor for 'MBR-class' objects.
#'
#' @param gexp A numerical matrix, typically with mRNA and/or miRNA expression values.
#' @param regulatoryElements1 A named vector with regulatory elements listed in 'gexp' rownames.
#' @param regulatoryElements2 A named vector with regulatory elements listed in 'gexp' rownames.
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to \code{\link[RTN:tni.preprocess]{tni.preprocess}} function.
#' @return A preprocessed 'MBR-class' object.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#'
#' @importFrom RTN tni.preprocess
#'
#' @import methods
#' @docType methods
#' @rdname mbr.preprocess-methods
#' @aliases mbr.preprocess
#' @export

##Regulons pre-processing method
setMethod("mbr.preprocess",
           "matrix",
           function(gexp, regulatoryElements1, regulatoryElements2, verbose=TRUE,...)
               {
                   ##---
                   mbr.checks(name="verbose", para=verbose)  
                   object <- new("MBR", gexp=gexp, regulatoryElements1=regulatoryElements1, regulatoryElements2=regulatoryElements2)
                   ##---pre-processing TNIs
                   if(verbose) cat("-Preprocessing TNI objects...\n\n")
                   object@TNI1 <- tni.preprocess(object@TNI1, verbose=verbose, ...=...)
                   object@TNI2 <- tni.preprocess(object@TNI2, verbose=verbose,...=...)
                   object@status["Preprocess"] <- "[x]"
                   mbr.checks(name="regulatoryElements", para=c(object@TNI1@transcriptionFactors, object@TNI2@transcriptionFactors))
                   return(object)
               }
          )
##------------------------------------------------------------------------------
#' Inference of transcriptional networks.
#'
#' This function takes an MBR object and computes two transcriptional networks inferred 
#' by mutual information (with multiple hypothesis testing corrections).
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}.
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to the \code{\link[RTN:tni.permutation]{tni.permutation}} function.
#' @return An \linkS4class{MBR} object with two mutual information matrices, one in each "TNI" slot.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbr.permutation
#' rmbr <- mbr.permutation(rmbr, nPermutations=10)
#'
#' @importFrom RTN tni.permutation
#'
#' @import methods
#' @docType methods
#' @rdname mbr.permutation-methods
#' @aliases mbr.permutation
#' @export

## permutation
setMethod("mbr.permutation",
           "MBR",
           function(object, verbose=TRUE, ...)
               {
                   ##---checks
                   mbr.checks(name="object", para=object)
                   mbr.checks(name="verbose", para=verbose)
                   ##---permutation TNIs
                   if(verbose) cat("-Performing permutation analysis for two TNI objects...\n\n")
                   object@TNI1 <- tni.permutation(object@TNI1, verbose=verbose, ...=...)
                   object@TNI2 <- tni.permutation(object@TNI2, verbose=verbose, ...=...)
                   object@status["Permutation"] <- "[x]"
                   object@para$TNIs$perm <- object@TNI1@summary$para$perm
                   ##summaryTNIs
                   object@summary$TNIs$TNI1 <- object@TNI1@summary$results$tnet
                   colnames(object@summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
                   object@summary$TNIs$TNI2 <- object@TNI2@summary$results$tnet
                   colnames(object@summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
                   return(object)
               }
          )

#' Inference of consensus transcriptional networks.
#'
#' This function takes an MBR object and computes two consensus transcriptional networks.
#'
#' @param object A processed objec of class \linkS4class{MBR} evaluated by the method \code{\link[RTNmotifs:mbr.permutation]{mbr.permutation}}.
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} function.
#' @return An \linkS4class{MBR} object with two consensus mutual information matrices, one in each "TNI" slot.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbr.permutation
#' rmbr <- mbr.permutation(rmbr, nPermutations=10)
#' ##---mbr.bootstrap
#' rmbr <- mbr.bootstrap(rmbr, nBootstrap=10)
#'
#' @importFrom RTN tni.bootstrap
#'
#' @import methods
#' @docType methods
#' @rdname mbr.bootstrap-methods
#' @aliases mbr.bootstrap
#' @export

##------------------------------------------------------------------------------
## bootstrap method
setMethod("mbr.bootstrap",
           "MBR",
           function(object, verbose=TRUE, ...)
               {
                   ##---checks
                   mbr.checks(name="object", para=object)
                   mbr.checks(name="verbose", para=verbose)
                   ##---bootstrap TNIs
                   if(verbose) cat("-Performing bootstrap analysis for two TNI objects...\n\n")
                   object@TNI1 <- tni.bootstrap(object@TNI1, verbose=verbose, ...=...)
                   object@TNI2 <- tni.bootstrap(object@TNI2, verbose=verbose, ...=...)
                   object@status["Bootstrap"] <- "[x]"
                   object@para$TNIs$boot <- object@TNI1@summary$para$boot
                   ##summaryTNIs
                   object@summary$TNIs$TNI1 <- object@TNI1@summary$results$tnet
                   colnames(object@summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
                   object@summary$TNIs$TNI2 <- object@TNI2@summary$results$tnet
                   colnames(object@summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
                   return(object)
               }
          )

#' A filter based on the Data Processing Inequality (DPI) algorithm.
#'
#' This function takes an MBR object and computes two transcriptional networks filtered by the data processing inequality algorithm.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the methods
#'  \code{\link[RTNmotifs:mbr.permutation]{mbr.permutation}} and \code{\link[RTNmotifs:mbr.bootstrap]{mbr.bootstrap}}.
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} function.
#' @return An \linkS4class{MBR} object with two DPI-filtered mutual information matrices, one in each "TNI" slot.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbr.permutation
#' rmbr <- mbr.permutation(rmbr, nPermutations=10)
#' ##---mbr.bootstrap
#' rmbr <- mbr.bootstrap(rmbr, nBootstrap=10)
#' ##---mbr.dpi.filter
#' rmbr <- mbr.dpi.filter(rmbr)
#'
#' @importFrom RTN tni.dpi.filter
#'
#' @import methods
#' @docType methods
#' @rdname mbr.dpi.filter-methods
#' @aliases mbr.dpi.filter
#' @export

##------------------------------------------------------------------------------
## dpi filter method
setMethod("mbr.dpi.filter",
           "MBR",
           function(object, verbose=TRUE, ...)
               {
                   ##---checks
                   mbr.checks(name="object", para=object)
                   mbr.checks(name="verbose", para=verbose)
                   ##---Dpi filter TNIs
                   if(verbose) cat("-Applying dpi filter for two TNI objects...\n")
                   object@TNI1 <- tni.dpi.filter(object@TNI1, verbose=verbose, ...=...)
                   object@TNI2 <- tni.dpi.filter(object@TNI2, verbose=verbose, ...=...)
                   object@status["DPI.filter"] <- "[x]"
                   object@para$TNIs$dpi <- object@TNI1@summary$para$dpi
                   ##summary TNIs
                   object@summary$TNIs$TNI1 <- object@TNI1@summary$results$tnet
                   colnames(object@summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
                   object@summary$TNIs$TNI2 <- object@TNI2@summary$results$tnet
                   colnames(object@summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
                   return(object)
               }
          )

#' Motifs analysis and inference of 'dual regulons'.
#'
#' This function takes an MBR object and compares the shared regulon 
#' targets in order to test whether regulon pairs agree on the predicted downstream effects.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' methods \code{\link[RTNmotifs:mbr.permutation]{mbr.permutation}}, \code{\link[RTNmotifs:mbr.bootstrap]{mbr.bootstrap}} 
#' and \code{\link[RTNmotifs:mbr.dpi.filter]{mbr.dpi.filter}}.
#' @param regulatoryElements1 An optional character vector specifying which 'TNI1' regulatory elements should be evaluated. If 'NULL' all regulatory elements will be evaluated.
#' @param regulatoryElements2 An optional character vector specifying which 'TNI2' regulatory elements should be evaluated. If 'NULL' all regulatory elements will be evaluated.
#' @param minRegulonSize A single integer or numeric value specifying the minimum number of elements in a regulon. Gene sets with fewer than this number are removed from the analysis.
#' @param prob A quantile filter applyed to the association metric used to infer 'dual regulons'.
#' @param estimator A character value specifying the estimator used in the association analysis. One of "spearman" (default), "kendall", or "pearson", can be abbreviated.
#' @param pAdjustMethod A single character value specifying the p-value adjustment method to be used (see 'p.adjust' for details).
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with two data.frames in the slot 'results' listing the inferred 'dual regulons' and a hypergeometric test for each 'dual regulon'.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbr.permutation
#' rmbr <- mbr.permutation(rmbr, nPermutations=10)
#' ##---mbr.bootstrap
#' rmbr <- mbr.bootstrap(rmbr, nBootstrap=10)
#' ##---mbr.dpi.filter
#' rmbr <- mbr.dpi.filter(rmbr)
#' ##---mbr.association
#' rmbr <- mbr.association(rmbr, prob=0.75)
#'
#' @importFrom RTN tni.get
#' @importFrom stats p.adjust phyper
#' @importFrom stats cor quantile
#'
#' @import methods
#' @docType methods
#' @rdname mbr.association-methods
#' @aliases mbr.association
#' @export

##------------------------------------------------------------------------------
##Inference of duals
setMethod("mbr.association",
          "MBR",
          function(object, regulatoryElements1=NULL, regulatoryElements2=NULL, minRegulonSize=30, prob=0.95, 
                   estimator='spearman', pAdjustMethod="BH", verbose=TRUE)
          {
            ##---checks
            mbr.checks(name="minRegulonSize", para=minRegulonSize)
            mbr.checks(name="prob", para=prob)
            mbr.checks(name="estimator", para=estimator)
            mbr.checks(name="pAdjustMethod", para=pAdjustMethod)
            mbr.checks(name="verbose", para=verbose)
            if(is.null(regulatoryElements1))
            {
              if(verbose) cat("-Selecting regulatory elements from TNI1 object...\n")
              regulatoryElements1 <- object@TNI1@transcriptionFactors
            } else
            {
              regulatoryElements1 <- .checkRegel(object@TNI1, regulatoryElements1)
            }
            if(is.null(regulatoryElements2))
            {
              if(verbose) cat("-Selecting regulatory elements from TNI2 object...\n")
              regulatoryElements2 <- object@TNI2@transcriptionFactors
            } else
            {
              regulatoryElements2 <- .checkRegel(object@TNI2, regulatoryElements2)
            }
            mbr.checks(name="regulatoryElements", para=c(regulatoryElements1, regulatoryElements2))
            mbr.checks(name="numberRegElements", para=c(regulatoryElements1, regulatoryElements2))
            
            ##-----get regulons
            what <- "refregulons.and.mode"
            ##if (tnet=="dpi") what <- "regulons.and.mode"
            regulons1 <- tni.get(object@TNI1, what=what)
            regulons2 <- tni.get(object@TNI2, what=what)
            
            ##-----get regulatory elements
            regulons1 <- regulons1[regulatoryElements1]
            regulons2 <- regulons2[regulatoryElements2]
            
            ##-----get regulons by min size
            size1 <- unlist(lapply(regulons1, length))
            idx <- size1 >= (minRegulonSize)
            regulons1 <- regulons1[idx]; regulatoryElements1 <- regulatoryElements1[idx]
            size1 <- size1[idx]
            ##---
            size2 <- unlist(lapply(regulons2, length))
            idx <- size2 >= (minRegulonSize)
            regulons2 <- regulons2[idx]; regulatoryElements2 <- regulatoryElements2[idx]
            size2 <- size2[idx]
            
            ##-----group regulons and regulatory elements
            regulons <- c(regulons1, regulons2)
            regel <- c(regulatoryElements1, regulatoryElements2)
            
            ##-----Correlation
            tnet <- .regMatrix(regulons, regel, getNames=FALSE)
            tbmi <- tnet
            ##-----
            tnet <- .tni.cor(object@TNI1@gexp, tnet, estimator=estimator, dg=0, asInteger=FALSE, mapAssignedAssociation=TRUE)
            regcor <- cor(tnet[, regulatoryElements1], tnet[, regulatoryElements2], method=estimator)
            ##-----
            rownames(regcor) <- names(regulatoryElements1)
            colnames(regcor) <- names(regulatoryElements2)
            
            
            ##-----select motifs based on 'prob' quantile
            if(verbose) cat(paste("-Getting duals in probs >", prob, "...\n", sep = ""))
            pvlist <- .motifsquantile(regcor=regcor, th=prob)
            
            if(nrow(pvlist) > 0)
            {
              ##-----Mutual Information
              pvlist <- .getMIorP(pvlist, regel, tbmi, size1, size2)
              ##----PadjustValue
              if(!is.null(object@TNI1@results$adjpv))
              {
                tbpValue <- object@TNI1@results$adjpv
                cutoff <- object@TNI1@para$perm$pValueCutoff
                pvlist <- .getMIorP(pvlist, regel, tbpValue, size1, size2, mutualInformation=FALSE, cutoff=cutoff)
              }
              ##-----Jaccard
              pvlist <- .jcOverlap(pvlist, regel, tbmi)
              
              ##-----Hypergeometric
              tni1 <- mbr.get(object, what="TNI1")
              universe <- rownames(tni1@gexp)
              hyperresults <- .mbr.hyper(pvlist=pvlist, regulons=regulons, regel=regel, universe=universe, pAdjustMethod=pAdjustMethod, verbose=verbose)
              pvlist$Hypergeometric.Adjusted.Pvalue <- hyperresults$Adjusted.Pvalue
              pvlist$Hypergeometric.Pvalue <- hyperresults$Pvalue
            }else
            {
              warning("No 'dual' has been found for the input parameters.")
            }
            ##-----organize pvlist
            pvlist <- pvlist[,c("Regulon1","Size.Regulon1","Regulon2","Size.Regulon2","Jaccard.coefficient", "Hypergeometric.Pvalue","Hypergeometric.Adjusted.Pvalue", "MI","MI.Adjusted.Pvalue","R","Quantile")]
            ##-----
            sum.info.par <- c(minRegulonSize, prob, estimator)
            object@para$MBR$association['Parameter', ] <- sum.info.par
            object@status ["Association"] <- "[x]"
            ##---
            info.summary.results <- c(nrow(pvlist), max(pvlist$R), min(pvlist$R))
            object@summary$MBR$Duals['duals',] <- info.summary.results
            ##-----
            object@testedElementsTNI1 <- regulatoryElements1
            object@testedElementsTNI2 <- regulatoryElements2
            object@dualsRegulons <- rownames(pvlist)
            object@results$motifsInformation <- pvlist
            object@results$hypergeometricResults <- hyperresults
            return(object)
          }
)


#' A summary for results from the MBR methods.
#'
#' This function lists the inferred 'dual regulons' and, if available, adds external evidences.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the method \code{\link[RTNmotifs:mbr.association]{mbr.association}}.
#' @param supplementary.table An optional 'data.frame' with three columns representing 
#' (1) regulatory elements of 'TNI1', (2) regulatory elements of 'TNI2', and 
#' (3) external evidences between the regulatory elements.
#' @param evidenceColname A single character value specifying a column in the 'supplementary.table'.
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with a data.frame in the slot 'results' listing the input additional evidences.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbr.permutation
#' rmbr <- mbr.permutation(rmbr, nPermutations=10)
#' ##---mbr.bootstrap
#' rmbr <- mbr.bootstrap(rmbr, nBootstrap=10)
#' ##---mbr.dpi.filter
#' rmbr <- mbr.dpi.filter(rmbr)
#' ##---mbr.association
#' rmbr <- mbr.association(rmbr, prob=0.75)
#' ##---a 'toy' table with supplementary evidences
#' motifsInformation <- mbr.get(rmbr, what="motifsInformation")
#' n <- nrow(motifsInformation)
#' idx <- sample(1:n, size=round(n/2))
#' supplementaryTable <- motifsInformation[idx,c("Regulon1","Regulon2")]
#' supplementaryTable$ToyEvidence <- rnorm(nrow(supplementaryTable))
#' ##---mbr.duals
#' rmbr <- mbr.duals(rmbr, supplementary.table = supplementaryTable, evidenceColname = "ToyEvidence")
#' ##---motifsInformation with 'Evidence'
#' motifsInformation <- mbr.get(rmbr, what="motifsInformation")
#' head(motifsInformation)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import methods
#' @docType methods
#' @rdname mbr.duals-methods
#' @aliases mbr.duals
#' @export

##------------------------------------------------------------------------------
##organize duals
setMethod( "mbr.duals",
            "MBR",
          function(object, supplementary.table=NULL, evidenceColname, verbose=TRUE)
              {
                  ##---checks
                  mbr.checks(name="object", para=object)
                  ##---
                  motifsInformation <- object@results$motifsInformation
                  if(is.null(dim(motifsInformation)))
                      stop("'motifsInformation' seems null!")
                  if(verbose) cat("-Sorting by the R value...\n")
                  idx <- sort(abs(motifsInformation[,"R"]), decreasing=TRUE, index.return=TRUE)
                  object@results$motifsInformation <- motifsInformation[idx$ix, ]
                  object@dualsRegulons <- rownames(object@results$motifsInformation)
                  if(!is.null(supplementary.table))
                      {
                          ##---checks
                          if(missing(evidenceColname)) stop("'evidenceColname' should be a character value present in colnames of supplementary.table!")
                          mbr.checks(name="supplementary.table", para=supplementary.table)
                          mbr.checks(name="uniqueInput", para=supplementary.table)
                          mbr.checks(name="evidenceColname", para=evidenceColname)
                          
                          ##---consistency
                          if(verbose) cat("-Checking the 'supplementary.table' consistency...\n")
                          supplementary.table <- .consisSuppTable(object, supplementary.table, evidenceColname, verbose=verbose)
                          ##---find duals
                          object <- .checkLoops(object, supplementary.table, evidenceColname, verbose=verbose)
                  }
                  return(object)
              }
          )

#' A preprocessing function for objects of class MBR.
#'
#' This function merges two TNI class objects and creates one MBR class object.
#'
#' @param tni1 A 'TNI' class object.
#' @param tni2 Another 'TNI' class object
#' @param verbose A single logical value specifying to display detailed messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbr.permutation
#' rmbr <- mbr.permutation(rmbr, nPermutations=10)
#' ##---mbr.bootstrap
#' rmbr <- mbr.bootstrap(rmbr, nBootstrap=10)
#' rmbr <- mbr.dpi.filter(rmbr)
#' ##---tni2mbr.preprocess
#' tni1 <- mbr.get(rmbr, what="TNI1")
#' tni2 <- mbr.get(rmbr, what="TNI2")
#' rmbr <- tni2mbr.preprocess(tni1, tni2)
#'
#' @import methods
#' @docType methods
#' @rdname tni2mbr.preprocess-methods
#' @aliases tni2mbr.preprocess
#' @export

##----------------------------------------------------------------
##Combine two TNIs produced separately
setMethod("tni2mbr.preprocess",
          "TNI",
          function (tni1,  tni2,  verbose=TRUE)
              {
                  mbr.checks (name='tni', para=tni1)
                  mbr.checks (name='tni', para=tni2)
                  .combineTNIs (tni1=tni1, tni2=tni2, verbose=verbose)
                  ##---- creates MBR object
                  object <- new("MBR", gexp=tni1@gexp,
                                regulatoryElements1=tni1@transcriptionFactors,
                                regulatoryElements2=tni2@transcriptionFactors)
                  object@TNI1 <- tni1
                  object@TNI2 <- tni2
                  status <- names(object@TNI1@status[object@TNI1@status=="[x]"])
                  object@status[status] <- "[x]"
                  ##---permutation
                  object@para$TNIs$perm <- tni1@summary$para$perm
                  object@summary$TNIs$TNI1 <- tni1@summary$results$tnet
                  colnames(object@summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
                  object@summary$TNIs$TNI2 <- tni2@summary$results$tnet
                  colnames(object@summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
                  ##---bootstrap
                  object@para$TNIs$boot <- tni1@summary$para$boot
                  object@summary$TNIs$TNI1 <- tni1@summary$results$tnet
                  colnames(object@summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
                  object@summary$TNIs$TNI2 <- tni2@summary$results$tnet
                  colnames(object@summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
                  ##---summary dpi.filter
                  object@para$TNIs$dpi <- tni1@summary$para$dpi
                  object@summary$TNIs$TNI1 <- tni1@summary$results$tnet
                  colnames(object@summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
                  object@summary$TNIs$TNI2 <- tni2@summary$results$tnet
                  colnames(object@summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
                  #object
                  return (object)
              }
          )

##----------------------------------------------------------------
##show summary information on screen
setMethod( "show",
            "MBR",
          function(object)
          {
            cat("an MBR (Motifs Between Regulons) object:\n")
            message("--status:")
            print(object@status, quote=FALSE)
          }
         )

#' Get information from individual slots in MBR object.
#' 
#' Get information from individual slots in an MBR object and any available results from previous analysis.
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}
#' @param what a single character value specifying which information should be retrieved from the slots. Options: "TNI1", "TNI2", "testedElementsTNI1", "testedElementsTNI2", "dualsRegulons", "results", "para", "summary", "status", "motifsInformation" and "hyperResults"
#' @return A slot content from a object of class 'MBR' \linkS4class{MBR} object
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbr.preprocess
#' rmbr <- mbr.preprocess(gexp=gexp, regulatoryElements1 = tfs1, regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---get the 'TNI1' slot using 'mbr.get'
#' tni1 <- mbr.get(rmbr, what="TNI1")
#' 
#' @import methods
#' @docType methods
#' @rdname mbr.get-methods
#' @aliases mbr.get
#' @export
##----------------------------------------------------------------
##get slots from MBR object
setMethod( "mbr.get",
           "MBR", 
           function(object, what="status")
           {
             ##---check input arguments
             mbr.checks(name="object", para=object)
             mbr.checks(name="mbr.get", para=what)
             ##---Association options any change needs update!
             optsAssoci <- c("testedElementsTNI1", "testedElementsTNI2", "dualsRegulons", "motifsInformation", "results", "hyperResults")
             ##---get query
             if(what=="TNI1")
             {
               query <- object@TNI1
             }
             else if(what=="TNI2")
             {
               query <- object@TNI2
             }
             else if(what=="para")
             {
               query <- object@para
             }
             else if(what=="summary")
             {
               query <- object@summary
             }
             else if(what=="status")
             {
               query <- object@status
             }
             else if(what%in%optsAssoci)
             {
               if(object@status["Association"] != "[x]")
               {
                 warning("NOTE: input 'object' needs 'mbr.association' evaluation!")
                 query <- NULL
               } else {
                 if(what=="testedElementsTNI1")
                 {
                   query <- object@testedElementsTNI1
                 }
                 else if(what=="testedElementsTNI2")
                 {
                   query <- object@testedElementsTNI2
                 }
                 else if(what=="dualsRegulons")
                 {
                   query <- object@dualsRegulons
                 }
                 else if(what=="results")
                 {
                   query <- object@results
                 }
                 else if(what=="motifsInformation")
                 {
                   query <- object@results$motifsInformation
                 }
                 else if(what=="hyperResults")
                 {
                    query <- object@results$hypergeometricResults
                 }
               }
             }
             return(query)
           }
          )
