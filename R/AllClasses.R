setClassUnion("TNInull", members = c("TNI", "NULL"))

#' MBR: an S4 class for co-regulation analysis and inference of 'dual regulons'.
#'
#' @slot TNI1 a 'TNI' object created by the RTN package.
#' @slot TNI2 another 'TNI' object created by the RTN package.
#' @slot testedElementsTNI1 regulatory elements listed in the TNI1.
#' @slot testedElementsTNI2 regulatory elements listed in the TNI2.
#' @slot dualsRegulons all possible 'duals regulons' computed by \code{\link[RTNmotifs:mbr.association]{mbr.association}}
#' @slot results a list, results from the MBR methods.
#' @slot para a list, parameters used in the MBR methods.
#' @slot summary a list, summary for 'para' and 'results'.
#' @slot status a character vector specifying the status of the MBR object based on the available methods.
#'
#' @exportClass MBR

##Class MBR (Motifs Between Regulons)
setClass(
    "MBR",
    slots=c(
        TNI1="TNInull",
        TNI2="TNInull",
        testedElementsTNI1="character",
        testedElementsTNI2="character",
        dualsRegulons="character",
        results="list",
        para='list',
        summary='list',
        status="character"
    ), prototype=list(
           TNI1=NULL,
           TNI2=NULL,
           testedElementsTNI1=character(),
           testedElementsTNI2=character(),
           dualsRegulons=character(),
           results=list(),
           para=list(),
           summary=list(),
           status=character()
    )
)
