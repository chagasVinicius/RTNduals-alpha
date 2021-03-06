\name{RTNduals-package}
\alias{RTNduals-package}
\alias{RTNduals}
\docType{package}

\title{
RTNduals: An R/Bioconductor package for analysis of co-regulatory network motifs and inference of 'dual regulons'.
}

\description{
RTNduals is a tool that searches for possible co-regulatory loops between regulon pairs generated by the RTN package. It compares the shared targets in order to infer "dual regulons", a new concept that tests whether regulon pairs agree on the predicted downstream effects.
}

\details{
\tabular{ll}{
Package: \tab RTNduals\cr
Type: \tab Package\cr
Depends: \tab R (>= 3.4.0), methods, RTN\cr
Imports: \tab grDevices, stats, utils\cr
Suggests: \tab knitr, rmarkdown, BiocStyle, RUnit, BiocGenerics\cr
License: \tab Artistic-2.0\cr
biocViews: \tab NetworkInference, NetworkEnrichment, GeneRegulation, 
GeneExpression, GraphAndNetwork\cr
}
}

\section{Index}{

\tabular{ll}{
\link[RTNduals:MBR-class]{MBR-class}: \tab an S4 class for co-regulation 
analysis and inference of 'dual regulons'.\cr
\link[RTNduals:mbrPreprocess]{mbrPreprocess}: \tab a preprocessing function 
for objects of class MBR.\cr
\link[RTNduals:mbrPermutation]{mbrPermutation}: \tab inference of 
transcriptional networks.\cr
\link[RTNduals:mbrBootstrap]{mbrBootstrap}: \tab inference of consensus 
transcriptional networks.\cr
\link[RTNduals:mbrDpiFilter]{mbrDpiFilter}: \tab a filter based on the 
Data Processing Inequality (DPI) algorithm.\cr
\link[RTNduals:mbrAssociation]{mbrAssociation}: \tab motifs analysis and 
inference of "dual regulons".\cr
\link[RTNduals:mbrDuals]{mbrDuals}: \tab a summary for results from the 
MBR methods.\cr
\link[RTNduals:tni2mbrPreprocess]{tni2mbrPreprocess}: \tab a preprocessing 
function for objects of class MBR.\cr
}

Further information is available in the vignettes by typing 
\code{vignette("RTNduals")}. Documented 
topics are also available in HTML by typing \code{help.start()} and 
selecting the RTNduals package from the menu.

}

\author{
Vinicius S. Chagas, Clarice S. Groeneveld, Kerstin B Meyer, Gordon Robertson, Mauro A. A. Castro
}

\references{
Fletcher M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast 
cancer risk.} Nature Communications, 4:2464, 2013.

Castro M.A.A. et al., \emph{Regulators of genetic risk of breast cancer 
identified by integrative network analysis.} Nature Genetics, 48:12-21, 2016.
}

\keyword{ package }
