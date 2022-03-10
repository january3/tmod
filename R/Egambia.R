#' Transcriptomic responses to vaccination
#'
#' Transcriptomic responses to vaccination
#' 
#' The data shows the time course of transcriptomic responses to influenza
#' vaccination in healthy volunteers. The source of the data is 
#' GEO project PRJNA515032, associated with the following paper:
#' 
#' Weiner, January, et al. "Characterization of potential
#' biomarkers of reactogenicity of licensed antiviral vaccines: randomized
#' controlled clinical trials conducted by the BIOVACSAFE consortium."
#' Scientific reports 9.1 (2019): 1-14.
#'
#' For the data set, 3000 genes with top variance were chosen.
#' @format
#' Data frame with one row per gene containing log fold changes and FDR (q
#' values) for the Fluad vaccine as compared to placebo on day 0, day 1, day 2
#' and day 3 after the vaccination.
#' @name vaccination
NULL




#' Gene expression in TB patients and Healthy controls
#'
#' Gene expression in TB patients and Healthy controls
#'
#' This data set has been constructed from the gene expression data set
#' accessible in the Gene Expression Omnibus (GEO) at the accession number
#' GSE28623. Ten healthy donors (NID, non-infected donors) and 10 tubercolosis patients
#' (TB) have been randomly selected from the full data set, and top 25%
#' genes with the highest IQR have been selected for further analysis. Genes without
#' an Entrez gene (EG) identifier have likewise been omitted. 
#'
#' The Egambia object is a data frame. The first three columns are the gene
#' symbol, gene name and Entrez gene (EG) identifier. The remaining columns
#' correspond to the gene expression data.
#'
#' @examples
#' \dontrun{
#' # The data set has been generated as follows:
#' # get the data set from GEO
#' library(GEOquery)
#' gambia <- getGEO("GSE28623")[[1]]
#'
#' # Convert to limma and normalize
#' library(limma)
#' e <- new("EListRaw", list(E= exprs(gambia), genes=fData(gambia), targets= pData(gambia)))
#' e.bg <- backgroundCorrect(e, method= "normexp")
#' en <- normalizeBetweenArrays(e.bg, method= "q")
#' en <- avereps(en, ID=en$genes$NAME)
#' en <- en[ en$genes$CONTROL_TYPE == "FALSE", ]
#' en$targets$group <- factor(gsub(" whole blood RNA *", "", en$targets$description))
#'
#' # Fill in Entrez Gene IDs
#' library(org.Hs.eg.db)
#' en$genes$EG <- ""
#' sel <- en$genes$REFSEQ %in% ls(org.Hs.egREFSEQ2EG)
#' en$genes$EG[sel] <- mget(as.character(en$genes$REFSEQ[sel]), org.Hs.egREFSEQ2EG)
#'
#' # Filter by IQR and missing EG's
#' iqrs <- apply(en$E, 1, IQR)
#' en2 <- en[ iqrs > quantile(iqrs, 0.75) & en$genes$EG != "", ]
#'
#' # Select 10 random samples from NID and TB groups
#' en2 <- en2[ , c(sample(which(en2$targets$group == "NID"), 10), 
#'                  sample(which(en2$targets$group == "TB"), 10)) ]
#' colnames(en2$E) <- en2$targets$group
#' Egambia <- cbind(en2$genes[ , c("GENE_SYMBOL", "GENE_NAME", "EG") ], en2$E)
#' }
#' @name Egambia
NULL
