################################################################################
# clip out the first 3 characters, and
# name according to the taxonomic rank
#' @keywords internal
parseGreenGenesPrefix <- function(char.vec){
	# Define the meaning of each prefix according to GreenGenes (and RDP?) taxonomy
	Tranks <- c(k="Kingdom", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species")
	taxvec        <- substr(char.vec, 4, 1000)
	names(taxvec) <- Tranks[substr(char.vec, 1, 1)]
	# Make sure order is same as Tranks
	taxvec <- taxvec[Tranks]; names(taxvec) <- Tranks
	return(taxvec)
}
################################################################################
################################################################################
#' Read a BIOM file
#'
#' New versions of QIIME produce a more-comprehensive and formally-defined
#' JSON file format. From the QIIME website:
#'
#' ``The biom file format (canonically pronounced `biome') is designed to be a 
#' general-use format for representing counts of observations in one or
#' more biological samples.''
#' 
#' \url{http://biom-format.org/} 
#'
#' @usage read_biom(BIOMfilename, taxaPrefix=NULL, version=0.9)
#'
#' @param BIOMfilename (Required). A character string indicating the 
#'  file location of the BIOM formatted file. This is a JSON formatted file,
#'  specific to biological datasets, as described in 
#' 
#'  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
#' 
#' @param taxaPrefix (Optional). Character string. What category of prefix precedes
#'  the taxonomic label at each taxonomic rank. Currently only ``greengenes'' is
#'  a supported option, and implies that the first letter indicates the 
#'  taxonomic rank, followed by two underscores and then the actual taxonomic
#'  assignment at that rank. The default value is \code{NULL}, meaning that
#'  no prefix or rank identifier will be interpreted. 
#' 
#' @param version (Optional). Numeric. The expected version number of the file.
#'  As the BIOM format evolves, version-specific importers will be available
#'  by adjusting the version value. Default is \code{1.0}. Not implemented.
#'  Has no effect (yet).
#'
#' @return An instance of the \code{biom-class}.
#'
#' @references \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
#'
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr ldply
#' @importFrom plyr laply
#' @export
#' @examples
#'  # # # import with default parameters, specify a file
#'  # read_biom(myBIOMfile)
read_biom <- function(BIOMfilename, taxaPrefix=NULL, version=0.9){
	
	# Read the data
	x <- fromJSON(BIOMfilename)
	
	########################################
	# OTU table:
	########################################
	# Check if sparse. Must parse differently than dense
	if( x$matrix_type == "sparse" ){
		otumat <- matrix(0, nrow=x$shape[1], ncol=x$shape[2])
		dummy <- sapply(x$data, function(i){otumat[(i[1]+1), (i[2]+1)] <<- i[3]})
	}
	# parse the dense matrix instead.
	if( x$matrix_type == "dense" ){
		# each row will be complete data values, should use laply
		otumat <- laply(x$data, function(i){i})
	}
	
	# Get row (OTU) and col (sample) names
	rownames(otumat) <- sapply(x$rows, function(i){i$id})
	colnames(otumat) <- sapply(x$columns, function(i){i$id})
	
	# Instantiates a "Matrix" class object, usually sparse.
	ab.mat <- Matrix(otumat)
	
	########################################
	# Taxonomy Table
	########################################
	# Need to check if taxonomy information is empty (minimal BIOM file)
	if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
		taxtab <- NULL
	} else {
		# taxdf <- laply(x$rows, function(i){i$metadata$taxonomy}, .parallel=parallel)
		# Figure out the max number of columns (could be jagged in BIOM format)
		ncols <- max(sapply(head(x$rows), function(i){length(i$metadata$taxonomy)}))
		# Initialize character matrix
		taxdf <- matrix(NA_character_, nrow=length(x$rows), ncol=ncols)
		# Fill in the matrix by row.
		for( i in 1:length(x$rows) ){
			if( sum(taxaPrefix %in% "greengenes") > 0 ){
				taxdf[i, 1:length(x$rows[[i]]$metadata$taxonomy)] <- parseGreenGenesPrefix(x$rows[[i]]$metadata$taxonomy)
			} else {
				taxdf[i, 1:length(x$rows[[i]]$metadata$taxonomy)] <- x$rows[[i]]$metadata$taxonomy
			}
		}
		# Now convert to matrix, name the rows as "id" (the taxa name), coerce to taxonomyTable
		# Note that there is no character matrix in the Matrix-package. standard R char matrix used instead.
		taxtab           <- as(taxdf, "matrix")
		rownames(taxtab) <- sapply(x$rows, function(i){i$id})
	}
	
	########################################
	# Sample Data ("columns" in QIIME/BIOM)
	########################################
	# If there is no metadata (all NULL), then set samdata <- NULL
	if(  all( sapply(sapply(x$columns, function(i){i$metadata}), is.null) )  ){
		samdata <- NULL
	} else {
		samdata <- ldply(x$columns, function(i){
			if( class(i$metadata) == "list"){
				return(i$metadata[[1]])
			} else {
				return(i$metadata)				
			}
		}, .parallel=parallel)
		rownames(samdata) <- sapply(x$columns, function(i){i$id})

	}
	
	########################################
	# Put together into a list. 
	# Should eventually be a BIOM object (special class)
	########################################
	# Add header (not empty list by default, read file.)
	# Add tree read/test if that is implemented.
	return( BIOM(ab.mat, header=list(), taxonomy=taxtab, sampleData=samdata, tree=NULL) )

}
################################################################################
################################################################################
################################################################################
################################################################################