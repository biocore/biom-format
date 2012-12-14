################################################################################
# clip out the first 3 characters, and
# name according to the taxonomic rank
#' @keywords internal
parseGreenGenesPrefix <- function(char.vec){
	# Define the meaning of each prefix according to GreenGenes taxonomy
	Tranks <- c(k="Kingdom", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species")
	taxvec        <- substr(char.vec, 4, 1000)
	names(taxvec) <- Tranks[substr(char.vec, 1, 1)]
	# Make sure order is same as Tranks
	taxvec <- taxvec[Tranks]; names(taxvec) <- Tranks
	return(taxvec)
}
################################################################################
################################################################################
#' Read a biom-format file, returning a \code{biom-class}.
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
#' @usage read_biom(biom_file, taxaPrefix=NULL, parallel=FALSE, version=1.0)
#'
#' @param biom_file (Required). A character string indicating the 
#'  file location of the biom formatted file. This is a JSON formatted file,
#'  specific to biological datasets, as described in 
#' 
#'  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
#' 
#' @param taxaPrefix (Optional). Character string. 
#'  What category of prefix precedes
#'  the taxonomic label at each taxonomic rank. 
#'  Currently only ``greengenes'' is a supported option, 
#'  and implies that the first letter indicates the taxonomic rank, 
#'  followed by two underscores and then the actual taxonomic
#'  assignment at that rank. The default value is \code{NULL}, meaning that
#'  no prefix or rank identifier will be interpreted. 
#'
#' @param parallel (Optional). Logical. Wrapper option for \code{.parallel}
#'  parameter in \code{plyr-package} functions. 
#'  If \code{TRUE},
#'  apply parsing functions in parallel using parallel backend
#'  declared with (probably) \code{foreach}-package and its supporting backend packages.
#'  This is implemented via plyr, so backend options could change depending on
#'  how plyr-package evolves.
#'  A second caveat,
#'  plyr-parallelization currently works most-cleanly with \code{multicore}-like
#'  backends (Mac OS X, Unix?), and may throw warnings for SNOW-like backends.
#'  See the example below for code invoking multicore-style backend within
#'  the \code{doParallel} package.
#' 
#' @param version (Optional). Numeric. The expected version number of the file.
#'  As the biom format evolves, version-specific importers will be available
#'  by adjusting the version value. Default is \code{1.0}. Not implemented.
#'  Has no effect (yet).
#'
#' @return An instance of the \code{biom-class}.
#'
#' @seealso 
#' The \code{\link{biom}} constructor function.
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
#'
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr ldply
#' @importFrom plyr laply
#' @export
#' @examples
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' read_biom(biom_file)
#' x <- read_biom(biom_file)
#' show(x)
#' print(x)
#' header(x)
#' abundance(x)
#' taxonomy(x)
#' sampleData(x)
#' tree(x, FALSE)
#' ## The previous example uses system.file() because of constraints in specifying a fixed
#' ##   path within a reproducible example in a package. 
#' ## In practice, however, you should simply provide "hard-link"
#' ## character string path to your file:
#' # mybiomfile <- "path/to/my/biomfile.biom"
#' # read_biom(mybiomfile)
read_biom <- function(biom_file, taxaPrefix=NULL, parallel=FALSE, version=1.0){
	
	# Read the data
	x <- fromJSON(biom_file)

	########################################
	# Header
	########################################
	# # Store header as everything up to the "rows" key.
	# # Protect against missing header keys causing an error
	header <- list() # Initialize header
	try(header <- x[1:(charmatch("rows", names(x)) - 1)], TRUE)
	if( identical(header, list()) ){
		warning("Header information missing from this file.")
	}
		
	########################################
	# OTU table:
	########################################
	# Check if sparse. Must parse differently than dense
	if( x$matrix_type == "sparse" ){
		otumat <- Matrix(0, nrow=x$shape[1], ncol=x$shape[2])
		# Loop through each sparse line and assign to relevant position in otumat.
		for( i in x$data ){
			otumat[(i[1]+1), (i[2]+1)] <- i[3]
		}
	} else if( x$matrix_type == "dense" ){ 
		# parse the dense matrix instead using plyr's laply
		otumat <- laply(x$data, function(i){i})
	}
	
	# Get row (OTU) and col (sample) names
	rownames(otumat) <- sapply(x$rows, function(i){i$id})
	colnames(otumat) <- sapply(x$columns, function(i){i$id})
	
	# Instantiates a "Matrix" daughter class, usually sparse,
	# but precise daughter class is chosen dynamically based
	# on properties of the data (if actually dense, it stays dense).
	ab.mat <- Matrix(otumat)
	
	########################################
	# Taxonomy Table
	########################################
	# Need to check if taxonomy information is empty (minimal biom file)
	if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
		taxdf <- NULL
	} else {
		# If GreenGenes, trim the prefix and use to name ultimately variable/column
		if( sum(taxaPrefix %in% "greengenes") > 0 ){
			taxdf <- ldply(x$rows, function(i){parseGreenGenesPrefix(i$metadata$taxonomy)}, .parallel=parallel)
		} else {
			taxdf <- ldply(x$rows, function(i){i$metadata$taxonomy}, .parallel=parallel)
		}
		# Add rownames to taxonomy data.frame.
		rownames(taxdf) <- sapply(x$rows, function(i){i$id})
	}
	
	########################################
	# Sample Data ("columns" in biom)
	########################################
	# If there is no metadata (all NULL), then set samdata to NULL, representing empty.
	if(  all( sapply(sapply(x$columns, function(i){i$metadata}), is.null) )  ){
		samdata <- NULL
	# Otherwise, parse it.
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
	# Use the biom() function to instantiate a biom-class with the data.
	########################################
	# Add header (not empty list by default, read file.)
	# Add tree read/test if that is implemented.
	return( biom(ab.mat, header=header, taxonomy=taxdf, sampleData=samdata, tree=NULL) )

}
################################################################################
################################################################################
################################################################################
################################################################################