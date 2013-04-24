# Run manual tests
#
# First check that biom is available. 
# If not, attempt to install
if( !"biom" %in% .packages(TRUE) ){
	# if( !"biom" %in% installed.packages()[, "Package"] ){
	system("R CMD INSTALL biom")	
}
# If still not available, throw error
if( !"biom" %in% .packages(TRUE) ){ 
	stop("Attempt to install biom-package prior to tests failed.\n
	 Please install package manually prior to testing.")
}

# Load biom package
library("biom")

# Testing package "testthat" is required for unit testing. It is in CRAN.
# Just install it in the usual way if it is missing.
if( !"testthat" %in% .packages(TRUE) ){ 
	install.packages("testthat")
}
library("testthat")

# Now run the package unit tests.
test_package("biom")
