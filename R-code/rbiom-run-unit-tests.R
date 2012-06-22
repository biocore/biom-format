# Run manual tests
#
# First check that rbiom is available. 
# If not, attempt to install
if( !"rbiom" %in% .packages(TRUE) ){ # if( !"rbiom" %in% installed.packages()[, "Package"] ){
	system("R CMD INSTALL rbiom")	
}
# If still not available, throw error
if( !"rbiom" %in% .packages(TRUE) ){ 
	stop("Attempt to install rbiom-package prior to tests failed.\n
	 Please install package manually prior to testing.")
}
# # # # # # try(test <- require("rbiom", warn.conflicts=FALSE, quietly=TRUE), TRUE)

# Load rbiom package
library("rbiom")

# Testing package "testthat" is required for unit testing. It is in CRAN.
# Just install it in the usual way if it is missing.
if( !"rbiom" %in% .packages(TRUE) ){ 
	install.packages("testthat")
}
library("testthat")

# Now run the package unit tests.
test_package("rbiom")
