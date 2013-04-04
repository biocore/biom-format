################################################################################
# Validity methods:
#
# These are delicate, because they are effectively at the S4 infrastructure
# level, in between "new" and the constructor. Some of the issues that might
# otherwise go here can be also be dealt with in the constructors, if needed.
################################################################################
validbiom = function(object){
	
	# All the required top-level keys are present
	biomreqdkeys = c("id", "format", "format_url", "type", "generated_by", "date",
									 "rows", "columns", "matrix_type", "matrix_element_type", "shape", "data")
	if( !all(biomreqdkeys %in% names(object)) ){
		missingkeysmessage = paste("", 
					"Not all required top-level keys are present in biom-object.",
					"Required keys are:", paste0(biomreqdkeys, collapse="\n"),
					sep="\n")
		return(missingkeysmessage)
	}
	
	# Matrix shape and number of row/col elements matches
	bshape = biomshape(object)
	if( bshape["ncol"] > length(object$columns) ){
		return("shape specifies more cols than are present in metadata")
	}
	if( bshape["nrow"] > length(object$rows) ){
		return("shape field specifies more rows than are present in metadata")
	}
	if( bshape["ncol"] < length(object$columns) ){
		return("more metadata cols than specified by shape field")
	}
	if( bshape["nrow"] < length(object$rows) ){
		return("more metadata rows than specified by shape field")
	}	
	
	# The type field has an acceptable value
	biomtypevals = c("OTU table",
		"Pathway table",
		"Function table",
		"Ortholog table",
		"Gene table",
		"Metabolite table",
		"Taxon table"
	)
	if( !object$type %in% biomtypevals ){
		return("type field has unsupported value")
	}
	
	# The matrix_type has an acceptable value
	biommatrixtypes = c("sparse", "dense")
	if( !object$matrix_type %in% biommatrixtypes ){
		return("matrix_type field has unsupported value")
	}
	
	# The matrix_element_type has an acceptable value
	biommatrixelemtypes = c("int", "float", "unicode")
	if( !object$matrix_element_type %in% biommatrixelemtypes ){
		return("matrix_element_type field has unsupported value")
	}	
	
	# If sparse, all data fields should have length of 3 (row, col, val)
	if( identical(object$matrix_type, "sparse") ){
		if( !all(sapply(object$data, length)==3) ){
			return("Some data fields for this sparse biom format do not have 3 elements")
		}
	}
	# If dense, data fields should have length equal to
	if( identical(object$matrix_type, "dense") ){
		if( !all(sapply(object$data, length)==bshape["ncol"]) ){
			return(paste("Some data fields for this dense biom format", 
									 "do not have the expected number columns:", bshape["ncol"], sep=" "))
		}
	}	
	
	# metadata ids -- rownames and colnames -- must be unique
	row_names = sapply(object$rows, function(i) i$id)
	dupr = duplicated(row_names)
	col_names = sapply(object$columns, function(i) i$id)
	dupc = duplicated(col_names)
	if( any(dupr) ){
		return(paste("The following row ids were duplicated: \n",
					paste0(row_names[dupr], collapse="\n"), sep=""))
	}
	if( any(dupc) ){
		return(paste("The following col ids were duplicated: \n",
								 paste0(col_names[dupc], collapse="\n"), sep=""))
	}	
	
	# If we get all the way here without returning a string,
	# then consider it a valid object.
	# In which case should return TRUE.
	return(TRUE)
}
################################################################################
## assign the function as the validity method for the otuTable class
setValidity("biom", validbiom)
################################################################################