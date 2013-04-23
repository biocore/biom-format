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

	met = matrix_element_type(object)
	if(!identical(length(met), 1L)){
		return("biom matrix_element_type field should have only 1 element")
	}
	if( !met %in% c("int", "float", "unicode") ){
		return("biom matrix_element_type value is unsupported.")
	}
	
	# Basic biom_shape value validity.
	bshape = biom_shape(object)
	if(!inherits(bshape, "integer")){
		return("problem with biom shape value type")
	}
	if(!identical(length(bshape), 2L)){
		return("problem with biom shape value length")
	}
	if(any(bshape < 0)){
		return("problem with biom shape value: negative value")
	}
	
	# Matrix shape and number of row/col elements matches
	if( ncol(object) > length(object$columns) ){
		return("shape field specifies more cols than are present in metadata")
	}
	if( nrow(object) > length(object$rows) ){
		return("shape field specifies more rows than are present in metadata")
	}
	if( ncol(object) < length(object$columns) ){
		return("more metadata cols than specified by shape field")
	}
	if( nrow(object) < length(object$rows) ){
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
		if( !all(sapply(object$data, length)==3L) ){
			return("Some data fields for this sparse biom format do not have 3 elements")
		}
	}
	# If dense, data fields should have length equal to ncol
	if( identical(object$matrix_type, "dense") ){
		if( !all(sapply(object$data, length)==ncol(object)) ){
			return(paste("Some data fields for this dense biom format", 
									 "do not have the expected number columns:", ncol(object), sep=" "))
		}
	}	
	
	# metadata ids -- rownames and colnames -- must be unique
	dupr = duplicated(rownames(object))
	dupc = duplicated(colnames(object))
	if( any(dupr) ){
		return(paste("The following row ids were duplicated: \n",
					paste0(rownames(object)[dupr], collapse="\n"), sep=""))
	}
	if( any(dupc) ){
		return(paste("The following col ids were duplicated: \n",
					paste0(colnames(object)[dupc], collapse="\n"), sep=""))
	}	
	
	# If we get all the way here without returning a string,
	# then consider it a valid object.
	# In which case should return TRUE.
	return(TRUE)
}
################################################################################
## assign the function as the validity method for the biom-class
setValidity("biom", validbiom)
################################################################################
