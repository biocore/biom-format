.. _converting:

===============================
Converting between file formats
===============================

The ``convert_biom.py`` script in the biom-format project can be used to convert between biom and tab-delimited table formats. This is useful for several reasons:
 - converting biom format to classic OTU tables for easy viewing in programs such as Excel
 - converting between sparse and dense biom formats

Usage examples
--------------

Convert a tab-delimited table to sparse biom format::

	convert_biom.py -i otu_table.txt -o otu_table.biom --biom_table_type="otu table"

Convert a tab-delimited table to dense biom format::

	convert_biom.py -i otu_table.txt -o otu_table.biom --biom_table_type="otu table" --biom_type=dense

Convert biom format to tab-delimited table format::

	convert_biom.py -i otu_table.biom -o otu_table.txt -b

Convert dense biom format to sparse biom format::

	convert_biom.py -i otu_table.dense.biom -o otu_table.sparse.biom --dense_biom_to_sparse_biom

Convert sparse biom format to dense biom format::

	convert_biom.py -i otu_table.sparse.biom -o otu_table.dense.biom --sparse_biom_to_dense_biom

Convert biom format to classic format, including the ``taxonomy`` observation metadata as the last column of the classic format table::

	convert_biom.py -i otu_table.biom -o otu_table.txt -b --header_key taxonomy

Convert biom format to classic format, including the ``taxonomy`` observation metadata as the last column of the classic format table, but renaming that column as ``ConsensusLineage``::

	convert_biom.py -i otu_table.biom -o otu_table.txt -b --header_key taxonomy --output_metadata_id "ConsensusLineage"

.. _sparse-or-dense:

