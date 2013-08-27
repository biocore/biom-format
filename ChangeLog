*********
ChangeLog
*********

BIOM-Format 1.1.2 - 1.2.0
=========================

New Features
------------

* Table.collapseObservationsByMetadata and Table.collapseSamplesByMetadata now have an additional argument, include_collapsed_metadata, which allows the user to either include or exclude collapsed metadata in the collapsed table.
* Table.collapseObservationsByMetadata and Table.collapseSamplesByMetadata now have an additional argument, one_to_many_mode, which allows the user to specify a collapsing strategy for one-to-many metadata relationships (currently supports adding and dividing counts).
* Table.binObservationsByMetadata, Table.binSamplesByMetadata, Table.collapseObservationsByMetadata, and Table.collapseSamplesByMetadata now have an additional argument, constructor, which allows the user to choose the return type of the binned/collapsed table(s).
* Table.delimitedSelf now has an additional argument, observation_column_name, which allows the user to specify the name of the first column in the output table (e.g. 'OTU ID', 'Taxon', etc.).
* Added new Table.transpose method.

Changes
-------

* Table.addSampleMetadata and Table.addObservationMetadata now support adding metadata to a subset of the samples/observations in a table that previously was without any sample/observation metadata. This used to result in an error.

Bug Fixes
---------

* Fixed performance issue with formatting BIOM tables for writing to a file.

BIOM-Format 1.1.1 - 1.1.2
=========================

New Features
------------

* Table.collapseObservationsByMetadata and Table.collapseSamplesByMetadata now
support one-to-many relationships on the metadata field to collapse on.

* added new script called print_biom_table_summary.py (and accompanying tutorial) that prints summary statistics of the input BIOM table as a whole and on a per-sample basis

Changes
-------

* SparseMat now uses cython for loops more efficiently

Bug Fixes
---------

* fixed serious performance issue with Table.transformSamples/Observations when using CSMat as the sparse backend

BIOM-Format 1.1.0 - 1.1.1
=========================

New Features
------------

Changes
-------

* added documentation for how to switch sparse backends via BIOM config file

Bug Fixes
---------

* performance issue on table creation with CSMat where an O(N) lookup was being performed

BIOM-Format 1.0.0 - 1.1.0
=========================

New Features
------------

* new default sparse matrix backend CSMat (COO/CSR/CSC) more efficient than SparseDict and SparseMat (pure python + numpy)
* support for biom config file, which allows specification of sparse backend to use. Currently supports CSMat (default), SparseMat, and SparseDict. Default can be found under support_files/biom_config, and can be copied to $HOME/.biom_config or located by setting $BIOM_CONFIG_FP
* new script called add_metadata.py with accompanying tutorial that allows users to add arbitrary sample and/or observation metadata to biom files
* new script called subset_biom.py that efficiently pulls out a subset of a biom table (either by samples or observations). Useful for very large tables where memory may be an issue

Changes
-------

* parser is more efficient for sparse tables and formatter is more efficient for both table types (less memory consumption)
* biom.Table objects are now immutable (except that metadata can still be added via addSampleMetadata/addObservationMetadata). __setitem__ and setValueByIds have been removed and SampleIds, ObservationIds, SampleMetadata, and ObservationMetadata members are now tuples as a result
* biom.Table object has a new method called getTableDensity()
* performance testing framework has been added for Table objects

Bug Fixes
---------

* convert_biom.py now converts dense tables to sparse tables (previously it didn't do anything)
* many misc. fixes to script help/documentation and docstrings (fixing typos, editing for clarity, etc.)

BIOM-Format 0.9.1 - 1.0.0
=========================

New Features
------------

* new default sparse matrix backend SparseMat (requires Cython) more efficient over existing SparseDict backend

Changes
-------

Bug Fixes
---------

BIOM-Format 0.9 - 0.9.1
=======================

* format now accepts unicode but does not accept str due to JSON parsing from Python
* specification for metadata is now either null or an object
* PySparse has been gutted, sparse matrix support is now through Table.SparseDict

New Features
------------

* more table types!

Changes
-------

* Table.getBioFormatJsonString() and similar methods now require a generatedby string

Bug Fixes
---------
