BIOM-Format ChangeLog
=====================

biom 2.1.8
----------

New features and bug fixes, released on 28 January 2020.

Important:

* Python 2.7 and 3.5 support has been dropped.
* Python 3.8 support has been added into Travis CI. 
* A change to the defaults for `Table.nonzero_counts` was performed such that the default now is to count the number of nonzero features. See [issue #685](https://github.com/biocore/biom-format/issues/685)
* We now require a SciPy >= 1.3.1. See [issue #816](https://github.com/biocore/biom-format/issues/816)

New Features:

* The detailed report is no longer part of the table validator. See [issue #378](https://github.com/biocore/biom-format/issues/378).
* `load_table` now accepts open file handles. See [issue #481](https://github.com/biocore/biom-format/issues/481).
* `biom export-metadata` has been added to export metadata as TSV. See [issue #820](https://github.com/biocore/biom-format/issues/820).
* `Table.to_tsv` has been modified to allow for `direct_io`. See [issue #836](https://github.com/biocore/biom-format/pull/836).

Bug fixes:

* `Table.to_dataframe(dense=False)` does now correctly produce sparse data frames (and not accidentally dense ones as before). See [issue #808](https://github.com/biocore/biom-format/issues/808).
* Order of error evaluations was unstable in Python versions without implicit `OrderedDict`. See [issue #813](https://github.com/biocore/biom-format/issues/813). Thanks @gwarmstrong for identifying this bug.
* `Table._extract_data_from_tsv` would fail if taxonomy was provided, and if the first row had the empty string for taxonomy. See [issue #827](https://github.com/biocore/biom-format/issues/827). Thanks @KasperSkytte for identifying this bug.

biom 2.1.7
----------

New features and bug fixes, released on 28 September 2018.

Important:

* Python 3.4 support has been dropped. We now only support Python 2.7, 3.5, 3.6 and 3.7.
* We will be dropping Python 2.7 support on the next release.
* Pandas >= 0.20.0 is now the minimum required version.
* pytest is now used instead of nose.

New Features:

* Massive performance boost to `Table.collapse` with the default collapse function. The difference was 10s of milliseconds vs. minutes stemming from prior use of `operator.add`. See [issue #761](https://github.com/biocore/biom-format/issues/761).
* `Table.align_to` for aligning one table to another. This is useful in multi-omic analyses where multiple preparations have been performed on the sample physical samples. This is essentially a helper method around `Table.sort_order`. See [issue #747](https://github.com/biocore/biom-format/issues/747).
* Added additional sanity checks when calling `Table.to_hdf5`, see [PR #769](https://github.com/biocore/biom-format/pull/769).
* `Table.subsample()` can optionally perform subsampling with replacement. See [issue #774](https://github.com/biocore/biom-format/issues/774).
* `Table.to_dataframe()` now supports a `dense` argument to return `pd.DataFrame`. See [issue #762](https://github.com/biocore/biom-format/issues/762).
* Parsing methods for BIOM-Format 1.0.0 tables now preserve dict ordering. See [issue #781](https://github.com/biocore/biom-format/issues/781).

Bug fixes:

* `Table.subsample(by_id=True, axis='observation')` did not subsample over the 'observations'. Because of the nature of the bug, an empty table was returned, so the scope of the issue is such that it should not have produced misleading results but instead triggered empty table errors, with the exception of the pathological case of the ID namespaces between features and samples not being disjoint. See [PR #759](https://github.com/biocore/biom-format/pull/759) for more information.
* Tables of shape `(0, n)` or `(n, 0)` were raising exceptions when being written out. See [issue #619](https://github.com/biocore/biom-format/issues/619).
* Tables loaded with a `list` of empty `dict`s will have their metadata attributes set to None. See [issue #594](https://github.com/biocore/biom-format/issues/594).

biom 2.1.6
----------

New features and bug fixes, released on 27 April 2017.

New Features:

* `Table.from_hdf5` now supports a rapid subset in the event that metadata is
   not needed. In benchmarking against the Earth Microbiome Project BIOM table,
   the reduction in runtime was multiple orders of magnitude while additionally
   preserving substantial memory.
* `Table.rankdata` has been added to convert values to ranked abundances on
  either axis. See [issue #645](https://github.com/biocore/biom-format/issues/639).
* Format of numbers in ``biom summarize-table`` output is now more readable and localized. See [issue #679](https://github.com/biocore/biom-format/issues/679).
* `Table.concat` has been added to the API and allows for concatenating multiple tables in which the IDs of one of the axes are known to be disjoint. This has substantial performance benefits over `Table.merge`.
* `Table.sort_order` was performing an implicit cast to dense, and not leveraging fancy indexing. A substantial performance gain was acheived. See [PR #720](https://github.com/biocore/biom-format/pull/720).
* `biom subset-table` now accepts a QIIME-like mapping file when subsetting by IDs [Issue #587](https://github.com/biocore/biom-format/issues/587)
* `Table.del_metadata` was added to support the removal of metadata entries from the table [Issue #708](https://github.com/biocore/biom-format/issues/708).
* `Table.to_dataframe` was added to cast the internal matrix data to a Pandas `SparseDataFrame` [Issue #622](https://github.com/biocore/biom-format/issues/622).
* `Table.metadata_to_dataframe` was added to cast axis metadata to a Pandas `DataFrame` [Issue #622](https://github.com/biocore/biom-format/issues/622).
* `test_table.py` and `test_util.py` now use a stable random seed. See issue [#728](https://github.com/biocore/biom-format/issues/728)
* Failure to cast a value when parsing a TSV will now print the associated line number which had the bad value. See [#284](https://github.com/biocore/biom-format/issues/284).
* `Table.remove_empty` has been added to remove zero'd samples, observations or both. See [#721](https://github.com/biocore/biom-format/issues/721).
* A subcommand of the command line interface was added to obtain a table's IDs: `table-ids`.

Bug fixes:

* ``-o`` is now a required parameter of ``biom from-uc``. This was not the case previously, which resulted in a cryptic error message if ``-o`` was not provided. See [issue #683](https://github.com/biocore/biom-format/issues/683).
* Matrices are now cast to csr on `Table` construction if the data evaluate as `isspmatrix`. This fixes [#717](https://github.com/biocore/biom-format/issues/717) where some API methods assumed the data were csc or csr.
* `Table.concat` was not handling tables without metadata, resulting in an exception due to mismatches metadata shape. See [#724](https://github.com/biocore/biom-format/issues/724).
* When validating a BIOM-Format 1.0.0 table, specifying the version string would trigger an error. See [#664](https://github.com/biocore/biom-format/issues/664). An explicit regression test was not added as this stemmed from an integration, and there currently is not support for script usage tests; see [#656](https://github.com/biocore/biom-format/issues/656).
* `Table.nnz` was not calling `eliminate_zeros()` on the underlying sparse matrix, resulting in the possibility of counting explicitly set zero values. See [#727](https://github.com/biocore/biom-format/issues/727).
* `Table.from_hdf5` was not properly turning `bytes` into `str` for the `table_id` and the `type` HDF5 attributes. See [#731](https://github.com/biocore/biom-format/issues/731).
* `Table.__init__` now always performs an `astype(float)` on the contained `spmatrix`. This type normalization is beneficial for underlying Cython code on the filtering and transform operations. It is possible this will introduce some performance overhead, however in _most_ cases the data should already be float. See [#718](https://github.com/biocore/biom-format/issues/718).
* `Table.to_hdf5` was not handling lists of str appropriately in the general case. Ssee [#638](https://github.com/biocore/biom-format/issues/638).
* `Table.to_hdf5` was not handling taxonomy as flat strings, which was a common mistake that was outside of expectations for the formatter. The formatter now attempts to split on semicolon if this scenario is encountered, and errors with a more informative error if a problem occurs. See [#530](https://github.com/biocore/biom-format/issues/530).

biom 2.1.5
----------

New features and bug fixes, released on 21 October 2015.

Changes:

* Codebase is now Python 2/3 compatible. It is currently tested with Python
  versions 2.7, 3.4 and 3.5.
* `biom-serve` and the accompanying html interface has been removed.

New Features:

* `Table.head` has been added to retrieve the first few rows and or columns
  from a table. This can be accessed through the new ``biom head`` command.
  See [issue #639](https://github.com/biocore/biom-format/issues/639).
* ``biom.parse.from_uc`` has been added to support creation of ``biom.Table``
  objects from vsearch/uclust/usearch ``.uc`` files. This can be accessed
  through the new ``biom from-uc`` command. See
  [issue #648](https://github.com/biocore/biom-format/issues/648).
* Codebase now uses [click](http://click.pocoo.org) instead of
  [pyqi](https://github.com/biocore/pyqi) for its command line interface.
  See [issue #631](https://github.com/biocore/biom-format/issues/631).

Bug fixes:

* `Table.update_ids` strict check was too aggressive. See
 [issue #633](https://github.com/biocore/biom-format/issues/633).
* `biom --version` now prints the software version (previously the individual
  commands did this, but not the base command).
* `Table.vlen_list_of_str_formatter` was considering a `str` to be valid for
  formatting resulting in an obscure error when a `str`, as opposed to a
  `list` of `str`, was used for taxonomy. See
  [issue #709](https://github.com/biocore/biom-format/issues/709).

biom 2.1.4
----------

Bug fixes, released on April 22nd 2015

Changes:

* Codebase updated to reflect pep8 1.6.x

New features:

* `Table.to_hdf5` and `Table.from_hdf5` now support custom parsers and
    formatters, see issue #608

Bug fixes:

* `Table.update_ids` was not updating the internal ID lookup caches, issue #599
* `--is-json` has been removed from the table validator as it was being ignored
    anyway, issue #591
* `biom summarize-table` can now properly interact with pipes. This previously
    worked on OSX but did not on Linux. Issue #598
* `biom convert` was recording the wrong version information from HDF5 -> JSON,
    issue #595
* `Table.collapse`, under `one_to_many` was not constructing the resulting
    matrix properly, issue #606
* Improve error message when trying to load an empty file, issue #614.
* Improve error handling when filtering tables, and return tables of shape
    `(0, n)` instead of `(0, 0)` when fully filtering out a table along an
    axis, issue #620
* Fix `Table.nonzero` to work on data that is not already in csr, issue #625.

biom 2.1.3
----------

Minor fixes, released on January 29, 2014

Bug fixes:

* Improve error message when trying to load an HDF5 file without h5py being
    installed.
* Allow validating json files when h5py is not installed.

biom 2.1.2
----------

Minor fixes, released on December 18, 2014

Bug fixes:

* Remove syntax error from `normalize_table.py`.
* `Table.to_json` was not serializing empty tables properly, see #571
* `biom summarize-table` could not handle empty tables, see #571

biom 2.1.1
----------

Minor fixes and performance improvements, released on November 19th 2014

Changes:

* The collapsing function to `Table.collapse` is now passed the entire table to
    allow for more complex collapses (e.g., median, random selection, etc). See
    #544, #545 and #547.
* Updated version strings in the project to be
    [Semantic Versioning](www.semver.org)-stlye. This better matches with other
    open source python projects, and plays nicer with pip.
* Conversion from TSV now takes less memory. See #551.
* Parameter header_mark has been removed from _extract_data_from_tsv()
    in table.py
* Order of magnitude improvement in parsing HDF5 BIOM tables, see #529
* Added `Table.length`, see #548
* Order of magnitude performance increase in `Table.nonzero`, see #538

Bug fixes:

* Ensure that a copy is performed in `Table.subsample`
* Avoided a memory leak when checking if a table is JSON or TSV, see #552.

biom 2.1
--------

Format finalization, released on August 7th 2014

New features:

* Group metadata (e.g., a phylogenetic tree) can now be stored within the HDF5
    representation. These data are available within the `Table` object
* Matrix data can now be accessed by the ``Table.matrix_data`` property
* ``Table`` IDs are now accessed via the ``Table.ids`` method
* ``Table`` metadata are now accessed via the ``Table.metadata`` method
* New method ``Table.update_ids``, which allows for updating the ids along
    either axis.
* added ``normalize-table`` option to optparse and HTML interfaces which
    utilizes the new TableNormalizer command from ``table_normalizer.py``

Changes:

* Metadata are now stored in individual datasets within HDF5. This resulted in
    a change to the BIOM-Format spec which has now been bumped to format
    version 2.1.
* ``Table.collapse`` ``min_group_size`` is now 1 by default, see #480
* General improvements to BIOM 2.x online documentation
* ``Table.pa`` now supports negative values
* dropped old, unused scripts
* added ``Table.iter_pairwise``
* added ``Table.min`` and ``Table.max``, see #459
* iter methods now support dense/sparse
* added ``Table.matrix_data`` property
* ``Table.filter`` yields a sparse vector, see #470
* ``Table.subsample`` can now sample by IDs (e.g., get a random subset of
    samples or observations from a ``Table``).
* ``biom.util.generate_subsamples`` will generate an infinite number of
    subsamples and can be used for rarefaction.
* ``biom summarize-table`` can now operate on observations.
* 10% performance boost in ``Table.subsample``, see #532

Bug fixes:

* ``Table.transform`` operates on full vectors now, see #476
* ``biom convert`` now handles taxonomy strings correctly, see #504
* ``Table.sort_order`` was not retaining ``Table.type``, see #474
* ``convert_biom_to_table`` now uses ``load_table``, see #478
* ``Table.pa`` now handles negative values, see #492
* ``Table.copy`` was not retaining ``Table.type``, see #494

biom 2.0.1
----------

Bug fix release, released on June 3rd 2014

Changes:

* Light weight loading mechanism (`biom.load_table`) added
* `Table.data` now has a default axis
* Convert documentation updated
* Quick start page added to documentation

Bug fixes:

* missing fields from JSON representation reintroduced
* `TableConverter` works as expected

biom 2.0.0
----------

Major release, released on May 15th 2014

Changes:

* NumPy 1.7 or above is required
* Support for HDF5
* Codebase is PEP-8 compliant
* CSMat has been removed and Scipy is now a required dependency
* Requires pyqi 0.3.2
* New HTML interface
* No longer dependent on dateutil
* `Table.bin_samples_by_metadata` and `Table.bin_observations_by_metadata` have
    been combined into `Table.partition`, which takes an axis argument
* `Table.collapse_samples_by_metadata` and
    `Table.collapse_observations_by_metadata` have been combined into
    `Table.collapse`, which now takes an axis argument
* `Table.filter_samples` and `Table.filter_observations` have been combined
    into `Table.filter`, which now takes an axis argument
* `Table.transform_samples` and `Table.transform_observations` have been
    combined into `Table.transform`, which now takes an axis argument
* `Table.norm_sample_by_observation` and `Table.norm_observation_by_sample`
    have been combined into `Table.norm`, which now takes an axis argument
* `Table.iter_samples` and `Table.iter_observations` have been combined into
    `Table.iter`, which now takes an axis argument
* `Table.iter_sample_data` and `Table.iter_observation_data` have been combined
    into `Table.iter_data`, which now takes an axis argument
* `Table.get_sample_index` and `Table.get_observation_index` have been combined
    into `Table.get_index`, which now takes an axis argument
* `Table.add_sample_metadata` and `Table.add_observation_metadata` have been
    combined into `Table.add_metadata`, which now takes an axis argument
* `Table.sample_data` and `Table.observation_data` have been combined into
    `Table.data`, which now takes an axis argument
* `Table.sample_exists` and `Table.observation_exists` have been combined into
    `Table.exists`, which now takes an axis argument
* `Table.sort_by_sample_ids` and `Table.sort_by_observation_ids` have been
    combined into `Table.sort`, which now takes an axis argument
* `Table.sort_sample_order` and `Table.sort_observation_order` have been
    combined into `Table.sort_order`, which now takes an axis argument
* `Table.norm_samples_by_metadata` and `Table.norm_observations_by_metadata`
    have been removed
* Added `Table.metadata` to allow fetching of metadata by an ID instead of just
    by index
* Added `Table.pa` for conversion to presence/absence
* Added `Table.subsample` for randomly subsampling data
* `Table` now embraces numpydoc

biom 1.3.1
----------

Documentation release, released on December 4th 2013

New Features:

* biom-format is now installable via pip! Simply run ``pip install biom-format``.

Changes:

* Fixed installation instructions to be clearer about the various ways of installing biom-format. Also fixed a couple of minor formatting issues.

biom 1.3.0
----------

Feature release, released on December 4th 2013

New Features:

* Added new sparse matrix backend ``ScipySparseMat``, which requires that [scipy](http://www.scipy.org/) is installed if this backend is in use. This backend will generally yield improvements in both runtime and memory consumption, especially with larger sparse tables. The default sparse matrix backend is still ``CSMat`` (this means that scipy is an optional dependency of the biom-format project).

Changes:

* Sparse backends ```SparseDict``` and ```SparseMat``` have been removed in favor of ```CSMat```. Cython is no longer a dependency.
* The BIOM Format project license is now Modified BSD (see COPYING.txt for more details) and is no longer GPL. To change the license, we obtained written permission (by email) from all past and present developers on the biom-format project. The core developers, including @gregcaporaso, @wasade, @jrrideout, and @rob-knight were included on these emails. For code that was derived from the QIIME and PyCogent projects, which are under the GPL license, written permission was obtained (by email) from the developers of the original code (tracing through the commit history, as necessary). @gregcaporaso, @wasade, @jrrideout, and @rob-knight were included on these emails.
* Removed the top-level ```python-code``` directory, moving all contents up one level. If you are installing the biom-format project by manually setting ```PYTHONPATH``` to ```<dir prefix>/biom-format/python-code```, you will need to change the path to ```<dir prefix>/biom-format``` instead. Please see the installation instructions for more details.
* Reorganized sparse backend code into a new subpackage, ```biom.backends```. This change should not affect client code.

biom 1.2.0
----------

New Features:

* ```Table.collapseObservationsByMetadata``` and ```Table.collapseSamplesByMetadata``` now have an additional argument, ```include_collapsed_metadata```, which allows the user to either include or exclude collapsed metadata in the collapsed table.
* ```Table.collapseObservationsByMetadata``` and ```Table.collapseSamplesByMetadata``` now have an additional argument, ```one_to_many_mode```, which allows the user to specify a collapsing strategy for one-to-many metadata relationships (currently supports adding and dividing counts).
* ```Table.binObservationsByMetadata```, ```Table.binSamplesByMetadata```, ```Table.collapseObservationsByMetadata```, and ```Table.collapseSamplesByMetadata``` now have an additional argument, ```constructor```, which allows the user to choose the return type of the binned/collapsed table(s).
* ```Table.delimitedSelf``` now has an additional argument, ```observation_column_name```, which allows the user to specify the name of the first column in the output table (e.g. 'OTU ID', 'Taxon', etc.).
* Added new ```Table.transpose``` method.
* ```Table.__init``` has change from ```__init__(self, data, sample_ids, observation_ids, sample_metadata=None,
observation_metadata=None, table_id=None, type=None, **kwargs)``` to ```__init__(self, data, observation_ids, sample_ids, observation_metadata=None, sample_metadata=None, table_id=None, type=None, **kwargs)``` This is for clarity, the data is in the same order as the arguments to the constructor.
*```table_factory``` has changed from ```table_factory(data, sample_ids, observation_ids, sample_metadata=None, observation_metadata=None, table_id=None, input_is_dense=False, transpose=False, **kwargs)``` to ```table_factory(data, observation_ids, sample_ids, observation_metadata=None, sample_metadata=None, table_id=None, input_is_dense=False, transpose=False, **kwargs)``` This is for clarity, the data is in the same order as the arguments to the function.

Changes:

* [pyqi](http://bipy.github.io/pyqi) 0.2.0 is now a required dependency. This changes the look-and-feel of the biom-format command-line interfaces and introduces a new executable, ```biom```, which can be used to see a list of all available biom-format command-line commands. The ```biom``` command is now used to run biom-format commands, instead of having a Python script (i.e., .py file) for each biom-format command. The old scripts (e.g., add_metadata.py, convert_biom.py, etc.) are still included but are deprecated. Users are pointed to the new ```biom``` command to run instead. Bash tab completion is now supported for all command and option names (see the biom-format documentation for instructions on how to enable this).
* The following scripts have had their names and options changed:
    * ```add_metadata.py``` is now ```biom add-metadata```. Changed option names:
        * ```--input_fp``` is now ```--input-fp```
        * ```--output_fp``` is now ```--output-fp```
        * ```--sample_mapping_fp``` is now ```--sample-metadata-fp```
        * ```--observation_mapping_fp``` is now ```--observation-metadata-fp```
        * ```--sc_separated``` is now ```--sc-separated```
        * ```--int_fields``` is now ```--int-fields```
        * ```--float_fields``` is now ```--float-fields```
        * ```--sample_header``` is now ```--sample-header```
        * ```--observation_header``` is now ```--observation-header```
        * New option ```--sc-pipe-separated```
    * ```biom_validator.py``` is now ```biom validate-table```. Changed option names:
        * ```-v```/```--verbose``` is now ```--detailed-report```
        * ```--biom_fp``` is now ```--input-fp```
    * ```convert_biom.py``` is now ```biom convert```. Changed option names:
        * ```--input_fp``` is now ```--input-fp```
        * ```--output_fp``` is now ```--output-fp```
        * ```--biom_type``` is now ```--matrix-type```
        * ```--biom_to_classic_table``` is now ```--biom-to-classic-table```
        * ```--sparse_biom_to_dense_biom``` is now ```--sparse-biom-to-dense-biom```
        * ```--dense_biom_to_sparse_biom``` is now ```--dense-biom-to-sparse-biom```
        * ```--sample_mapping_fp``` is now ```--sample-metadata-fp```
        * ```--observation_mapping_fp``` is now ```--observation-metadata-fp```
        * ```--header_key``` is now ```--header-key```
        * ```--output_metadata_id``` is now ```--output-metadata-id```
        * ```--process_obs_metadata``` is now ```--process-obs-metadata```
        * ```--biom_table_type``` is now ```--table-type```
    * ```print_biom_python_config.py``` is now ```biom show-install-info```.
    * ```print_biom_table_summary.py``` is now ```biom summarize-table```. Changed option names:
        * ```--input_fp``` is now ```--input-fp```
        * ```--output_fp``` is now ```--output-fp```. This is now a required option (output is no longer printed to stdout).
        * ```--num_observations``` is now ```--qualitative```
        * ```--suppress_md5``` is now ```--suppress-md5```
    * ```subset_biom.py``` is now ```biom subset-table```. Changed option names:
        * ```--biom_fp``` is now ```--input-fp```
        * ```--output_fp``` is now ```--output-fp```
        * ```--ids_fp``` is now ```--ids```
* ```biom.parse.parse_mapping``` has been replaced by ```biom.parse.MetadataMap```. ```biom.parse.MetadataMap.from_file``` can be directly substituted in place of ```biom.parse.parse_mapping```.

Bug Fixes:

* Fixed performance issue with formatting BIOM tables for writing to a file.
* Fixed issue with ```Table.addSampleMetadata``` and ```Table.addObservationMetadata``` when adding metadata to a subset of the samples/observations in a table that previously was without any sample/observation metadata.
* Fixed issue with ```Table.addSampleMetadata``` and ```Table.addObservationMetadata``` when updating a table's existing metadata, including the case where there are sample/observation IDs that are in the metadata file but not in the table.

biom 1.1.2
----------

New Features:

* ```Table.collapseObservationsByMetadata``` and ```Table.collapseSamplesByMetadata``` now
support one-to-many relationships on the metadata field to collapse on.

* added new script called ```print_biom_table_summary.py``` (and accompanying tutorial) that prints summary statistics of the input BIOM table as a whole and on a per-sample basis

Changes:

* ```SparseMat``` now uses cython for loops more efficiently

Bug Fixes:

* fixed serious performance issue with ```Table.transformSamples/Observations``` when using ```CSMat``` as the sparse backend

biom 1.1.1
----------

Changes:

* added documentation for how to switch sparse backends via BIOM config file

Bug Fixes:

* performance issue on table creation with ```CSMat``` where an ```O(N)``` lookup was being performed

biom 1.1.0
----------

New Features:

* new default sparse matrix backend ```CSMat``` (COO/CSR/CSC) more efficient than ```SparseDict``` and ```SparseMat``` (pure python + numpy)
* support for biom config file, which allows specification of sparse backend to use. Currently supports ```CSMat``` (default), ```SparseMat```, and ```SparseDict```. Default can be found under ```support_files/biom_config```, and can be copied to ```$HOME/.biom_config``` or located by setting ```$BIOM_CONFIG_FP```
* new script called ```add_metadata.py``` with accompanying tutorial that allows users to add arbitrary sample and/or observation metadata to biom files
* new script called ```subset_biom.py``` that efficiently pulls out a subset of a biom table (either by samples or observations). Useful for very large tables where memory may be an issue

Changes:

* parser is more efficient for sparse tables and formatter is more efficient for both table types (less memory consumption)
* ```biom.Table``` objects are now immutable (except that metadata can still be added via ```addSampleMetadata```/```addObservationMetadata```). ```__setitem__``` and ```setValueByIds``` have been removed and ```SampleIds```, ```ObservationIds```, ```SampleMetadata```, and ```ObservationMetadata``` members are now tuples as a result
* ```biom.Table``` object has a new method called ```getTableDensity()```
* performance testing framework has been added for ```Table``` objects

Bug Fixes:

* ```convert_biom.py``` now converts dense tables to sparse tables (previously it didn't do anything)
* many misc. fixes to script help/documentation and docstrings (fixing typos, editing for clarity, etc.)

biom 1.0.0
----------

New Features:

* new default sparse matrix backend ```SparseMat``` (requires Cython) more efficient over existing ```SparseDict``` backend

biom 0.9.1
----------

* format now accepts unicode but does not accept str due to JSON parsing from Python
* specification for metadata is now either ```null``` or an object
* PySparse has been gutted, sparse matrix support is now through ```Table.SparseDict```

New Features:

* more table types!

Changes:

* ```Table.getBioFormatJsonString()``` and similar methods now require a ```generatedby``` string
