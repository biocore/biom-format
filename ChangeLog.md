BIOM-Format ChangeLog
=====================

biom 2.1.4
----------

Bug fixes, released on April 22nd 2015

Changes:

* Codebase updated to reflect pep8 1.6.x

Changes:

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
