BIOM-Format ChangeLog
=====================

biom 1.3.1-dev
--------------

Changes:

* NumPy 1.7 or above is required
* Support for HDF5
* Codebase is PEP-8 compliant
* CSMat has been removed and Scipy is now a required dependency
* Requires PyQi 0.3.1-dev
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
