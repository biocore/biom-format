.. _biom_format:

===========================================
The biom file format
===========================================
    
The ``biom`` format is based on `JSON <http://www.json.org>`_ to provide the overall structure for the format. JSON is a widely supported format with native parsers available within many programming languages. 

Required top-level fields::

    id                  : <string or null> a field that can be used to id a table (or null)
    format              : <string> The name and version of the current biom format
    format_url          : <url> A string with a static URL providing format details
    type                : <string> Table type (a controlled vocabulary)
                          Acceptable values:
                           "OTU table"
                           "Pathway table"
                           "Function table"
                           "Ortholog table"
                           "Gene table"
                           "Metabolite table"
                           "Taxon table"
    generated_by        : <string> Package and revision that built the table
    date                : <datetime> Date the table was built (ISO 8601 format)
    rows                : <list of objects> An ORDERED list of obj describing the rows 
                          (explained in detail below)
    columns             : <list of objects> An ORDERED list of obj  describing the columns 
                          (explained in detail below)
    matrix_type         : <string> Type of matrix data representation (a controlled vocabulary)
                          Acceptable values:
                           "sparse" : only non-zero values are specified
                           "dense" : every element must be specified
    matrix_element_type : Value type in matrix (a controlled vocabulary)
                          Acceptable values:
                           "int" : integer
                           "float" : floating point
                           "str" : string
    shape               : <list of ints>, the number of rows and number of columns in data
    data                : <list of lists>, counts of observations by sample
                           if matrix_type is "sparse", [[row, column, value],
                                                        [row, column, value],
                                                        ...]
                           if matrix_type is "dense",  [[value, value, value, ...],
                                                        [value, value, value, ...],
                                                        ...]

Optional top-level fields::

    comment             : <string> A free text field containing any information that you
                           feel is relevant (or just feel like sharing)

The rows value is an ORDERED list of objects where each object corresponds to a single
row in the matrix. Each object can currently store arbitrary keys, although
this might become restricted based on table type. Each object must provide, 
at the minimum::
    
    id                  : <string> an arbitrary UNIQUE identifier
    metadata            : <list of objects or null> A object containing key, value metadata pairs
  
The columns value is an ORDERED list of objects where each object corresponds to a single
column in the matrix. Each object can currently store arbitrary keys, although
this might become restricted based on table type. Each object must provide, 
at the minimum::
    
    id                  : <string> an arbitrary UNIQUE identifier
    metadata            : <list of objects or null> A object containing key, value metadata pairs

Example biom files
==================

Below are examples of minimal and rich biom files in both sparse and dense formats. To decide which of these you should generate for new data types, see the section on :ref:`sparse-or-dense`.

Minimal sparse OTU table
------------------------

::

    {
        "id":null,
        "format": "Biological Observation Matrix v0.9",
        "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "type": "OTU table",
        "generated_by": "QIIME revision 1.4.0-dev",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],  
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [5, 6], 
        "data":[[0,2,1],
                [1,0,5],
                [1,1,1],
                [1,3,2],
                [1,4,3],
                [1,5,1],
                [2,2,1],
                [2,3,4],
                [2,4,2],
                [3,0,2],
                [3,1,1],
                [3,2,1],
                [3,5,1],
                [4,1,1],
                [4,2,1]
               ]
    }

Minimal dense OTU table
-----------------------

::

    {
        "id":null,
        "format": "Biological Observation Matrix v0.9",
        "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "type": "OTU table",
        "generated_by": "QIIME revision 1.4.0-dev",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],  
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],  
        "matrix_type": "dense",
        "matrix_element_type": "int",
        "shape": [5,6],
        "data":  [[0,0,1,0,0,0], 
                  [5,1,0,2,3,1],
                  [0,0,1,4,2,0],
                  [2,1,1,0,0,1],
                  [0,1,1,0,0,0]]
    }

Rich sparse OTU table
---------------------

::

    {
     "id":null,
     "format": "Biological Observation Matrix v0.9",
     "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
     "type": "OTU table",
     "generated_by": "QIIME revision 1.4.0-dev",
     "date": "2011-12-19T19:00:00",
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_2", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichospermum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaerobium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                 "BarcodeSequence":"CGCTTATCGAGA",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"gut",
                                 "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                 "BarcodeSequence":"CATACCAGTAGC",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"gut",
                                 "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                 "BarcodeSequence":"CTCTCTACCTGT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"gut",
                                 "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                 "BarcodeSequence":"CTCTCGGCCTGT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"skin",
                                 "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                 "BarcodeSequence":"CTCTCTACCAAT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"skin",
                                 "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                 "BarcodeSequence":"CTAACTACCAAT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"skin",
                                 "Description":"human skin"}}
                ],
     "matrix_type": "sparse",
     "matrix_element_type": "int",
     "shape": [5, 6], 
     "data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,5,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
            ]
    }


Rich dense OTU table
--------------------

::

    {
     "id":null,
     "format": "Biological Observation Matrix v0.9",
     "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
     "type": "OTU table",
     "generated_by": "QIIME revision 1.4.0-dev",
     "date": "2011-12-19T19:00:00",  
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_2", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichospermum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaerobium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}}
        ],  
     "columns":[
        {"id":"Sample1", "metadata":{
                                 "BarcodeSequence":"CGCTTATCGAGA",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"gut",
                                 "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                 "BarcodeSequence":"CATACCAGTAGC",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"gut",
                                 "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                 "BarcodeSequence":"CTCTCTACCTGT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"gut",
                                 "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                 "BarcodeSequence":"CTCTCGGCCTGT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"skin",
                                 "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                 "BarcodeSequence":"CTCTCTACCAAT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"skin",
                                 "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                 "BarcodeSequence":"CTAACTACCAAT",
                                 "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                 "BODY_SITE":"skin",
                                 "Description":"human skin"}}
                ],
     "matrix_type": "dense",
     "matrix_element_type": "int",
     "shape": [5,6],
     "data":  [[0,0,1,0,0,0], 
               [5,1,0,2,3,1],
               [0,0,1,4,2,0],
               [2,1,1,0,0,1],
               [0,1,1,0,0,0]]
    }


.. _converting:

Converting between file formats
===============================

The ``convert_biom.py`` script in the biom-format project can be used to convert between biom and tab-delimited table formats. This is useful for several reasons:
 - converting biom format to classic OTU tables for easy viewing in programs such as Excel
 - converting between sparse and dense biom formats

Usage examples
--------------

Convert a tab-delimited table to sparse biom format::

	convert_biom.py -i otu_table.txt -o otu_table.biom

Convert a tab-delimited table to dense biom format::

	convert_biom.py -i otu_table.txt -o otu_table.biom -t dense

Convert biom format to tab-delimited table format::

	convert_biom.py -i otu_table.biom -o otu_table.txt -b

Convert dense biom format to sparse biom format::

	convert_biom.py -i otu_table.dense.biom -o otu_table.sparse.biom --dense_biom_to_sparse_biom

Convert sparse biom format to dense biom format::

	convert_biom.py -i otu_table.sparse.biom -o otu_table.dense.biom --sparse_biom_to_dense_biom

.. _sparse-or-dense:

Should I generate sparse or dense biom files?
=============================================

In general, we recommend using the sparse format for your biom files. These will be a lot smaller than the dense format biom files when your data is sparse (i.e., more than 85% of your counts are zero). This is common for OTU tables and metagenome tables, and you'll want to investigate whether it's true for your data. If you currently format your data in tab-separated tables where observations are rows and samples are columns, you can format that file to be convertible to biom format with the ``convert_biom.py``. Here you can create dense and sparse formats, and see which file size is smaller. See the section on :ref:`converting`. 

Motivation for the BIOM format: reasons for the switch in QIIME
===============================================================

The BIOM format was motivation by several goals. First, to facilitate efficient handling and storage of large, sparse biological contingency tables; second, to support encapsulation of core study data (contingency table data and sample/observation metadata) in a single file; and third, to facilitate the use of these tables between tools that support this format (e.g., passing of data between QIIME, MG-RAST, and VAMPS).

Efficient handling and storage of very large tables
-------------------------------------------------------

In QIIME, we began hitting limitations with OTU table objects when working with thousands of samples and hundreds of thousands of OTUs. In the near future we expect that we'll be dealing with hundreds of thousands of samples in single analyses.

The OTU table format up to QIIME 1.4.0 involved a dense matrix: if an OTU was not observed in a given sample, that would be indicated with a zero. We now primarily represent OTU tables in a sparse format: if an OTU is not observed in a sample, there is no count for that OTU. The two ways of representing this data are exemplified here. 

A dense representation of an OTU table:: 

   OTU ID PC.354  PC.355  PC.356  
   OTU0   0   0   4
   OTU1   6   0   0
   OTU2   1   0   7
   OTU3   0   0   3

A sparse representation of an OTU table::

    PC.354 OTU1 6
    PC.354 OTU2 1
    PC.356 OTU0 4
    PC.356 OTU2 7
    PC.356 OTU3 3

OTU table data tends to be sparse (e.g., greater than 90% of counts are zero, and frequently as many as 99% of counts are zero) in which case the latter format is more convenient to work with as it has a smaller memory footprint. Both of these representations are supported in the biom-format project via dense and sparse Table types. Generally if less than 85% of your counts are zero, a dense representation will be more efficient.

Encapsulation of core study data (OTU table data and sample/OTU metadata) in a single file
------------------------------------------------------------------------------------------

The JSON-format OTU table allow for storage of arbitrary amounts of sample and OTU metadata in a single file. Sample metadata corresponds to what is generally found in QIIME mapping files. At this stage inclusion of this information in the OTU table file is optional, but it may be useful for sharing these files with other QIIME users and for publishing or archiving results of analyses. OTU metadata (generally a taxonomic assignment for an OTU) is also optional. In contrast to the previous OTU table format, you can now store more than one OTU metadata value in this field, so for example you can score taxonomic assignments based on two different taxonomic assignment approaches.

Facilitating the use of tables between tools that support this format
---------------------------------------------------------------------

Different tools, such as `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_ work with similar data structures that represent different types of data. An example of this is a `metagenome` table that could be generated by MG-RAST (where for example, columns are metagenomes and rows are functional categories). Exporting this data from MG-RAST in a suitable format will allow for the application of many of the QIIME tools to this data (such as generation of alpha rarefaction plots or beta diversity ordination plots). This new format is far more general than previous formats, so will support adoption by groups working with different data types and is already being integrated to support transfer of data between `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_.
