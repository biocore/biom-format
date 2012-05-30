.. _biom-1.0:

===========================================
The biom file format: Version 1.0
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
                           "unicode" : unicode string
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
    metadata            : <an object or null> A object containing key, value metadata pairs
  
The columns value is an ORDERED list of objects where each object corresponds to a single
column in the matrix. Each object can currently store arbitrary keys, although
this might become restricted based on table type. Each object must provide, 
at the minimum::
    
    id                  : <string> an arbitrary UNIQUE identifier
    metadata            : <an object or null> A object containing key, value metadata pairs

Example biom files
==================

Below are examples of minimal and rich biom files in both sparse and dense formats. To decide which of these you should generate for new data types, see the section on :ref:`sparse-or-dense`.

Minimal sparse OTU table
------------------------

::

    {
        "id":null,
        "format": "Biological Observation Matrix 0.9.1-dev",
        "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
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
        "format": "Biological Observation Matrix 0.9.1-dev",
        "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
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
     "format": "Biological Observation Matrix 0.9.1-dev",
     "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
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
     "format": "Biological Observation Matrix 0.9.1-dev",
     "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
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

