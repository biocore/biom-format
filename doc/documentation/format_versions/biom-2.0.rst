.. _biom-2.0:

===========================================
The biom file format: Version 2.0
===========================================
    
The ``biom`` format is based on `HDF5 <http://www.hdfgroup.org>`_ to provide the overall structure for the format. HDF5 is a widely supported format with native parsers available within many programming languages. 

Required top-level attributes::

    id                  : <string or null> a field that can be used to id a table (or null)
    format              : <string> The name and version of the current biom format
    format-url          : <url> A string with a static URL providing format details
    type                : <string> Table type (a controlled vocabulary)
                          Acceptable values:
                           "OTU table"
                           "Pathway table"
                           "Function table"
                           "Ortholog table"
                           "Gene table"
                           "Metabolite table"
                           "Taxon table"
    generated-by        : <string> Package and revision that built the table
    creation-date       : <datetime> Date the table was built (ISO 8601 format)
    nnz                 : <int> The number of non-zero elements in the table
    shape               : <list of ints>, the number of rows and number of columns in data

Required groups::

    observation/        : The HDF5 group that contains observation specific information and an observation oriented view of the data
    observation/matrix  : The HDF5 group that contains matrix data oriented for observation-wise operations (e.g., in compressed sparse row format)
    sample/             : The HDF5 group that contains sample specific information and a sample oriented data oriented view of the data
    sample/matrix       : The HDF5 group that contains matrix data oriented for sample-wise operations (e.g., in compressed sparse column format)

Required datasets::

    observation/ids            : <string> or <variable length string> A (N,) dataset of the observation IDs, where N is the total number of IDs
    observation/matrix/data    : <float64> A (nnz,) dataset containing the actual matrix data
    observation/matrix/indices : <int32> A (nnz,) dataset containing the column indices (e.g., maps into samples/ids)
    observation/matrix/indptr  : <int32> A (M+1,) dataset containing the compressed row offsets
    sample/ids                 : <string> or <variable length string> A (M,) dataset of the sample IDs, where M is the total number of IDs
    sample/matrix/data         : <float64> A (nnz,) dataset containing the actual matrix data
    sample/matrix/indices      : <int32> A (nnz,) dataset containing the row indices (e.g., maps into observation/ids)
    sample/matrix/indptr       : <int32> A (N+1,) dataset containing the compressed column offsets

Optional datasets::

    observation/metadata       : <variable length string or null> If specified, a (1,) dataset containing a JSON-string representation of the metadata
    sample/metadata            : <variable length string or null> If specified, a (1,) dataset containing a JSON-string representation of the metadata


The metadata for each axis (observation and sample) are described with JSON. The required structure, if the metadata are specified, is a list of objects, where the list is in index order with respect to the axis (e.g, the object at element 0 corresponds to ID 0 for the given axis). Any metadata that corresponds to the ID, such as taxonomy, can be represented in the object. For instance, the following JSON string describes taxonomy for three IDs:

Metadata description::

    [
        {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}},
        {"taxonomy": ["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichospermum", "s__"]}},
        {"taxonomy": ["k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "g__Methanosarcina", "s__"]}}
    ]

Example biom files
==================

Below are examples of minimal and rich biom files in both sparse and dense formats. To decide which of these you should generate for new data types, see the section on :ref:`sparse-or-dense`.

BIOM 2.0 OTU table in the HDF5 data description langauge (DDL)
--------------------------------------------------------------

::

    HDF5 "rich_sparse_otu_table_hdf5.biom" {
    GROUP "/" {
       ATTRIBUTE "creation-date" {
          DATATYPE  H5T_STRING {
             STRSIZE H5T_VARIABLE;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_ASCII;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
          DATA {
          (0): "2014-05-13T14:50:32.052446"
          }
       }
       ATTRIBUTE "format-url" {
          DATATYPE  H5T_STRING {
             STRSIZE H5T_VARIABLE;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_ASCII;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
          DATA {
          (0): "http://biom-format.org"
          }
       }
       ATTRIBUTE "format-version" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
          DATA {
          (0): 2, 0
          }
       }
       ATTRIBUTE "generated-by" {
          DATATYPE  H5T_STRING {
             STRSIZE H5T_VARIABLE;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_ASCII;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
          DATA {
          (0): "example"
          }
       }
       ATTRIBUTE "id" {
          DATATYPE  H5T_STRING {
             STRSIZE H5T_VARIABLE;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_ASCII;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
          DATA {
          (0): "No Table ID"
          }
       }
       ATTRIBUTE "nnz" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
          DATA {
          (0): 15
          }
       }
       ATTRIBUTE "shape" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
          DATA {
          (0): 5, 6
          }
       }
       ATTRIBUTE "type" {
          DATATYPE  H5T_STRING {
             STRSIZE H5T_VARIABLE;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_ASCII;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
          DATA {
          (0): "otu table"
          }
       }
       GROUP "observation" {
          DATASET "ids" {
             DATATYPE  H5T_STRING {
                STRSIZE H5T_VARIABLE;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SIMPLE { ( 5 ) / ( 5 ) }
             DATA {
             (0): "GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"
             }
          }
          GROUP "matrix" {
             DATASET "data" {
                DATATYPE  H5T_IEEE_F64LE
                DATASPACE  SIMPLE { ( 15 ) / ( 15 ) }
                DATA {
                (0): 1, 5, 1, 2, 3, 1, 1, 4, 2, 2, 1, 1, 1, 1, 1
                }
             }
             DATASET "indices" {
                DATATYPE  H5T_STD_I32LE
                DATASPACE  SIMPLE { ( 15 ) / ( 15 ) }
                DATA {
                (0): 2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2
                }
             }
             DATASET "indptr" {
                DATATYPE  H5T_STD_I32LE
                DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
                DATA {
                (0): 0, 1, 6, 9, 13, 15
                }
             }
          }
          DATASET "metadata" {
             DATATYPE  H5T_STRING {
                STRSIZE H5T_VARIABLE;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
             DATA {
             (0): "[{"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}, {"taxonomy": ["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichospermum", "s__"]}, {"taxonomy": ["k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "g__Methanosarcina", "s__"]}, {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaerobium", "s__Halanaerobiumsaccharolyticum"]}, {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}]"
             }
          }
       }
       GROUP "sample" {
          DATASET "ids" {
             DATATYPE  H5T_STRING {
                STRSIZE H5T_VARIABLE;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
             DATA {
             (0): "Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
             (5): "Sample6"
             }
          }
          GROUP "matrix" {
             DATASET "data" {
                DATATYPE  H5T_IEEE_F64LE
                DATASPACE  SIMPLE { ( 15 ) / ( 15 ) }
                DATA {
                (0): 5, 2, 1, 1, 1, 1, 1, 1, 1, 2, 4, 3, 1, 2, 1
                }
             }
             DATASET "indices" {
                DATATYPE  H5T_STD_I32LE
                DATASPACE  SIMPLE { ( 15 ) / ( 15 ) }
                DATA {
                (0): 1, 3, 1, 3, 4, 0, 2, 3, 4, 1, 2, 1, 1, 2, 3
                }
             }
             DATASET "indptr" {
                DATATYPE  H5T_STD_I32LE
                DATASPACE  SIMPLE { ( 7 ) / ( 7 ) }
                DATA {
                (0): 0, 2, 5, 9, 11, 12, 15
                }
             }
          }
          DATASET "metadata" {
             DATATYPE  H5T_STRING {
                STRSIZE H5T_VARIABLE;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
             DATA {
             (0): "[{"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CGCTTATCGAGA", "Description": "human gut", "BODY_SITE": "gut"}, {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CATACCAGTAGC", "Description": "human gut", "BODY_SITE": "gut"}, {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTCTCTACCTGT", "Description": "human gut", "BODY_SITE": "gut"}, {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTCTCGGCCTGT", "Description": "human skin", "BODY_SITE": "skin"}, {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTCTCTACCAAT", "Description": "human skin", "BODY_SITE": "skin"}, {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTAACTACCAAT", "Description": "human skin", "BODY_SITE": "skin"}]"
             }
          }
       }
    }
    } 
