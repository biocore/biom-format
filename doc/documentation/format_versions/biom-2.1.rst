.. _biom-2.1:

===========================================
The biom file format: Version 2.1
===========================================
    
The ``biom`` format is based on `HDF5 <http://www.hdfgroup.org>`_ to provide the overall structure for the format. HDF5 is a widely supported format with native parsers available within many programming languages. 

Required top-level attributes::

    id                   : <string or null> a field that can be used to id a table (or null)
    type                 : <string> Table type (a controlled vocabulary)
                           Acceptable values:
                            "OTU table"
                            "Pathway table"
                            "Function table"
                            "Ortholog table"
                            "Gene table"
                            "Metabolite table"
                            "Taxon table"
    format-url           : <url> A string with a static URL providing format details
    format-version       : <tuple> The version of the current biom format, major and minor
    generated-by         : <string> Package and revision that built the table
    creation-date        : <datetime> Date the table was built (ISO 8601 format)
    shape                : <list of ints>, the number of rows and number of columns in data
    nnz                  : <int> The number of non-zero elements in the table

Required groups::

    observation/               : The HDF5 group that contains observation specific information and an observation oriented view of the data
    observation/matrix         : The HDF5 group that contains matrix data oriented for observation-wise operations (e.g., in compressed sparse row format)
    observation/metadata       : The HDF5 group that contains observation specific metadata information
    observation/group-metadata : The HDF5 group that contains observation specific group metadata information (e.g., phylogenetic tree)
    sample/                    : The HDF5 group that contains sample specific information and a sample oriented data oriented view of the data
    sample/matrix              : The HDF5 group that contains matrix data oriented for sample-wise operations (e.g., in compressed sparse column format)
    sample/metadata            : The HDF5 group that contains sample specific metadata information
    sample/group-metadata      : The HDF5 group that contains sample specific group metadata information (e.g., relationships between samples)

Required datasets::

    observation/ids            : <string> or <variable length string> A (N,) dataset of the observation IDs, where N is the total number of IDs
    observation/matrix/data    : <float64> A (nnz,) dataset containing the actual matrix data
    observation/matrix/indices : <int32> A (nnz,) dataset containing the column indices (e.g., maps into samples/ids)
    observation/matrix/indptr  : <int32> A (M+1,) dataset containing the compressed row offsets
    sample/ids                 : <string> or <variable length string> A (M,) dataset of the sample IDs, where M is the total number of IDs
    sample/matrix/data         : <float64> A (nnz,) dataset containing the actual matrix data
    sample/matrix/indices      : <int32> A (nnz,) dataset containing the row indices (e.g., maps into observation/ids)
    sample/matrix/indptr       : <int32> A (N+1,) dataset containing the compressed column offsets

Under the ``observation/metadata`` and ``sample/metadata`` groups, the user can specify an arbitrary number of datasets that represents a metadata category for that axis. The expected structure for each of these metadata datasets is a list of atomic type objects (int, float, str, ...) where the index order of the list corresponds to the index order of the relevant axis IDs. Special complex metadata fields have been defined, and they are stored in a specific way. Currently, the available special metadata fields are::

    observation/metadata/taxonomy      : <string> or <variable length string> A (N, ?) dataset containing the taxonomy names assigned to the observation
    observation/metadata/KEGG_Pathways : <string> or <variable length string> A (N, ?) dataset containing the KEGG Pathways assigned to the observation
    observation/metadata/collapsed_ids : <string> or <variable length string> A (N, ?) dataset containing the observation ids of the original table that have been collapsed in the given observation
    sample/metadata/collapsed_ids      : <string> or <variable length string> A (M, ?) dataset containing the sample ids of the original table that have been collapsed in the given sample

Under the ``observation/group-metadata`` and ``sample/group-metadata`` groups, the user can specify an arbitrary number of datasets that represents a relationship between the ids for that axis. The expected structure for each of these group metadata datasets is a single string or variable length string. Each of these datasets should have defined an attribute called ``data_type``, which specifies how the string should be interpreted. One example of such group metadata dataset is ``observation/group-metadata/phylogeny``, with the attribute ``observation/group-metadata/phylogeny.attrs['data_type'] = "newick"``, which stores a single string with the newick format of the phylogenetic tree for the observations.


Example biom files
==================

Below is an examples of a rich biom file. To decide which of these you should generate for new data types, see the section on :ref:`sparse-or-dense`.

BIOM 2.1 OTU table in the HDF5 data description langauge (DDL)
--------------------------------------------------------------

::

    HDF5 "examples/rich_sparse_otu_table_hdf5.biom" {
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
          (0): "2014-07-29T16:16:36.617320"
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
          (0): 2, 1
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
          GROUP "group-metadata" {
          }
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
          GROUP "metadata" {
             DATASET "taxonomy" {
                DATATYPE  H5T_STRING {
                   STRSIZE H5T_VARIABLE;
                   STRPAD H5T_STR_NULLTERM;
                   CSET H5T_CSET_ASCII;
                   CTYPE H5T_C_S1;
                }
                DATASPACE  SIMPLE { ( 5, 7 ) / ( 5, 7 ) }
                DATA {
                (0,0): "k__Bacteria", "p__Proteobacteria",
                (0,2): "c__Gammaproteobacteria", "o__Enterobacteriales",
                (0,4): "f__Enterobacteriaceae", "g__Escherichia", "s__",
                (1,0): "k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae",
                (1,3): "o__Nostocales", "f__Nostocaceae", "g__Dolichospermum",
                (1,6): "s__",
                (2,0): "k__Archaea", "p__Euryarchaeota", "c__Methanomicrobia",
                (2,3): "o__Methanosarcinales", "f__Methanosarcinaceae",
                (2,5): "g__Methanosarcina", "s__",
                (3,0): "k__Bacteria", "p__Firmicutes", "c__Clostridia",
                (3,3): "o__Halanaerobiales", "f__Halanaerobiaceae",
                (3,5): "g__Halanaerobium", "s__Halanaerobiumsaccharolyticum",
                (4,0): "k__Bacteria", "p__Proteobacteria",
                (4,2): "c__Gammaproteobacteria", "o__Enterobacteriales",
                (4,4): "f__Enterobacteriaceae", "g__Escherichia", "s__"
                }
             }
          }
       }
       GROUP "sample" {
          GROUP "group-metadata" {
          }
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
          GROUP "metadata" {
             DATASET "BODY_SITE" {
                DATATYPE  H5T_STRING {
                   STRSIZE H5T_VARIABLE;
                   STRPAD H5T_STR_NULLTERM;
                   CSET H5T_CSET_UTF8;
                   CTYPE H5T_C_S1;
                }
                DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
                DATA {
                (0): "gut", "gut", "gut", "skin", "skin", "skin"
                }
             }
             DATASET "BarcodeSequence" {
                DATATYPE  H5T_STRING {
                   STRSIZE H5T_VARIABLE;
                   STRPAD H5T_STR_NULLTERM;
                   CSET H5T_CSET_UTF8;
                   CTYPE H5T_C_S1;
                }
                DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
                DATA {
                (0): "CGCTTATCGAGA", "CATACCAGTAGC", "CTCTCTACCTGT",
                (3): "CTCTCGGCCTGT", "CTCTCTACCAAT", "CTAACTACCAAT"
                }
             }
             DATASET "Description" {
                DATATYPE  H5T_STRING {
                   STRSIZE H5T_VARIABLE;
                   STRPAD H5T_STR_NULLTERM;
                   CSET H5T_CSET_UTF8;
                   CTYPE H5T_C_S1;
                }
                DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
                DATA {
                (0): "human gut", "human gut", "human gut", "human skin",
                (4): "human skin", "human skin"
                }
             }
             DATASET "LinkerPrimerSequence" {
                DATATYPE  H5T_STRING {
                   STRSIZE H5T_VARIABLE;
                   STRPAD H5T_STR_NULLTERM;
                   CSET H5T_CSET_UTF8;
                   CTYPE H5T_C_S1;
                }
                DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
                DATA {
                (0): "CATGCTGCCTCCCGTAGGAGT", "CATGCTGCCTCCCGTAGGAGT",
                (2): "CATGCTGCCTCCCGTAGGAGT", "CATGCTGCCTCCCGTAGGAGT",
                (4): "CATGCTGCCTCCCGTAGGAGT", "CATGCTGCCTCCCGTAGGAGT"
                }
             }
          }
       }
    }
    }
