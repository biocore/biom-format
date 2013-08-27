.. _biom_format:

===========================================
The biom file format
===========================================

The BIOM project consists of two independent tools: the `biom-format` software package, which contains software tools for working with BIOM-formatted files and the tables they represent; and the BIOM file format. As of the 1.0.0 software version and the 1.0 file format version, the version of the software and the file format are independent of one another. Version specific documentation of the file formats can be found on the following pages.

.. toctree::
   :maxdepth: 2

   ./format_versions/biom-1.0.rst

Release versions contain three integers in the following format: ``major-version.incremental-version.minor-version``. When ``-dev`` is appended to the end of a version string that indicates a development (or between-release version). For example, ``1.0.0-dev`` would refer to the development version following the 1.0.0 release. 


.. _sparse-or-dense:

============================================
Tips and FAQs regarding the BIOM file format
============================================

Should I generate sparse or dense biom files?
=============================================

In general, we recommend using the sparse format for your biom files. These will be a lot smaller than the dense format biom files when your data is sparse (i.e., more than 85% of your counts are zero). This is common for OTU tables and metagenome tables, and you'll want to investigate whether it's true for your data. If you currently format your data in tab-separated tables where observations are rows and samples are columns, you can format that file to be convertible to biom format with the ``biom convert`` command. Here you can create dense and sparse formats, and see which file size is smaller. See the section on :ref:`converting`. 

Motivation for the BIOM format
==============================

The BIOM format was motivation by several goals. First, to facilitate efficient handling and storage of large, sparse biological contingency tables; second, to support encapsulation of core study data (contingency table data and sample/observation metadata) in a single file; and third, to facilitate the use of these tables between tools that support this format (e.g., passing of data between `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_.).

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


File extension
==============
We recommend that BIOM files use the ``.biom`` extension.
