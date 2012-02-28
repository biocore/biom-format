.. BIOM documentation master file

The Biological Observation Matrix (BIOM) format
===============================================

The `BIOM file format <http://www.biom-format.org>`_ (canonically pronounced `biome`) is designed to be a general-use format for representing biological sample by observation contingency tables. We are currently seeking recognition of the BIOM format by the `Genomics Standards Consortium <http://gensc.org/>`_.

The `BIOM format <http://www.biom-format.org>`_ is designed for general use in broad areas of comparative -omics. For example, in marker-gene surveys, the primary use of this format is to represent OTU tables: the observations in this case are OTUs and the matrix contains counts corresponding to the number of times each OTU is observed in each sample. With respect to metagenome data, this format would be used to represent metagenome tables: the observations in this case might correspond to SEED subsystems, and the matrix would contain counts corresponding to the number of times each subsystem is observed in each metagenome. Similarly, with respect to genome data, this format may be used to represent a set of genomes: the observations in this case again might correspond to SEED subsystems, and the counts would correspond to the number of times each subsystem is observed in each genome.

There are two components to the BIOM project: first is `definition of the BIOM format <./documentation/biom_format.html>`_, and second is `development of support objects <./documentation/table_objects.html>`_ in multiple programming languages to support the use of BIOM in diverse bioinformatics applications. A 1.0 release of the biom-format project will coincide with publication of a manuscript describing the project.

Contents
========

.. toctree::
   :maxdepth: 3

   ./documentation/index.rst

BIOM version
============

The latest official version of the biom-format project and file format is |release|. 

The biom-format project and file format version are always the same. Official versions contain three integers in the following format: ``major-version.incremental-version.minor-version``. When ``-dev`` is appended to the end of a version string that indicates a development (or between-release version). For example, ``1.0.0-dev`` would refer to the development version following the 1.0.0 official release. 

Citing the BIOM project
=======================

While our manuscript is under review you can cite the BIOM project with this URL: `http://www.biom-format.org <http://www.biom-format.org>`_.

Development team
================

The biom-format project was conceived of and developed by the `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_ development groups to support interoperability of our software packages. If you have questions about the biom-format project you can contact gregcaporaso@gmail.com.


