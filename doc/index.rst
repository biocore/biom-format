.. BIOM documentation master file

.. image:: _static/biom-format.png

The Biological Observation Matrix (BIOM) format
===============================================

The `BIOM file format <http://www.biom-format.org>`_ (canonically pronounced `biome`) is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for the `Earth Microbiome Project <http://www.earthmicrobiome.org>`_ and is a `Genomics Standards Consortium <http://gensc.org/>`_ supported project.

The `BIOM format <http://www.biom-format.org>`_ is designed for general use in broad areas of comparative -omics. For example, in marker-gene surveys, the primary use of this format is to represent OTU tables: the observations in this case are OTUs and the matrix contains counts corresponding to the number of times each OTU is observed in each sample. With respect to metagenome data, this format would be used to represent metagenome tables: the observations in this case might correspond to SEED subsystems, and the matrix would contain counts corresponding to the number of times each subsystem is observed in each metagenome. Similarly, with respect to genome data, this format may be used to represent a set of genomes: the observations in this case again might correspond to SEED subsystems, and the counts would correspond to the number of times each subsystem is observed in each genome.

The BIOM project consists of the following components:

* `definition of the BIOM file format <./documentation/biom_format.html>`_;
* command line interface (CLI) for working with BIOM files, including `converting between file formats <./documentation/biom_conversion.html>`_, `adding metadata to BIOM files <./documentation/adding_metadata.html>`_, and `summarizing BIOM files <./documentation/summarizing_biom_tables.html>`_ (run ``biom`` to see the full list of commands);
* application programming interface (API) for working with BIOM files in multiple programming languages (including Python and R).

The ``biom-format`` package provides a command line interface and Python API for working with BIOM files. The rest of this site contains details about the BIOM file format (which is independent of the API) and the Python ``biom-format`` package. For more details about the R API, please see the `bioconductor biomformat package <https://bioconductor.org/packages/release/bioc/html/biomformat.html>`_.

Projects using the BIOM format
==============================

* `QIIME <http://www.qiime.org>`_
* `MG-RAST <http://metagenomics.anl.gov>`_
* `PICRUSt <http://picrust.github.io/picrust>`_
* `Mothur <http://www.mothur.org/wiki/Make.biom>`_
* `phyloseq <http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html>`_
* `MEGAN <http://ab.inf.uni-tuebingen.de/software/megan5/>`_
* `VAMPS <http://vamps.mbl.edu/>`_
* `metagenomeSeq <http://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html>`_
* `Phinch <http://phinch.org>`_
* `RDP Classifier <https://github.com/rdpstaff/classifier>`_
* `USEARCH <http://www.drive5.com/usearch/>`_
* `PhyloToAST <http://phylotoast.org/>`_
* `EBI Metagenomics <https://www.ebi.ac.uk/metagenomics>`_
* `GCModeller <http://gcmodeller.org>`_
* `MetaPhlAn 2 <http://segatalab.cibio.unitn.it/tools/metaphlan2/>`__

If you are using BIOM in your project, and would like your project to be listed, please submit a `pull request <https://github.com/biocore/biom-format/pulls>`_ to the BIOM project. More information on `submitting pull requests can be found here <https://help.github.com/articles/using-pull-requests>`_.

Contents
========

.. toctree::
   :maxdepth: 3

   ./documentation/index.rst
   ./BIOM_LICENSE.rst

BIOM version
============

The latest official version of the biom-format project is |release| and of the BIOM file format is 2.0. Details on the `file format can be found here <./documentation/biom_format.html>`_.

Installing the ``biom-format`` Python package
=============================================

To install the latest release of the ``biom-format`` Python package::

    pip install numpy
    pip install biom-format

To work with BIOM 2.0+ files::

    pip install h5py

To see a list of all ``biom`` commands, run::

    biom 

To enable Bash tab completion of ``biom`` commands, add the following line to ``$HOME/.bashrc`` (if on Linux) or ``$HOME/.bash_profile`` (if on Mac OS X)::

    eval "$(_BIOM_COMPLETE=source biom)"

Citing the BIOM project
=======================

You can cite the BIOM format as follows (`link <http://www.gigasciencejournal.com/content/1/1/7>`_):

| The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome.
| Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso.
| GigaScience 2012, 1:7. doi:10.1186/2047-217X-1-7

Development team
================

The biom-format project was conceived of and developed by the `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_ development groups to support interoperability of our software packages. If you have questions about the biom-format project please post them on the `QIIME Forum <http://forum.qiime.org>`_.
