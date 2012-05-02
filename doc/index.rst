.. BIOM documentation master file

The Biological Observation Matrix (BIOM) format
===============================================

The `BIOM file format <http://www.biom-format.org>`_ (canonically pronounced `biome`) is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for the `Earth Microbiome Project <http://www.earthmicrobiome.org>`_ and is a `Genomics Standards Consortium <http://gensc.org/>`_ candidate project.

The `BIOM format <http://www.biom-format.org>`_ is designed for general use in broad areas of comparative -omics. For example, in marker-gene surveys, the primary use of this format is to represent OTU tables: the observations in this case are OTUs and the matrix contains counts corresponding to the number of times each OTU is observed in each sample. With respect to metagenome data, this format would be used to represent metagenome tables: the observations in this case might correspond to SEED subsystems, and the matrix would contain counts corresponding to the number of times each subsystem is observed in each metagenome. Similarly, with respect to genome data, this format may be used to represent a set of genomes: the observations in this case again might correspond to SEED subsystems, and the counts would correspond to the number of times each subsystem is observed in each genome.

There are two components to the BIOM project: first is `definition of the BIOM format <./documentation/biom_format.html>`_, and second is `development of support objects <./documentation/table_objects.html>`_ in multiple programming languages to support the use of BIOM in diverse bioinformatics applications. A 1.0.0 release of the biom-format project will coincide with publication of a manuscript describing the project.

Contents
========

.. toctree::
   :maxdepth: 3

   ./documentation/index.rst

BIOM version
============

The latest official version of the biom-format project and file format is |release|. 

The biom-format project and file format version are currently the same, but these will be decoupled with the 1.0.0 release. Release versions contain three integers in the following format: ``major-version.incremental-version.minor-version``. When ``-dev`` is appended to the end of a version string that indicates a development (or between-release version). For example, ``1.0.0-dev`` would refer to the development version following the 1.0.0 release. 

Installing the BIOM project
===========================

To install the BIOM project, you can download the release version `biom-format-0.9.3 <https://github.com/downloads/biom-format/biom-format/biom-format-0.9.3.tgz>`_, or work with the development version. Generally we recommend working with the release version as it will be more stable, but if you want access to the latest features (and can tolerate some instability) you should work with the development version. 

To pull the development version from our svn repository, first ``cd`` to the directory where you'd like to install the code. We'll call this ``$HOME/code``:: 

	cd $HOME/code

To install the release version, download from `biom-format-0.9.3 <https://github.com/downloads/biom-format/biom-format/biom-format-0.9.3.tgz>`_, and then run ``tar -xvzf`` on the resulting file to unzip it. To install the development version, run the following command::

	git clone git://github.com/biom-format/biom-format.git
	
To install (either the development or release version), follow these steps::

	cd $HOME/code/biom-format
	sudo python setup.py install

You should then have access to the biom-format project. You can test this by running the following command::
	
	python -c "from biom import __version__; print __version__"

You should see the current version of the biom-format project.

Next you can run::

	which convert_biom.py

You should get ``$HOME/code/biom-format/scripts/convert_biom.py`` printed to your screen if it is installed correctly.


Citing the BIOM project
=======================

While our manuscript is under review you can cite the BIOM project with this URL: `http://www.biom-format.org <http://www.biom-format.org>`_.

Development team
================

The biom-format project was conceived of and developed by the `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_ development groups to support interoperability of our software packages. If you have questions about the biom-format project you can contact gregcaporaso@gmail.com.


