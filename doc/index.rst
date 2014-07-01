.. BIOM documentation master file

.. image:: _static/biom-format.png

The Biological Observation Matrix (BIOM) format
===============================================

The `BIOM file format <http://www.biom-format.org>`_ (canonically pronounced `biome`) is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for the `Earth Microbiome Project <http://www.earthmicrobiome.org>`_ and is a `Genomics Standards Consortium <http://gensc.org/>`_ supported project.

The `BIOM format <http://www.biom-format.org>`_ is designed for general use in broad areas of comparative -omics. For example, in marker-gene surveys, the primary use of this format is to represent OTU tables: the observations in this case are OTUs and the matrix contains counts corresponding to the number of times each OTU is observed in each sample. With respect to metagenome data, this format would be used to represent metagenome tables: the observations in this case might correspond to SEED subsystems, and the matrix would contain counts corresponding to the number of times each subsystem is observed in each metagenome. Similarly, with respect to genome data, this format may be used to represent a set of genomes: the observations in this case again might correspond to SEED subsystems, and the counts would correspond to the number of times each subsystem is observed in each genome.

There are two components to the BIOM project: first is the `definition of the BIOM format <./documentation/biom_format.html>`_, and second is `development of support objects <./documentation/table_objects.html>`_ in multiple programming languages to support the use of BIOM in diverse bioinformatics applications. The version of the BIOM file format is independent of the version of the `biom-format` software.

There are official implementations of BIOM format support objects (APIs) in the Python and R programming languages. The rest of this site contains details about the BIOM file format (which is independent of the API) and the Python ``biom-format`` API. For more details about the R API, please see the `CRAN biom package <http://cran.r-project.org/web/packages/biom/index.html>`_.

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

Installing the biom-format project
==================================

To install the ``biom-format`` project, you can download the `latest version here <https://pypi.python.org/pypi/biom-format/>`_, or work with the development version. Generally we recommend working with the release version as it will be more stable, but if you want access to the latest features (and can tolerate some instability) you should work with the development version.

The biom-format project has the following dependencies:

	*  `Python <http://www.python.org/>`_ >= 2.7 and < 3.0
	* `numpy <http://www.numpy.org/>`_ >= 1.7.0
	* `pyqi <http://pyqi.readthedocs.org>`_ 0.3.2
	* `scipy <http://www.scipy.org/>`_ >= 0.13.0 
	* `h5py <http://www.h5py.org/>`_ >= 2.20.0 (optional; must be installed if creating or reading HDF5 formatted files)

The easiest way to install the latest version of the biom-format project and its required dependencies is via pip::

	pip install numpy
	pip install biom-format

That's it!

If you decided not to install biom-format using pip, it is also possible to manually install the latest release. We'll illustrate the install process in the ``$HOME/code`` directory. You can either work in this directory on your system (creating it, if necessary, by running ``mkdir $HOME/code``) or replace all occurrences of ``$HOME/code`` in the following instructions with your working directory. Please note that ``numpy`` must be in your installed prior to installing ``biom-format``. Change to this directory to start the install process::

	cd $HOME/code

Download the `latest release, which can be found here <https://pypi.python.org/pypi/biom-format>`_. After downloading, unpack and install (note: x.y.z refers to the downloaded version)::

	tar xzf biom-format-x.y.z.tar.gz
	cd $HOME/code/biom-format-x.y.z

Alternatively, to install the development version, pull it from GitHub, and change to the resulting directory::

	git clone git://github.com/biocore/biom-format.git
	cd $HOME/code/biom-format

To install (either the development or release version), follow these steps::

	sudo python setup.py install

If you do not have sudo access on your system (or don't want to install the ``biom-format`` project in the default location) you'll need to install the library code and scripts in specified directories, and then tell your system where to look for those files. You can do this as follows::

	echo "export PATH=$HOME/bin/:$PATH" >> $HOME/.bashrc
	echo "export PYTHONPATH=$HOME/lib/:$PYTHONPATH" >> $HOME/.bashrc
	mkdir -p $HOME/bin $HOME/lib/
	source $HOME/.bashrc
	python setup.py install --install-scripts=$HOME/bin/ --install-purelib=$HOME/lib/ --install-lib=$HOME/lib/

You should then have access to the biom-format project. You can test this by running the following command::
	
	python -c "from biom import __version__; print __version__"

You should see the current version of the biom-format project.

Next you can run::

	which biom

You should get a file path ending with ``biom`` printed to your screen if it is installed correctly. Finally, to see a list of all ``biom`` commands, run::

	biom

Enabling tab completion of biom commands
----------------------------------------

The ``biom`` command referenced in the previous section is a driver for commands in biom-format, powered by `the pyqi project <http://biocore.github.io/pyqi>`_. You can enable tab completion of biom command names and command options (meaning that when you begin typing the name of a command or option you can auto-complete it by hitting the *tab* key) by following a few simple steps from the pyqi documentation. While this step is optional, tab completion is very convenient so it's worth enabling.

To enable tab completion, follow the steps outlined under `Configuring bash completion <http://biocore.github.io/pyqi/doc/tutorials/defining_your_command_driver.html#configuring-bash-completion>`_ in the pyqi install documentation, substituting ``biom`` for ``my-project`` and ``my_project`` in all commands. After completing those steps and closing and re-opening your terminal, auto-completion should be enabled.

BIOM format in R
================

There is also a BIOM format package for R, called ``biom``. This package includes basic tools for reading biom-format files, accessing and subsetting data tables from a biom object, as well as limited support for writing a biom-object back to a biom-format file. The design of this API is intended to match the python API and other tools included with the biom-format project, but with a decidedly "R flavor" that should be familiar to R users. This includes S4 classes and methods, as well as extensions of common core functions/methods.

To install the latest stable release of the ``biom`` package enter the following command from within an R session::

	install.packages("biom")

To install the latest development version of the ``biom`` package, enter the following lines in an R session::

	install.packages("devtools") # if not already installed
	library("devtools")
	install_github("biom", "joey711")

Please post any support or feature requests and bugs to `the biom issue tracker <https://github.com/joey711/biom/issues>`_.

See `the biom project on GitHub <https://github.com/joey711/biom/>`_ for further details, or if you would like to contribute.

Note that the licenses between the ``biom`` R package (GPL-2) and the other biom-format software (Modified BSD) are different.

Citing the BIOM project
=======================

You can cite the BIOM format as follows (`link <http://www.gigasciencejournal.com/content/1/1/7>`_):

| The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome.
| Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso.
| GigaScience 2012, 1:7. doi:10.1186/2047-217X-1-7

Development team
================

The biom-format project was conceived of and developed by the `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_ development groups to support interoperability of our software packages. If you have questions about the biom-format project you can contact gregcaporaso@gmail.com.
