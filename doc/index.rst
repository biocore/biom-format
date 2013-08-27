.. BIOM documentation master file

The Biological Observation Matrix (BIOM) format
===============================================

The `BIOM file format <http://www.biom-format.org>`_ (canonically pronounced `biome`) is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for the `Earth Microbiome Project <http://www.earthmicrobiome.org>`_ and is a `Genomics Standards Consortium <http://gensc.org/>`_ candidate project.

The `BIOM format <http://www.biom-format.org>`_ is designed for general use in broad areas of comparative -omics. For example, in marker-gene surveys, the primary use of this format is to represent OTU tables: the observations in this case are OTUs and the matrix contains counts corresponding to the number of times each OTU is observed in each sample. With respect to metagenome data, this format would be used to represent metagenome tables: the observations in this case might correspond to SEED subsystems, and the matrix would contain counts corresponding to the number of times each subsystem is observed in each metagenome. Similarly, with respect to genome data, this format may be used to represent a set of genomes: the observations in this case again might correspond to SEED subsystems, and the counts would correspond to the number of times each subsystem is observed in each genome.

There are two components to the BIOM project: first is `definition of the BIOM format <./documentation/biom_format.html>`_, and second is `development of support objects <./documentation/table_objects.html>`_ in multiple programming languages to support the use of BIOM in diverse bioinformatics applications. The version of the BIOM file format is independent of the version of the `biom-format` software.

There are official implementations of BIOM format support objects (APIs) in the Python and R programming languages. The rest of this site contains details about the BIOM file format (which is independent of the API) and the Python ``biom-format`` API. For more details about the R API, please see the `CRAN biom package <http://cran.r-project.org/web/packages/biom/index.html>`_.

Contents
========

.. toctree::
   :maxdepth: 3

   ./documentation/index.rst
   ./BIOM_LICENSE.rst

BIOM version
============

The latest official version of the biom-format project is |release| and of the BIOM file format is 1.0. Details on the file format can be found `here <./documentation/biom_format.html>`_.

Installing the biom-format project
==================================

To install the ``biom-format`` project, you can download the release version `biom-format-1.2.0 <ftp://thebeast.colorado.edu/pub/biom-format-releases/biom-format-1.2.0.tar.gz>`_, or work with the development version. Generally we recommend working with the release version as it will be more stable, but if you want access to the latest features (and can tolerate some instability) you should work with the development version.

The biom-format project has the following dependencies:
	* Python >= 2.7 and < 3.0
	* numpy >= 1.3.0
	* `pyqi <http://bipy.github.io/pyqi>`_ 0.2.0
	* gcc >= 4.1.2 (optional; used for more efficient sparse table representations)
	* cython >= 0.14.1 (optional; used for more efficient sparse table representations)
	* dateutil (optional; must be installed if using the ``biom validate-table`` command)

We'll illustrate the install process in the ``$HOME/code`` directory. You can either work in this directory on your system (creating it, if necessary, by running ``mkdir $HOME/code``) or replace all occurrences of ``$HOME/code`` in the following instructions with your working directory. Change to this directory to start the install process::

	cd $HOME/code

To install the release version, download from `biom-format-1.2.0 <ftp://thebeast.colorado.edu/pub/biom-format-releases/biom-format-1.2.0.tar.gz>`_, uncompress the file, and change to the resulting directory::

	wget ftp://thebeast.colorado.edu/pub/biom-format-releases/biom-format-1.2.0.tar.gz
	tar -xvzf biom-format-1.2.0.tar.gz
	cd $HOME/code/biom-format-1.2.0

Alternatively, to install the development version, pull it from github, and change to the resulting directory::

	git clone git://github.com/biom-format/biom-format.git
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

The ``biom`` command referenced in the previous section is a driver for commands in biom-format, powered by `pyqi <http://bipy.github.io/pyqi>`_. You can enable tab completion of biom command names and command options (meaning that when you begin typing the name of a command or option you can auto-complete it by hitting the *tab* key) by following a few simple steps from the pyqi documentation. While this step is optional, tab completion is very convenient so it's worth enabling.

To enable tab completion, follow the steps outlined under `Configuring bash completion <http://bipy.github.io/pyqi/doc/tutorials/defining_your_command_driver.html#configuring-bash-completion>`_ in the pyqi install documentation, substituting ``biom`` for ``my-project``/``my_project`` in all commands. After completing those steps and closing and re-opening your terminal, auto-completion should be enabled.

Citing the BIOM project
=======================

You can cite the BIOM format as follows:

| The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome.
| Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso.
| GigaScience, June 2012.

Development team
================

The biom-format project was conceived of and developed by the `QIIME <http://www.qiime.org>`_, `MG-RAST <http://metagenomics.anl.gov>`_, and `VAMPS <http://vamps.mbl.edu/>`_ development groups to support interoperability of our software packages. If you have questions about the biom-format project you can contact gregcaporaso@gmail.com.


