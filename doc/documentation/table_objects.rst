.. _table_objects:

===========================================
biom-format ``Table`` objects
===========================================

The biom-format project provides rich ``Table`` objects to support use of the BIOM file format. The objects encapsulate matrix data (such as OTU counts) and abstract the interaction away from the programmer. This provides the immediate benefit of the programmer not having to worry about what the underlying data object is, and in turn allows for different data representations to be supported. Currently, biom-format supports a ``dense`` object built off of ``numpy.array`` (`NumPy <http://numpy.scipy.org/>`_) and a ``sparse`` object built off of Python dictionaries. 

biom-format ``table_factory`` method
=========================================

Generally, construction of a ``Table`` subclass will be through the ``table_factory`` method. This method facilitates any necessary data conversions and supports a wide variety of input data types.

.. autofunction:: biom.table.table_factory

Description of available ``Table`` objects
==========================================

There are multiple objects available but some of them are unofficial abstract base classes (does not use the ``abc`` module for historical reasons). In practice, the objects used should be the derived Tables such as ``SparseOTUTable`` or ``DenseGeneTable``. 

Abstract base classes
---------------------

Abstract base classes establish standard interfaces for subclassed types and provide common functionality for derived types. 

``Table``
^^^^^^^^^

``Table`` is a container object and an abstract base class that provides a common and required API for subclassed objects. Through the use of private interfaces, it is possible to create public methods that operate on the underlying datatype without having to implement each method in each subclass. For instance, ``Table.iter_samplesData`` will return a generator that always yields ``numpy.array`` vectors for each sample regardless of how the table data is actually stored. This functionality results from derived classes implementing private interfaces, such as ``Table._conv_to_np``.

.. autoclass:: biom.table.Table
   :members:

``OTUTable``
^^^^^^^^^^^^

The ``OTUTable`` base class provides functionality specific for OTU tables. Currently, it only provides a static private member variable that describes its ``BIOM`` type. This object was stubbed out incase future methods are developed that do not make sense with the context of, say, an MG-RAST metagenomic abundance table. It is advised to always use an object that subclasses ``OTUTable`` if the analysis is on OTU data.

.. autoclass:: biom.table.OTUTable
   :members:

``PathwayTable``
^^^^^^^^^^^^^^^^

A table type to represent gene pathways.

.. autoclass:: biom.table.PathwayTable
   :members:

``FunctionTable``
^^^^^^^^^^^^^^^^^

A table type to represent gene functions.

.. autoclass:: biom.table.FunctionTable
   :members:

``OrthologTable``
^^^^^^^^^^^^^^^^^

A table type to represent gene orthologs.

.. autoclass:: biom.table.OrthologTable
   :members:

``GeneTable``
^^^^^^^^^^^^^

A table type to represent genes.

.. autoclass:: biom.table.GeneTable
   :members:

``MetaboliteTable``
^^^^^^^^^^^^^^^^^^^

A table type to represent metabolite profiles.

.. autoclass:: biom.table.MetaboliteTable
   :members:

``TaxonTable``
^^^^^^^^^^^^^^

A table type to represent taxonomies. 

.. autoclass:: biom.table.TaxonTable
   :members:

Container classes
-----------------
The container classes implement required private member variable interfaces as defined by the ``Table`` abstract base class. Specifically, these objects define the ways in which data is moved into and out of the contained data object. These are fully functional and usable objects, however they do not implement table type specifc functionality.

``SparseTable``
^^^^^^^^^^^^^^^

The subclass ``SparseTable`` can be derived for use with table data. This object implemented all of the required private interfaces specified by the ``Table`` base class. The object contains a ``_data`` private member variable that is an instance of the current sparse backend. It is advised to used derived objects of SparseTable if the data being operated on is sparse.

.. autoclass:: biom.table.SparseTable
   :members:

``DenseTable``
^^^^^^^^^^^^^^

The ``DenseTable`` object fulfills all private member methods stubbed out by the ``Table`` base class. The dense table contains a private member variable that is an instance of ``numpy.array``. The ``array`` object is a matrix that contains all values including zeros. It is advised to use this table only if the number of samples and observations is reasonble. Unfortunately, it isn't reasonable to define reasonable in this context. However, if either the number of observations or the number of samples is > 1000, it would probably be a good idea to rely on a ``SparseTable``.

.. autoclass:: biom.table.DenseTable
   :members:

Table type objects
------------------

The table type objects define variables and methods specific to a table type. These inherit from a ``Container Class`` and a table type base class, and are therefore instantiable. Generally you'll instantiate tables with ``biom.table.table_factory``, and one of these will be passed as the ``constructor`` argument.

``DenseOTUTable``
^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.DenseOTUTable
   :members:

``SparseOTUTable``
^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.SparseOTUTable
   :members:

``DensePathwayTable``
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.DensePathwayTable
   :members:

``SparsePathwayTable``
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.SparsePathwayTable
   :members:

``DenseFunctionTable``
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.DenseFunctionTable
   :members:

``SparseFunctionable``
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.SparseFunctionTable
   :members:

``DenseOrthologTable``
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.DenseOrthologTable
   :members:

``SparseOrthologTable``
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.SparseOrthologTable
   :members:

``DenseGeneTable``
^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.DenseGeneTable
   :members:

``SparseGeneTable``
^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.SparseGeneTable
   :members:

``DenseMetaboliteTable``
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.DenseMetaboliteTable
   :members:

``SparseMetaboliteTable``
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.SparseMetaboliteTable
   :members:
