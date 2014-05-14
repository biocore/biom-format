.. _table_objects:

===========================================
biom-format ``Table`` objects
===========================================

The biom-format project provides rich ``Table`` objects to support use of the BIOM file format. The objects encapsulate matrix data (such as OTU counts) and abstract the interaction away from the programmer.

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

``Table`` is a container object and an abstract base class that provides a common and required API for subclassed objects. Through the use of private interfaces, it is possible to create public methods that operate on the underlying datatype without having to implement each method in each subclass. For instance, ``Table.iter_data`` will return a generator that always yields ``numpy.array`` vectors for each sample regardless of how the table data is actually stored. This functionality results from derived classes implementing private interfaces, such as ``Table._to_dense``.

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
