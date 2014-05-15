.. _table_objects:

===========================================
biom-format ``Table`` objects
===========================================

The biom-format project provides rich ``Table`` objects to support use of the BIOM file format. The objects encapsulate matrix data (such as OTU counts) and abstract the interaction away from the programmer.

Description of available ``Table`` objects
==========================================

There are multiple objects available but some of them are unofficial abstract base classes (does not use the ``abc`` module for historical reasons). In practice, the objects used should be the derived Tables such as ``SparseOTUTable`` or ``DenseGeneTable``. 

Abstract base classes
---------------------

Abstract base classes establish standard interfaces for subclassed types and provide common functionality for derived types. 

``Table``
^^^^^^^^^

``Table`` is a container object and an abstract base class that provides a common and required API for subclassed objects. Through the use of private interfaces, it is possible to create public methods that operate on the underlying datatype without having to implement each method in each subclass. For instance, ``Table.iter_data`` will return a generator that always yields ``numpy.array`` vectors for each sample regardless of how the table data is actually stored. This functionality results from derived classes implementing private interfaces, such as ``Table._to_dense``.

``Table``
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: biom.table.Table
   :members:

