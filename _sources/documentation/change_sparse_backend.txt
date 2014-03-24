.. _change_sparse_backend:

======================
Sparse matrix backends
======================

As the BIOM project evolves, so do the underlying data structures, leading to potential runtime and memory trade-offs between implementations. The BIOM project currently supports two sparse matrix backends: ``CSMat`` and ``ScipySparseMat``. ``CSMat`` is the default sparse matrix backend used in the BIOM project.

How to check what sparse matrix backend is in use
=================================================

To check what sparse matrix backend is in use, simply execute ``biom show-install-info``. The last line shows the SparseObj type. For instance::

 $ biom show-install-info

 System information
 ==================
           Platform:	darwin
 Python/GCC version:	2.7.1 (r271:86832, Aug 30 2012, 10:07:33)  [GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)]
  Python executable:	/Users/caporaso/.virtualenvs/biom/bin/python
 
 Dependency versions
 ===================
     pyqi version:	0.2.0
    NumPy version:	1.7.1
 dateutil version:	1.5
 
 biom-format package information
 ===============================
 biom-format version:	1.2.0
      SparseObj type:	CSMat

The last line indicates that the ``CSMat`` sparse matrix backend is in use.

Changing the sparse backend
===========================

There are two methods to change the backend that is used. The first method is by copying the ``biom_config`` file located under ``support_files/`` and placing it in your home directory as ``~/.biom_config``. Then, edit ``~/.biom_config`` and replace the current backend type with the desired type.

The second method is to set the environment variable ``$BIOM_CONFIG_FP`` to a file path of your choice, and place the following line into that file::

	python_code_sparse_backend	<BACKEND TYPE>

Where ``<BACKEND TYPE>`` is replaced by the specific backend implementation to use.

Sparse matrix backend descriptions
==================================

Different sparse matrix backends have different performance characteristics. As the BIOM project changes over time, additional backends may be added that address specific performance concerns.

CSMat
-----

The default sparse matrix backend is ``CSMat``. ``CSMat`` is a pure-Python implementation of coordinate list, compressed sparse row and compressed sparse column formats and facilitates interaction with these representations. ``CSMat`` is recommended if you are unable to install `scipy <http://www.scipy.org/>`_, which is a required dependency when using ``ScipySparseMat``.

ScipySparseMat
--------------

``ScipySparseMat`` is a relatively lightweight Python wrapper around `scipy's <http://www.scipy.org/>`_ sparse matrix library (``scipy.sparse``). You will need to have scipy installed if you plan to use this backend. In general, expect to see performance improvements in both runtime and overall memory consumption when using this backend, especially with larger sparse BIOM tables (e.g., thousands of samples and/or observations).
