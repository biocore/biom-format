#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Antonio Gonzalez",
               "Yoshiki Vazquez Baeza", "Jose Navas", "Adam Robbins-Pianka",
               "Rob Knight", "Joshua Shorenstein", "Emily TerAvest"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__version__ = "1.3.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

__all__ = ['table', 'parse', 'util', 'unit_test', 'exception']

from sys import modules, stderr

from biom.exception import InvalidSparseBackendException
from biom.util import load_biom_config

biom_config = load_biom_config()

sparse_backends = ['CSMat', 'scipy_sparse_mat']


def set_sparse_backend(sparse_backend, warn=True):
    """Sets the sparse matrix backend to use in biom.table.

    Use this function to programmatically override the backend specified in the
    biom config file. Similar to matplotlib.use(). Has no effect if biom.table
    has already been imported.

    Arguments:
        sparse_backend - sparse backend string identifier
        warn - if True, will print out a warning if biom.table has already been
        loaded.
    """
    if 'biom.table' in modules:
        if warn:
            print ("Warning: biom.table has already been loaded. This call to "
                   "biom.set_sparse_backend() has no effect. It must be "
                   "called before biom.table is imported for the first time.")
    else:
        biom_config['python_code_sparse_backend'] = sparse_backend


def get_sparse_backend():
    """Returns the constructor and functions needed to use the current backend.

    Will look at whatever the current backend is in the loaded biom config
    dict. If one isn't specified, will default to scipy_sparse_mat. Will raise
    an InvalidSparseBackendException if the current sparse backend isn't
    supported or cannot be used for whatever reason.
    """
    backend = biom_config['python_code_sparse_backend']
    if backend is None:
        backend = 'scipy_sparse_mat'

    if backend not in sparse_backends:
        raise InvalidSparseBackendException("Unrecognized sparse backend "
                                            "'%s'. Choose from %s." %
                                            (backend,
                                             ', '.join(sparse_backends)))

    valid_backend = False
    if backend == 'scipy_sparse_mat':
        try:
            from biom.backends.scipysparse import ScipySparseMat, to_scipy, \
                dict_to_scipy, list_dict_to_scipy, list_nparray_to_scipy, \
                nparray_to_scipy, list_list_to_scipy
            sparse_obj = ScipySparseMat
            to_sparse = to_scipy
            dict_to_sparse_obj = dict_to_scipy
            list_dict_to_sparse_obj = list_dict_to_scipy
            list_nparray_to_sparse_obj = list_nparray_to_scipy
            nparray_to_sparse_obj = nparray_to_scipy
            list_list_to_sparse_obj = list_list_to_scipy
            valid_backend = True
        except ImportError:
            valid_backend = False
            stderr.write("Cannot load scipy_sparse_mat (requires that scipy "
                         "is installed). Using CSMat sparse backend.\n")

    if backend == 'CSMat' or (not valid_backend):
        try:
            from biom.backends.csmat import CSMat, to_csmat, dict_to_csmat, \
                list_dict_to_csmat, list_nparray_to_csmat, nparray_to_csmat, \
                list_list_to_csmat
            sparse_obj = CSMat
            to_sparse = to_csmat
            dict_to_sparse_obj = dict_to_csmat
            list_dict_to_sparse_obj = list_dict_to_csmat
            list_nparray_to_sparse_obj = list_nparray_to_csmat
            nparray_to_sparse_obj = nparray_to_csmat
            list_list_to_sparse_obj = list_list_to_csmat
            valid_backend = True
        except ImportError:
            valid_backend = False
            stderr.write('Cannot load CSMat sparse backend.\n')

    if not valid_backend:
        raise InvalidSparseBackendException("The sparse matrix backend '%s' "
                                            "could not be loaded. Please "
                                            "check your biom-format "
                                            "installation." % backend)

    return sparse_obj, to_sparse, dict_to_sparse_obj, list_dict_to_sparse_obj, \
        list_nparray_to_sparse_obj, nparray_to_sparse_obj, \
        list_list_to_sparse_obj
