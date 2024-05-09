#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2011-2020, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import inspect
from contextlib import contextmanager
import io
import codecs
import functools

from collections import defaultdict
from os import getenv
from os.path import abspath, dirname, exists
import re
from hashlib import md5
from gzip import open as gzip_open
import h5py

try:
    H5PY_VLEN_STR = h5py.special_dtype(vlen=str)
    H5PY_VLEN_UNICODE = h5py.special_dtype(vlen=str)

except ImportError:
    H5PY_VLEN_STR = None
    H5PY_VLEN_UNICODE = None

from numpy import mean, median, min, max

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2020, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso",
               "Jose Clemente", "Justin Kuczynski", "Jorge Cañardo Alastuey"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__format_version__ = (2, 1)
__version__ = "2.1.16"


def generate_subsamples(table, n, axis='sample', by_id=False):
    """Indefinitely generate random subsamples

    Parameters
    ----------
    table : Table
        The table to subsample from
    n : int
        The size of the subsample
    axis : {'sample', 'observation'}, optional
        The axis to operate on, defaults to 'sample'.
    by_id : bool, optional
        If `True`, operate on IDs, if `False`, operate on values. Defaults to
        `False`.

    Returns
    -------
    GeneratorType
        Each subsequent call to .next() will yield a new randomly subsampled
        table

    Examples
    --------

    Randomly subsample the samples in the example_table 10 times:

    >>> from biom import example_table
    >>> gen = generate_subsamples(example_table, 2, by_id=True)
    >>> observed_ids = []
    >>> for _, table in zip(range(100), gen):
    ...     observed_ids.append(tuple(table.ids()))
    >>> print(sorted(set(observed_ids)))
    [('S1', 'S2'), ('S1', 'S3'), ('S2', 'S3')]

    """
    while True:
        yield table.subsample(n, axis, by_id)


def get_biom_format_version_string(version=None):
    """Returns the current Biom file format version.

    Parameters
    ----------
    version : tuple
        a tuple containing the version number of the biom table
    """
    if version is None:
        return "Biological Observation Matrix 1.0.0"
    else:
        return f"Biological Observation Matrix {version[0]}.{version[1]}.0"


def get_biom_format_url_string():
    """Returns the current Biom file format description URL."""
    return __url__


def unzip(items):
    """Performs the reverse of zip, i.e. produces separate lists from tuples.

    items should be list of k-element tuples. Will raise exception if any
    tuples contain more items than the first one.

    Conceptually equivalent to transposing the matrix of tuples.

    Returns list of lists in which the ith list contains the ith element of
    each tuple.

    Note: zip expects *items rather than items, such that unzip(zip(*items))
    returns something that compares equal to items.

    Always returns lists: does not check original data type, but will accept
    any sequence.

    This function is ported from PyCogent (http://www.pycogent.org). PyCogent
    is a GPL project, but we obtained permission from the authors of this
    function to port it to the BIOM Format project (and keep it under BIOM's
    BSD license).
    """
    if items:
        return list(map(list, zip(*list(items))))
    else:
        return []


def flatten(items):
    """Removes one level of nesting from items.

    items can be any sequence, but flatten always returns a list.

    This function is ported from PyCogent (http://www.pycogent.org). PyCogent
    is a GPL project, but we obtained permission from the authors of this
    function to port it to the BIOM Format project (and keep it under BIOM's
    BSD license).
    """
    result = []
    for i in items:
        try:
            result.extend(i)
        except:  # noqa
            result.append(i)
    return result


def _natsort_key(item):
    """Provides normalized version of item for sorting with digits.

    This function is ported from QIIME (http://www.qiime.org) and is based on:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    QIIME is a GPL project, but we obtained permission from the authors of this
    function to port it to the BIOM Format project (and keep it under BIOM's
    BSD license).
    """
    item = str(item)
    try:
        chunks = re.split(r'(\d+(?:\.\d+)?)', item)
    except TypeError:
        # if item is a tuple or list (i.e., indexable, but not a string)
        # work with the first element
        chunks = re.split(r'(\d+(?:\.\d+)?)', item[0])
    for ii in range(len(chunks)):
        if chunks[ii] and chunks[ii][0] in '0123456789':
            if '.' in chunks[ii]:
                numtype = float
            else:
                numtype = int
            # wrap in tuple with '0' to explicitly specify numbers come first
            chunks[ii] = (0, numtype(chunks[ii]))
        else:
            chunks[ii] = (1, chunks[ii])
    return (chunks, item)


def natsort(seq):
    """Sort a sequence of text strings in a reasonable order.

    This function is ported from QIIME (http://www.qiime.org) and is based on:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    QIIME is a GPL project, but we obtained permission from the authors of this
    function to port it to the BIOM Format project (and keep it under BIOM's
    BSD license).
    """
    alist = list(seq)
    alist.sort(key=_natsort_key)
    return alist


def prefer_self(x, y):
    """Merge metadata method, return X if X else Y"""
    return x if x is not None else y


def index_list(item):
    """Takes a list and returns {l[idx]:idx}"""
    return {id_: idx for idx, id_ in enumerate(item)}


def load_biom_config():
    """Returns biom-format configuration read in from file.

    This function is ported from QIIME (http://www.qiime.org), previously named
    load_qiime_config. QIIME is a GPL project, but we obtained permission from
    the authors of this function to port it to the BIOM Format project (and
    keep it under BIOM's BSD license).
    """
    biom_config_fps = []
    biom_project_dir = get_biom_project_dir()
    biom_config_fps.append(biom_project_dir + '/support_files/biom_config')

    biom_config_env_fp = getenv('BIOM_CONFIG_FP')
    if biom_config_env_fp:
        biom_config_fps.append(biom_config_env_fp)

    home_dir = getenv('HOME')
    if home_dir:
        biom_config_home_fp = home_dir + '/.biom_config'
        biom_config_fps.append(biom_config_home_fp)

    biom_config_files = []
    for biom_config_fp in biom_config_fps:
        if exists(biom_config_fp):
            biom_config_files.append(open(biom_config_fp))

    return parse_biom_config_files(biom_config_files)


def get_biom_project_dir():
    """Returns the top-level biom-format directory.

    This function is ported from QIIME (http://www.qiime.org), previously named
    get_qiime_project_dir. QIIME is a GPL project, but we obtained permission
    from the authors of this function to port it to the BIOM Format project
    (and keep it under BIOM's BSD license).
    """
    # Get the full path of util.py.
    current_fp = abspath(__file__)

    # Get the directory containing util.py.
    current_dir = dirname(current_fp)

    # Get the directory containing that directory.
    current_dir = dirname(current_dir)

    # Return the directory containing that directory.
    return dirname(current_dir)


def parse_biom_config_files(biom_config_files):
    """Parses files in (ordered!) list of biom_config_files.

    The order of files must be least important to most important. Values
    defined in earlier files will be overwritten if the same values are defined
    in later files.

    This function is ported from QIIME (http://www.qiime.org), previously named
    parse_qiime_config_files. QIIME is a GPL project, but we obtained
    permission from the authors of this function to port it to the BIOM Format
    project (and keep it under BIOM's BSD license).
    """
    # The biom_config object is a default dict: if keys are not present, None
    # is returned.
    def return_none():
        return None

    results = defaultdict(return_none)
    for biom_config_file in biom_config_files:
        try:
            results.update(parse_biom_config_file(biom_config_file))
        except OSError:
            pass

    return results


def parse_biom_config_file(biom_config_file):
    """Parses lines in a biom_config file.

    This function is ported from QIIME (http://www.qiime.org), previously named
    parse_qiime_config_file. QIIME is a GPL project, but we obtained permission
    from the authors of this function to port it to the BIOM Format project
    (and keep it under BIOM's BSD license).
    """
    result = {}
    for line in biom_config_file:
        line = line.strip()

        # Ignore blank lines or lines beginning with '#'.
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        param_id = fields[0]
        param_value = ' '.join(fields[1:]) or None
        result[param_id] = param_value
    return result


def compute_counts_per_sample_stats(table, binary_counts=False):
    """Return summary statistics on per-sample observation counts

        table: a BIOM table object
        binary_counts: count the number of unique observations per
         sample, rather than the sum of the total counts (i.e., counts
         are qualitative rather than quantitative)

    This function is ported from QIIME (http://www.qiime.org), previously named
    compute_seqs_per_library_stats. QIIME is a GPL project, but we obtained
    permission from the authors of this function to port it to the BIOM Format
    project (and keep it under BIOM's BSD license).
    """
    sample_counts = {}
    for count_vector, sample_id, metadata in table.iter():
        if binary_counts:
            sample_counts[sample_id] = (count_vector != 0).sum()
        else:
            sample_counts[sample_id] = float(count_vector.sum())
    counts = list(sample_counts.values())

    if len(counts) == 0:
        return (0, 0, 0, 0, sample_counts)
    else:
        return (min(counts),
                max(counts),
                median(counts),
                mean(counts),
                sample_counts)


def safe_md5(open_file, block_size=2 ** 20):
    """Computes an md5 sum without loading the file into memory

    This method is based on the answers given in:
    http://stackoverflow.com/questions/1131220/get-md5-hash-of-a-files-without-open-it-in-python

    This function is ported from PyCogent (http://www.pycogent.org). PyCogent
    is a GPL project, but we obtained permission from the authors of this
    function to port it to the BIOM Format project (and keep it under BIOM's
    BSD license).
    """
    data = True
    result = md5()

    # While a little hackish, this allows this code to
    # safely work either with a file object or a list of lines.
    if hasattr(open_file, 'read'):
        data_getter = open_file.read
        data_getter_i = block_size
    elif isinstance(open_file, list):
        def f(i):
            try:
                return open_file.pop(i)
            except IndexError:
                return None
        data_getter = f
        data_getter_i = 0
    else:
        raise TypeError("safe_md5 can only handle a file handle or list of "
                        "lines but recieved %r." % type(open_file))

    while data:
        data = data_getter(data_getter_i)
        if data:
            result.update(data.encode('utf-8'))
    return result.hexdigest()


def is_gzip(fp):
    """Checks the first two bytes of the file for the gzip magic number

    If the first two bytes of the file are 1f 8b (the "magic number" of a
    gzip file), return True; otherwise, return false.

    This function is ported from QIIME (http://www.qiime.org). QIIME is a GPL
    project, but we obtained permission from the authors of this function to
    port it to the BIOM Format project (and keep it under BIOM's BSD license).
    """
    with open(fp, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


@contextmanager
def biom_open(fp, permission='r'):
    """Wrapper to allow opening of gzipped or non-compressed files

    Read or write the contents of a file

    Parameters
    ----------
    file_fp : file path or pathlib.Path
    permission : str, {'r', 'w', 'wb', 'rb', 'U'}

    Returns
    -------
    [h5py.File, file, gzip.GzipFile]

    Notes
    -----
    If the file is binary, be sure to pass in a binary mode (append 'b' to
    the mode); opening a binary file in text mode (e.g., in default mode 'U')
    will have unpredictable results.

    If h5py is available on the system, you cannot use biom_open to create a
    writable ASCII file handle. You can use it to create writable GZIP handles
    and HDF5 handles, however.

    This function is ported from QIIME (http://www.qiime.org), previously named
    qiime_open. QIIME is a GPL project, but we obtained permission from the
    authors of this function to port it to the BIOM Format project (and keep it
    under BIOM's BSD license).

    Raises
    ------
    RuntimeError
        If the user tries to parse an HDF5 file without having h5py installed.
    ValueError
        If the user tries to read an empty file.

    """
    if permission not in ['r', 'w', 'U', 'rb', 'wb']:
        raise OSError("Unknown mode: %s" % permission)

    mode = permission

    # don't try to open an HDF5 file if H5PY is not installed, this can only
    # happen if we are reading a file
    if mode in {'r', 'rb', 'U'}:
        if os.path.getsize(fp) == 0:
            raise ValueError("The file '%s' is empty and can't be parsed" % fp)

    if mode in ['U', 'r', 'rb'] and is_gzip(fp):
        def opener(fp, mode):  # noqa
            return codecs.getreader('utf-8')(gzip_open(fp, mode))
        mode = 'rb' if permission in ['U', 'r'] else permission
    elif mode in ['w', 'wb'] and str(fp).endswith('.gz'):
        def opener(fp, mode):  # noqa
            codecs.getwriter('utf-8')(gzip_open(fp, mode))
    elif mode in ['U', 'r', 'rb'] and h5py.is_hdf5(fp):
        opener = h5py.File
        mode = 'r' if permission == 'U' else permission
    elif mode == 'w':
        opener = h5py.File
    else:
        opener = functools.partial(io.open, encoding='utf-8')

    f = opener(fp, mode)
    try:
        yield f
    finally:
        f.close()


def get_data_path(fn):
    """Return path to filename ``fn`` in the data folder.

    During testing it is often necessary to load data files. This
    function returns the full path to files in the ``data`` subfolder.

    Parameters
    ----------
    fn : str
        File name.

    Returns
    -------
    str
        Inferred absolute path to the test data for the module where
        ``get_data_path(fn)`` is called.

    Notes
    -----
    The requested path may not point to an existing file, as its
    existence is not checked.

    This method was adapted from scikit-bio, specifically `skbio.util.testing`.
    """
    # getouterframes returns a list of tuples: the second tuple
    # contains info about the caller, and the second element is its
    # filename
    callers_filename = inspect.getouterframes(inspect.currentframe())[1][1]
    path = os.path.dirname(os.path.abspath(callers_filename))
    data_path = os.path.join(path, 'test_data', fn)
    return data_path


def is_hdf5_file(fp):
    """Guess if file is HDF5.

    Parameters
    ----------
    fn : str
        File name

    Returns
    -------
    bool
        Whether the file is an HDF5 file
    """
    with open(fp, 'rb') as f:
        # from the HDF5 documentation about format signature
        return f.read(8) == b'\x89HDF\r\n\x1a\n'
