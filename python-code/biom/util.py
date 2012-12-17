#!/usr/bin/env python

from collections import defaultdict
from os import getenv
from os.path import abspath, dirname, exists
import re

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout", "Greg Caporaso", 
               "Jose Clemente", "Justin Kuczynski"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

def get_biom_format_version_string():
    """Returns the current Biom file format version."""
    return "Biological Observation Matrix 1.0.0"
 
def get_biom_format_url_string():
    """Returns the current Biom file format description URL."""
    return __url__

def unzip(items):
    """Performs the reverse of zip, i.e. produces separate lists from tuples.

    items should be list of k-element tuples. Will raise exception if any tuples
    contain more items than the first one.

    Conceptually equivalent to transposing the matrix of tuples.

    Returns list of lists in which the ith list contains the ith element of each
    tuple.

    Note: zip expects *items rather than items, such that unzip(zip(*items))
    returns something that compares equal to items.

    Always returns lists: does not check original data type, but will accept
    any sequence.

    Method pulled from PyCogent (http://pycogent.sourceforge.net)
    """
    if items:
        return map(list, zip(*items))
    else:
        return []

def flatten(items):
    """Removes one level of nesting from items.

    items can be any sequence, but flatten always returns a list.

    Method pulled from PyCogent (http://pycogent.sourceforge.net) 
    """
    result = [] 
    for i in items:
        try: 
            result.extend(i)
        except:
            result.append(i)
    return result

def _natsort_key(item):
    """Provides normalized version of item for sorting with digits.

    Method pulled from QIIME (http://qiime.org), based on:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    item = str(item)
    try:
        chunks = re.split('(\d+(?:\.\d+)?)', item)
    except TypeError:
        # if item is a tuple or list (i.e., indexable, but not a string)
        # work with the first element
        chunks = re.split('(\d+(?:\.\d+)?)', item[0])
    for ii in range(len(chunks)):
        if chunks[ii] and chunks[ii][0] in '0123456789':
            if '.' in chunks[ii]: numtype = float
            else: numtype = int 
            # wrap in tuple with '0' to explicitly specify numbers come first
            chunks[ii] = (0, numtype(chunks[ii]))
        else:
            chunks[ii] = (1, chunks[ii])
    return (chunks, item)

def natsort(seq):
    """Sort a sequence of text strings in a reasonable order.

    Method pulled from QIIME (http://qiime.org), based on:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    alist = list(seq)
    alist.sort(key=_natsort_key)
    return alist

def prefer_self(x,y):
    """Merge metadata method, return X if X else Y"""
    return x if x is not None else y

def index_list(l):
    """Takes a list and returns {l[idx]:idx}"""
    return dict([(id_,idx) for idx,id_ in enumerate(l)])

def load_biom_config():
    """Returns biom-format configuration read in from file.

    Method pulled from QIIME (http://qiime.org).
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

    Method pulled from QIIME (http://qiime.org).
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

    Method pulled from QIIME (http://qiime.org).
    """
    # The biom_config object is a default dict: if keys are not present, None
    # is returned.
    def return_none():
        return None

    results = defaultdict(return_none)
    for biom_config_file in biom_config_files:
        try:
            results.update(parse_biom_config_file(biom_config_file))
        except IOError:
            pass

    return results

def parse_biom_config_file(biom_config_file):
    """Parses lines in a biom_config file.

    Method pulled from QIIME (http://qiime.org).
    """
    result = {}
    for line in biom_config_file:
        line = line.strip()

        # Ignore blank lines or lines beginning with '#'.
        if not line or line.startswith('#'): continue
        fields = line.split()
        param_id = fields[0]
        param_value = ' '.join(fields[1:]) or None
        result[param_id] = param_value
    return result
