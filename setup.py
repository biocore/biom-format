#!/usr/bin/env python

from setuptools import find_packages
from distutils.core import setup
from glob import glob
from os import chdir, getcwd
from os.path import join, abspath
from subprocess import call

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The BIOM Format"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente"]
__license__ = "GPL"
__version__ = "0.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

long_description = """BIOM: Biological Observation Matrix
http://www.biom-format.org

The Biological Observation Matrix (BIOM) Format or: How I Learned To Stop Worrying and Love the Ome-ome
Daniel McDonald, Jose C Clemente, Justin Kuczynski, Jai Ram Rideout,
Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John
Hufnagle, Folker Meyer, Rob Knight and J Gregory Caporaso.
GigaScience, accepted 2012.
"""

doc_imports_failed = False
try:
    import sphinx
except ImportError:
    doc_imports_failed = True

def build_html():
    """ Build the sphinx documentation 
    
    The code for building sphinx documentation is based on 
    PyCogent's setup.py.
    
    """
    cwd = getcwd()
    doc_dir = join(cwd,'doc')
    chdir(doc_dir)
    call(["make", "html"])
    chdir(cwd)
    index_html_path = join(abspath(doc_dir),'_build','html','index.html')
    print "Local documentation built with Sphinx. "+\
          "Open to following path with a web browser:\n%s" %\
            index_html_path

setup(name='biom-format',
      version=__version__,
      description='Biological Observation Matrix',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://www.biom-format.org',
      packages=find_packages("python-code"),
      scripts=glob('scripts/*py'),
      package_dir={'':'python-code'},
      package_data={},
      data_files={},
      long_description=long_description,
)

if doc_imports_failed:
    print "Sphinx not installed, so cannot build local html documentation."
else:
    build_html()

