#!/usr/bin/env python

#from setuptools import find_packages
from distutils.core import setup
from distutils.extension import Extension
from os.path import join, split
from os import getcwd
from glob import glob

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, The BIOM Format"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente"]
__license__ = "GPL"
__version__ = "0.9.3-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"


try:
    import numpy
except ImportError:
    raise ImportError, "numpy cannot be found. Can't continue."

try:
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError, "Cython cannot be found. Can't continue."


support_code_path = join(getcwd(), 'python-code', 'support-code')
library_path = split(numpy.__file__)[0]

long_description = """BIOM: Biological Observation Matrix
http://www.biom-format.org

The Biological Observation Matrix (BIOM) Format or: How I Learned To Stop Worrying and Love the Ome-ome
Daniel McDonald, Jose C Clemente, Justin Kuczynski, Jai Ram Rideout,
Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John
Hufnagle, Folker Meyer, Rob Knight and J Gregory Caporaso.
GigaScience, accepted 2012.
"""

from os import system
system('cython %s/_sparsemat.pyx -o %s/_sparsemat.cpp --cplus' % (support_code_path, support_code_path))

from subprocess import Popen, PIPE
gcc_version = map(int, Popen('gcc -dumpversion', shell=True, stdout=PIPE).stdout.read().split('.'))
if gcc_version[0] == 4 and gcc_version[1] > 2:
    extra_compile_args = ['-std=c++0x']
else:
    extra_compile_args = []

setup(name='biom-format',
      version=__version__,
      description='Biological Observation Matrix',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://www.biom-format.org',
      packages=['biom'],
      scripts=glob('scripts/*py'),
      package_dir={'':'python-code'},
      package_data={},
      data_files={},
      long_description=long_description,
      ext_modules=[Extension(
                   "biom._sparsemat",
                   sources=['python-code/support-code/_sparsemat.cpp',
                            'python-code/support-code/sparsemat_lib.cpp'],
                   language="c++",
                   extra_compile_args=extra_compile_args,
                   library_dirs=[library_path],
                   include_dirs=[])],
      cmdclass={'build_ext': build_ext}
)



#setup()
