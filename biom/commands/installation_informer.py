#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from sys import platform, version as python_version, executable
from pyqi.core.command import Command, CommandOut, ParameterCollection

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Daniel McDonald", "Jose Clemente", "Greg Caporaso",
               "Jai Ram Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", "Sue Huse"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class InstallationInformer(Command):
    BriefDescription = "Provide information about the biom-format installation"
    LongDescription = ("Provide information about the biom-format "
                       "installation, including settings pulled from the "
                       "configuration file. For more details, see "
                       "http://biom-format.org")
    CommandIns = ParameterCollection([])
    CommandOuts = ParameterCollection([
        CommandOut(Name='install_info_lines',
                   DataType='str',
                   Description='Installation info')
    ])

    def run(self, **kwargs):
        lines = []

        lines.extend(self.getFormattedSystemInfo())
        lines.extend(self.getFormattedDependencyVersionInfo())
        lines.extend(self.getFormattedPackageInfo())
        lines.append('')

        return {'install_info_lines': lines}

    def getFormattedSystemInfo(self):
        return self._format_info(self.getSystemInfo(), 'System information')

    def getFormattedDependencyVersionInfo(self):
        return self._format_info(self.getDependencyVersionInfo(),
                                 'Dependency versions')

    def getFormattedPackageInfo(self):
        return self._format_info(self.getPackageInfo(),
                                 'biom-format package information')

    def getSystemInfo(self):
        return (("Platform", platform),
                ("Python/GCC version", python_version.replace('\n', ' ')),
                ("Python executable", executable))

    def getDependencyVersionInfo(self):
        not_installed_msg = "Not installed"

        try:
            from pyqi import __version__ as pyqi_lib_version
        except ImportError:
            pyqi_lib_version = not_installed_msg

        try:
            from numpy import __version__ as numpy_lib_version
        except ImportError:
            numpy_lib_version = ("ERROR: Not installed - this is required! "
                                 "(This will also cause the BIOM library to "
                                 "not be importable.)")

        try:
            from scipy import __version__ as scipy_lib_version
        except ImportError:
            scipy_lib_version = not_installed_msg

        try:
            from dateutil import __version__ as dateutil_lib_version
        except ImportError:
            dateutil_lib_version = not_installed_msg

        return (("pyqi version", pyqi_lib_version),
                ("NumPy version", numpy_lib_version),
                ("SciPy version", scipy_lib_version),
                ("dateutil version", dateutil_lib_version))

    def getPackageInfo(self):
        import_error_msg = ("ERROR: Can't find the BIOM library code (or "
                            "numpy) - is it installed and in your "
                            "$PYTHONPATH?")
        try:
            from biom import __version__ as biom_lib_version
        except ImportError:
            biom_lib_version = import_error_msg

        try:
            from biom.exception import InvalidSparseBackendException
            from biom.table import SparseObj
            backend_name = SparseObj(0, 0).__class__.__name__
        except ImportError:
            backend_name = import_error_msg
        except InvalidSparseBackendException as e:
            backend_name = "ERROR: %s" % e

        return (("biom-format version", biom_lib_version),
                ("SparseObj type", backend_name))

    def _format_info(self, info, title):
        max_len = self._get_max_length(info)

        lines = ['']
        lines.append(title)
        lines.append('=' * len(title))
        for e in info:
            lines.append("%*s:\t%s" % (max_len, e[0], e[1]))

        return lines

    def _get_max_length(self, info):
        return max([len(e[0]) for e in info])

CommandConstructor = InstallationInformer
