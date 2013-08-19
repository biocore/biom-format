#!/usr/bin/env python

from __future__ import division
from sys import platform, version as python_version, executable
from pyqi.core.command import Command, Parameter, ParameterCollection

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BIOM-Format project"
__credits__ = ["Daniel McDonald", "Jose Clemente", "Greg Caporaso",
               "Jai Ram Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", "Sue Huse"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

class InstallationInformer(Command):
    BriefDescription = "Provide information about the biom-format installation"
    LongDescription = ("Provide information about the biom-format "
                       "installation, including settings pulled from the "
                       "configuration file. For more details, see "
                       "http://biom-format.org")
    Parameters = ParameterCollection([])

    def run(self, **kwargs):
        lines = []

        system_info = self.getSystemInfo()
        max_len =  self._get_max_length(system_info)

        lines.append('')
        lines.append("System information")
        lines.append("==================")
        for e in system_info:
            lines.append("%*s:\t%s" % (max_len, e[0], e[1]))

        version_info = self.getDependencyVersionInfo()
        max_len =  self._get_max_length(version_info)

        lines.append('')
        lines.append("Dependency versions")
        lines.append("===================")
        for e in version_info:
            lines.append("%*s:\t%s" % (max_len, e[0], e[1]))

        package_info = self.getPackageInfo()
        max_len =  self._get_max_length(package_info)

        lines.append('')
        lines.append("biom-format package information")
        lines.append("===============================")
        for e in package_info:
            lines.append("%*s:\t%s" % (max_len, e[0], e[1]))
        lines.append('')

        return {'install_info_lines': lines}

    def getSystemInfo(self):
        return (("Platform", platform),
                ("Python/GCC version", python_version.replace('\n', ' ')),
                ("Python executable", executable))

    def getDependencyVersionInfo(self):
        try:
            from pyqi import __version__ as pyqi_lib_version
        except ImportError:
            pyqi_lib_version = "Not installed"

        try:
            from numpy import __version__ as numpy_lib_version
        except ImportError:
            numpy_lib_version = ("ERROR: Not installed - this is required! "
                                 "(This will also cause the BIOM library to "
                                 "not be importable.)")

        return (("pyqi version", pyqi_lib_version),
                ("NumPy version", numpy_lib_version))

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
        except ImportError:
            SparseObj = import_error_msg
        except InvalidSparseBackendException as e:
            SparseObj = "ERROR: %s" % e

        return (("biom-format version", biom_lib_version),
                ("SparseObj type", SparseObj))

    def _get_max_length(self, items):
        return max([len(e[0]) for e in items])

CommandConstructor = InstallationInformer
