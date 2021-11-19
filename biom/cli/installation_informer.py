# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


import sys

import click

from biom.cli import cli


@cli.command(name='show-install-info')
def show_install_info():
    """Provide information about the biom-format installation.

    Provide information about the biom-format installation, including settings
    pulled from the configuration file. For more details, see
    http://biom-format.org

    Example usage:

    Display biom-format installation information:

    $ biom show-install-info

    """
    click.echo(_show_install_info())


def _show_install_info():
    lines = []
    lines.extend(_get_formatted_system_info())
    lines.extend(_get_formatted_dependency_version_info())
    lines.extend(_get_formatted_package_info())
    lines.append('')
    return '\n'.join(lines)


def _get_formatted_system_info():
    return _format_info(_get_system_info(), 'System information')


def _get_formatted_dependency_version_info():
    return _format_info(_get_dependency_version_info(), 'Dependency versions')


def _get_formatted_package_info():
    return _format_info(_get_package_info(), 'biom-format package information')


def _get_system_info():
    return (("Platform", sys.platform),
            ("Python version", sys.version.replace('\n', ' ')),
            ("Python executable", sys.executable))


def _get_dependency_version_info():
    not_installed_msg = "Not installed"

    try:
        from click import __version__ as click_lib_version
    except ImportError:
        click_lib_version = not_installed_msg

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
        from h5py import __version__ as h5py_lib_version
    except ImportError:
        h5py_lib_version = ("WARNING: Not installed - this is an optional "
                            "dependency. It is strongly recommended for "
                            "large datasets.")

    return (("click version", click_lib_version),
            ("NumPy version", numpy_lib_version),
            ("SciPy version", scipy_lib_version),
            ("h5py version", h5py_lib_version))


def _get_package_info():
    import_error_msg = ("ERROR: Can't find the BIOM library code (or "
                        "numpy) - is it installed and in your "
                        "$PYTHONPATH?")
    try:
        from biom import __version__ as biom_lib_version
    except ImportError:
        biom_lib_version = import_error_msg

    return (("biom-format version", biom_lib_version),)


def _format_info(info, title):
    max_len = _get_max_length(info)

    lines = ['']
    lines.append(title)
    lines.append('=' * len(title))
    for e in info:
        lines.append("%*s:\t%s" % (max_len, e[0], e[1]))

    return lines


def _get_max_length(info):
    return max(len(e[0]) for e in info)
