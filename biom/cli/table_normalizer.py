#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division


def normalize_table(table, relative_abund=False, presence_absence=False,
                    axis='sample'):
    if relative_abund is False and presence_absence is False:
        raise ValueError("Must specifiy a normalization type")
    elif relative_abund is True and presence_absence is True:
        raise ValueError("Must specify only one normalization type")

    if relative_abund:
        table.norm(axis=axis)
    else:
        table.pa()

    return table
