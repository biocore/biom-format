#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2011-2013, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011-2013, The BIOM Format Development Team"
__credits__ = ["Jai Ram Rideout"]
__license__ = "BSD"
__url__ = "http://biom-format.org"
__version__ = "1.2.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from biom.commands.installation_informer import InstallationInformer
from biom.unit_test import TestCase, main

class InstallationInformerTests(TestCase):
    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = InstallationInformer()

    def test_default(self):
        """Correctly returns info about the biom-format installation."""
        # Not really sure what to specifically test here, as this information
        # will change on a per-install basis. Just make sure the code is being
        # exercised and we have some output.
        obs = self.cmd()
        self.assertEqual(obs.keys(), ['install_info_lines'])
        self.assertTrue(len(obs['install_info_lines']) > 0)


if __name__ == "__main__":
    main()
