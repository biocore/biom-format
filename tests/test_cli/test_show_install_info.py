# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import unittest

from biom.cli.installation_informer import _show_install_info


class TestShowInstallInfo(unittest.TestCase):
    def test_default(self):
        # Not really sure what to specifically test here, as this information
        # will change on a per-install basis. Just make sure the code is being
        # exercised and we have some output.
        obs = _show_install_info()
        self.assertTrue(len(obs) > 0)


if __name__ == "__main__":
    unittest.main()
