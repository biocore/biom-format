#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import tempfile
from unittest import TestCase, main

import numpy as np

import biom
from biom.cli.uc_processor import _from_uc

class TestUcProcessor(TestCase):

    def setUp(self):
        """Set up data for use in unit tests."""
        self.cmd = _from_uc
        self.uc_minimal = uc_minimal.split('\n')
        self.uc = uc.split('\n')
        self.rep_set = rep_set.split('\n')
        self.rep_set_no_mapping = rep_set_no_mapping.split('\n')
        self.rep_set_missing_id = rep_set_missing_id.split('\n')

    def test_basic(self):
        obs = self.cmd(self.uc_minimal)
        expected = biom.Table(np.array([[1.0]]),
                              observation_ids=['f2_1539'],
                              sample_ids=['f2'])
        self.assertEqual(obs, expected)

    def test_basic_w_mapping(self):
        obs = self.cmd(self.uc_minimal, self.rep_set)
        expected = biom.Table(np.array([[1.0]]),
                              observation_ids=['otu1'],
                              sample_ids=['f2'])
        self.assertEqual(obs, expected)

    def test_rep_set_no_mapping(self):
        self.assertRaises(ValueError, self.cmd, self.uc_minimal,
                          self.rep_set_no_mapping)

    def test_rep_set_missing_id(self):
        self.assertRaises(ValueError, self.cmd, self.uc_minimal,
                          self.rep_set_missing_id)

    def test_uc(self):
        obs = self.cmd(self.uc)
        expected = biom.Table(np.array([[1.0, 1.0], [0.0, 1.0]]),
                              observation_ids=['f2_1539', 'f3_1540'],
                              sample_ids=['f2', 'f3'])
        self.assertEqual(obs, expected)

    def test_uc_w_mapping(self):
        obs = self.cmd(self.uc, self.rep_set)
        expected = biom.Table(np.array([[1.0, 1.0], [0.0, 1.0]]),
                              observation_ids=['otu1', 'otu2'],
                              sample_ids=['f2', 'f3'])
        self.assertEqual(obs, expected)

uc_minimal = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	133	*	*	*	*	*	f2_1539	*
"""

uc = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	133	*	*	*	*	*	f2_1539	*
S	0	133	*	*	*	*	*	f3_1540	*
H	0	141	100.0	+	0	0	133M8D	f3_42	f2_1539
"""

rep_set = """>otu1 f2_1539
ACGT
>otu2 f3_1540
ACCT
"""

rep_set_no_mapping = """>otu1
ACGT
>otu2
ACCT
"""

rep_set_missing_id = """>otu1 f99_1539
ACGT
>otu2 f99_1539
ACCT
"""

if __name__ == '__main__':
    main()
