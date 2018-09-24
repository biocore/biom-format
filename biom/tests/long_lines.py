# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


# no hits or library seeds
uc_empty = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
"""

# label not in qiime post-split-libraries format
uc_invalid_id = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	133	*	*	*	*	*	1539	*
"""

# contains single new (de novo) seed hit
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

# contains single library (reference) seed hit
uc_lib_minimal = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	f2_1539	295053
"""

# contains new seed (de novo) hits only
uc_seed_hits = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	133	*	*	*	*	*	f2_1539	*
H	0	141	100.0	+	0	0	133M8D	f3_42	f2_1539
H	0	141	100.0	+	0	0	133M8D	f2_43	f2_1539
S	0	133	*	*	*	*	*	f3_44	*
"""

# contains library (reference) and new seed (de novo) hits
uc_mixed_hits = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	133	*	*	*	*	*	f2_1539	*
H	0	141	100.0	+	0	0	133M8D	f3_42	f2_1539
H	0	141	100.0	+	0	0	133M8D	f2_43	f2_1539
S	0	133	*	*	*	*	*	f3_44	*
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	f2_1539	295053
"""

# contains library (reference) and new seed (de novo) hits
# and sample ids contain underscores
uc_underscores_in_sample_id = """# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	133	*	*	*	*	*	_f_2__1539	*
H	0	141	100.0	+	0	0	133M8D	f_3_42	_f_2__1539
H	0	141	100.0	+	0	0	133M8D	_f_2__43	_f_2__1539
S	0	133	*	*	*	*	*	f_3_44	*
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	_f_2__1539	295053
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
