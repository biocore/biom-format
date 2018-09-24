#!/bin/bash
set -xe

table=examples/min_sparse_otu_table_hdf5.biom
obsmd=examples/obs_md.txt
if [[ ! -f ${table} ]];
then
    echo "This script expects to operate in the base repository directory"
    exit 1
fi

biom show-install-info
biom validate-table -i $table
biom add-metadata -i $table -o test_add_metadata.biom --observation-metadata-fp $obsmd --sc-separated taxonomy
biom convert -i $table -o test_json.biom --to-json
biom convert -i $table -o test_tsv.txt --to-tsv
biom validate-table -i test_json.biom
biom convert -i test_json.biom -o test_hdf5.biom --to-hdf5
biom validate-table -i test_hdf5.biom
biom convert -i test_json.biom -o test_tsv.txt --to-tsv
biom convert -i test_tsv.txt -o test_hdf5.biom --to-hdf5 --table-type "OTU table"
biom validate-table -i test_hdf5.biom
biom convert -i test_tsv.txt -o test_json.biom --to-json --table-type "OTU table"
biom validate-table -i test_json.biom
biom head -i $table
biom table-ids -i $table
biom normalize-table -i $table -o test_norm.biom -p
cat << EOF >> ids.txt
GG_OTU_3
GG_OTU_5
EOF
biom subset-table -i $table -o test_subset.biom -a observation -s ids.txt
biom summarize-table -i $table

cat << EOF >> uctest
# uclust --input /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T/UclustExactMatchFilterrW47Ju.fasta --id 0.97 --tmpdir /var/folders/xq/0kh93ng53bs6zzk091w_bbsr0000gn/T --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc dn-otus/uclust_picked_otus/seqs_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S   0   133 *   *   *   *   *   f2_1539 *
S   0   133 *   *   *   *   *   f3_1540 *
H   0   141 100.0   +   0   0   133M8D  f3_42   f2_1539
EOF

biom from-uc -i uctest -o test_uc.biom

echo "If this message was shown, then all commands executed without failure"
