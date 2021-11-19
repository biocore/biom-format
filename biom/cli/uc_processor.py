# ----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


import click

from biom.cli import cli
from biom.cli.util import write_biom_table
from biom.parse import parse_uc
from biom.exception import TableException


@cli.command('from-uc')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input uc filepath.')
@click.option('-o', '--output-fp', required=True,
              type=click.Path(writable=True),
              help='The output BIOM filepath')
@click.option('--rep-set-fp', type=click.Path(exists=True, dir_okay=False),
              help="Fasta file containing representative sequences with "
                   "where sequences are labeled with OTU identifiers, and "
                   "description fields contain original sequence identifiers. "
                   "This output is created, for example, by vsearch with the "
                   "--relabel_sha1 --relabel_keep options.",
              required=False)
def from_uc(input_fp, output_fp, rep_set_fp):
    """Create a BIOM table from a vsearch/uclust/usearch BIOM file.

    Example usage:

    Simple BIOM creation:

    $ biom from-uc -i in.uc -o out.biom

    BIOM creation with OTU re-naming:

    $ biom from-uc -i in.uc -o out.biom --rep-set-fp rep-set.fna

    """
    input_f = open(input_fp)
    if rep_set_fp is not None:
        rep_set_f = open(rep_set_fp)
    else:
        rep_set_f = None
    table = _from_uc(input_f, rep_set_f)
    write_biom_table(table, 'hdf5', output_fp)


def _id_map_from_fasta(fasta_lines):
    result = {}
    for line in fasta_lines:
        if line.startswith('>'):
            try:
                obs_id, seq_id = line.split()[:2]
            except ValueError:
                raise ValueError('Sequence identifiers in fasta file '
                                 'must contain at least two space-'
                                 'separated fields.')
            result[seq_id] = obs_id[1:]
        else:
            pass
    return result


def _from_uc(input_f, rep_set_f=None):
    table = parse_uc(input_f)

    if rep_set_f is not None:
        obs_id_map = _id_map_from_fasta(rep_set_f)
        try:
            table.update_ids(obs_id_map, axis='observation', strict=True,
                             inplace=True)
        except TableException:
            raise ValueError('Not all sequence identifiers in the input BIOM '
                             'file are present in description fields in the '
                             'representative sequence fasta file.')

    return table
