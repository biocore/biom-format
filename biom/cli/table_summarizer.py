# -----------------------------------------------------------------------------
# Copyright (c) 2011-2017, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


from operator import itemgetter
import locale

import click
from numpy import std

from biom import load_table
from biom.cli import cli
from biom.util import compute_counts_per_sample_stats


@cli.command(name='summarize-table')
@click.option('-i', '--input-fp', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='The input BIOM table')
@click.option('-o', '--output-fp', default=None,
              type=click.Path(writable=True, dir_okay=False),
              help='An output file-path')
@click.option('--qualitative', default=False, is_flag=True,
              help="Present counts as number of unique observation ids per"
                   " sample, rather than counts of observations per sample.")
@click.option('--observations', default=False, is_flag=True,
              help="Summarize over observations")
def summarize_table(input_fp, output_fp, qualitative, observations):
    """Summarize sample or observation data in a BIOM table.

    Provides details on the observation counts per sample, including summary
    statistics, as well as metadata categories associated with samples and
    observations.

    Example usage:

    Write a summary of table.biom to table_summary.txt:

    $ biom summarize-table -i table.biom -o table_summary.txt

    """
    table = load_table(input_fp)
    result = _summarize_table(table, qualitative, observations)
    if output_fp:
        with open(output_fp, 'w') as fh:
            fh.write(result)
    else:
        click.echo(result)


def _summarize_table(table, qualitative=False, observations=False):
    lines = []
    locale.setlocale(locale.LC_ALL, '')

    if observations:
        table = table.transpose()

    min_counts, max_counts, median_counts, mean_counts, counts_per_samp =\
        compute_counts_per_sample_stats(table, qualitative)
    num_observations = len(table.ids(axis='observation'))

    counts_per_sample_values = list(counts_per_samp.values())

    if table.metadata() is None:
        sample_md_keys = ["None provided"]
    else:
        sample_md_keys = table.metadata()[0].keys()

    if table.metadata(axis='observation') is None:
        observation_md_keys = ["None provided"]
    else:
        observation_md_keys = table.metadata(axis='observation')[0].keys()

    num_samples = len(table.ids())

    if observations:
        # as this is a transpose of the original table...
        lines.append('Num samples: ' + locale.format_string('%d',
                     num_observations, grouping=True))
        lines.append('Num observations: ' + locale.format_string('%d',
                     num_samples, grouping=True))
    else:
        lines.append('Num samples: ' + locale.format_string('%d',
                     num_samples, grouping=True))
        lines.append('Num observations: ' + locale.format_string('%d',
                     num_observations, grouping=True))

    if not qualitative:
        total_count = sum(counts_per_sample_values)
        lines.append('Total count: ' + locale.format_string('%d',
                     total_count, grouping=True))
        lines.append('Table density (fraction of non-zero values): %1.3f' %
                     table.get_table_density())

    lines.append('')

    if qualitative:
        if observations:
            lines.append('Sample/observations summary:')
        else:
            lines.append('Observations/sample summary:')
    else:
        lines.append('Counts/sample summary:')

    lines.append(' Min: ' + locale.format_string('%1.3f',
                 min_counts, grouping=True))
    lines.append(' Max: ' + locale.format_string('%1.3f',
                 max_counts, grouping=True))
    lines.append(' Median: ' + locale.format_string('%1.3f',
                 median_counts, grouping=True))
    lines.append(' Mean: ' + locale.format_string('%1.3f',
                 mean_counts, grouping=True))
    lines.append(' Std. dev.: ' + locale.format_string('%1.3f',
                 std(counts_per_sample_values), grouping=True))

    if observations:
        # since this is a transpose...
        lines.append(
            ' Sample Metadata Categories: %s' %
            '; '.join(observation_md_keys))
        lines.append(
            ' Observation Metadata Categories: %s' %
            '; '.join(sample_md_keys))
        lines.append('')
    else:
        lines.append(
            ' Sample Metadata Categories: %s' %
            '; '.join(sample_md_keys))
        lines.append(
            ' Observation Metadata Categories: %s' %
            '; '.join(observation_md_keys))
        lines.append('')

    if qualitative:
        lines.append('Observations/sample detail:')
    else:
        lines.append('Counts/sample detail:')

    for k, v in sorted(counts_per_samp.items(), key=itemgetter(1)):
        lines.append('%s: ' % k + locale.format_string('%1.3f',
                     v, grouping=True))

    return "\n".join(lines)
