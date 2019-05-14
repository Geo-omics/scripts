# Copyright 2019 Regents of The University of Michigan.

# This file is part of geo-omics-scripts.

# Geo-omics-scripts is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# Geo-omics-scripts is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with Geo-omics-scripts.  If not, see <https://www.gnu.org/licenses/>.

"""
Module to compile statistics for bins
"""

import argparse
from collections import OrderedDict
import csv
from pathlib import Path
import sys

from . import OmicsArgParser


SUMMARY = 'summary'
TEST = 'test'

SUBCOMMANDS = [SUMMARY, TEST]


def make_statistics(args, bin_stats):
    ...


def bin_stats_convert(infile):
    """
    Convert the bin_stats.analyse.tsv file into a proper dict
    """
    data = OrderedDict()

    for line in infile:
        line = line.strip()
        try:
            bin_id, row = line.split('\t')
        except Exception:
            raise RuntimeError('Failed parsing {}: expected two elements '
                               'separated by a tab:\n{}'
                               ''.format(infile.name, line))

        row = row.lstrip('{').rstrip('}').split(', ')

        data[bin_id] = OrderedDict()
        for i in row:
            col_id, value = i.split(': ')
            col_id = col_id.strip("'")

            if 'scaffold' in col_id:
                # FIXME: for now skip scaffold info
                continue

            data[bin_id][col_id] = type_conv(value)

    return data


class DirectoryAction(argparse.Action):
    def __init__(self, option_strings, dest, default=Path.cwd(),
                 metavar='PATH', **kwargs):
        super().__init__(option_strings, dest, default=default,
                         metavar=metavar, **kwargs)

    def __call__(self, _, np, vals, opts=None):
        setattr(np, self.dest, Path(vals))


def type_conv(val):
    """ Helper function to convert numerical values """
    for t in int, float:
        try:
            val = t(val)
        except ValueError:
            pass
        else:
            break
    return val


def load_table(file, sep='\t', header=True, to_dict=True, row_id_col=0):
    """ Load a csv/tsv into dict or array"""
    if header and to_dict:
        id_col_name = file.buffer.peek().decode().splitlines()[0].split(sep)[0]
    else:
        id_col_name = None

    if to_dict:
        csv_reader_func = csv.DictReader
        data = OrderedDict()
    else:
        csv_reader_func = csv.reader
        data = []

    csv_reader = csv_reader_func(file, delimiter=sep)
    for row0 in csv_reader:
        if to_dict:
            # preserve order and convert type
            row = OrderedDict()
            for k in csv_reader.fieldnames:
                # NOTE: also converts row ids
                row[k] = type_conv(row0[k])

            try:
                rowid = row[id_col_name]
            except KeyError:
                raise RuntimeError(
                    'Expected row identifier "{}" but that is not in csv '
                    'parsed data: {}'.format(id_col_name, row)
                )

            data[rowid] = row
            del row[id_col_name]
        else:
            # FIXME: (or not) converts header row too, unlike if to_dict
            data.append([type_conv(i) for i in row])

    return data


def pick_columns(data, *cols):
    """
    Picks given columns from table data

    Returns a new dict of dicts
    """
    ret = OrderedDict()
    for row_id, row_dat in data.items():
        ret[row_id] = OrderedDict()
        for col_id in cols:
            ret[row_id][col_id] = row_dat[col_id]

    return ret


def join_tables(*tables):
    """ Join two or more tables in dict format """
    data = OrderedDict()
    for i in tables:
        for k, v in i.items():
            if k not in data:
                data[k] = OrderedDict()
            data[k].update(v)

    return data


def sort_table(table, key, rev=False):
    return OrderedDict(sorted(
        table.items(),
        key=lambda x: x[1][key],
        reverse=rev
    ))


def make_summary(stats_a, stats_checkm, args):
    """ Print summary table, selected columns from input tables """
    cols1 = ['# contigs', 'Genome size', '# predicted genes']
    cols2 = ['Marker lineage', 'Completeness', 'Contamination',
             'Strain heterogeneity']
    stats_a = pick_columns(stats_a, *cols1)
    stats_checkm = pick_columns(stats_checkm, *cols2)
    data = join_tables(stats_a, stats_checkm)
    data = sort_table(data, 'Genome size', rev=True)
    pprint_table_dict(data)


def get_args(argv=None, namespace=None):
    prog = __loader__.name.replace('.', ' ').replace('_', '-')
    argp = OmicsArgParser(prog=prog, description=__doc__, threads=False)
    argp.add_argument(
        'command',
        choices=SUBCOMMANDS,
        help='The action that will be performed.',
    )
    argp.add_argument(
        '-c', '--checkm-table',
        metavar='FILE',
        type=argparse.FileType(),
        help='CheckM output file as obtained with checkm\'s -f and '
             '--tab_table options',
    )
    argp.add_argument(
        '-a', '--bin-stats-analyse',
        metavar='FILE',
        type=argparse.FileType(),
        help='Path to CheckM\'s bin_stats_.analyse.tsv file',
    )
    argp.add_argument(
        '-b', '--bin-dir',
        action=DirectoryAction,
        help='Directory with bins.  There should be one fasta file per bin.'
             'The current working directory by default.',
    )
    argp.add_argument(
        '-o', '--out',
        metavar='FILE',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Name of output file.  By default output is printed to stdout.',
    )
    return argp.parse_args(args=argv, namespace=namespace)


def pprint_table_dict(data):
    """ Pretty print a dict data structure """

    def _format(values):
        """ float formatting helper function """
        ret = []
        for i in values:
            if isinstance(i, float):
                i = '{:.4f}'.format(i)
            ret.append(i)
        return ret

    # print header
    if data:
        print('', *list(list(data.values())[0].keys()), sep='\t')
        for k, v in data.items():
            print(k, *_format(v.values()), sep='\t')


def main(argv=None, namespace=None):
    args = get_args(argv=argv, namespace=namespace)
    if args.command == TEST:
        print(args, file=sys.stderr)

    if args.bin_stats_analyse is None:
        bin_stats_tab = []
    else:
        bin_stats_tab = bin_stats_convert(args.bin_stats_analyse)

    checkm_tab = load_table(args.checkm_table)

    make_statistics(
        args,
        bin_stats=bin_stats_tab,
    )

    if args.command == SUMMARY:
        make_summary(bin_stats_tab, checkm_tab, args)
    elif args.command == TEST:
        print('Printing bin stats table:')
        pprint_table_dict(bin_stats_tab)
        print()
        print()
        pprint_table_dict(checkm_tab)
    else:
        pass


if __name__ == '__main__':
    main()
