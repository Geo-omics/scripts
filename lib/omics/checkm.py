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
Utilities to work with CheckM output
"""

import argparse
from collections import OrderedDict
import sys

from . import OmicsArgParser

ACTIONS = ['convert']
DEFAULT_ACTION = ACTIONS[0]


def load_tsv(file):
    """
    Load data from CheckM-style tsv output file

    :param file: Filehandle
    :return: List of OrderedDicts, one per rows
    """
    cols = None

    data = []
    for line in file:
        line = line.split('\t')
        dat = line[1].strip().lstrip('{').rstrip('}')
        dat = dat.split(', ')
        dat = map(lambda x: x.split(': '), dat)

        dd = OrderedDict()
        dd['bin'] = line[0]

        for k, v in dat:
            dd[k.strip("'")] = v

        if cols is None:
            cols = dd.keys()
        else:
            if cols != dd.keys():
                raise RuntimeError('data keys, inconsistency: {} vs. {}'
                                   ''.format(cols, dd.keys()))

        for k, v in dd.items():
            try:
                v = int(v)
            except ValueError:
                try:
                    v = float(v)
                except ValueError:
                    pass
            dd[k] = v

        data.append(dd)
    return data


def main(argv=None, namespace=None):
    prog = __loader__.name.replace('.', ' ')
    argp = OmicsArgParser(prog=prog, description=__doc__, threads=False)
    argp.add_argument(
        '-c', '--convert',
        action='store_true',
        help='Convert input to real tab-separated table',
    )
    argp.add_argument(
        'inputfile',
        metavar='FILE',
        type=argparse.FileType(),
        default=sys.stdin,
        help='Input file, usually a .tsv file written by CheckM',
    )
    args = argp.parse_args(args=argv, namespace=namespace)

    if args.convert:
        data = load_tsv(args.inputfile)
        print(*data[0].keys(), sep='\t')
        for row in data:
            print(*row.values(), sep='\t')
    else:
        argp.error('no action specified, e.g. --convert')


if __name__ == '__main__':
    main()
