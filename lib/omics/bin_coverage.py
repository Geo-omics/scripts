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

"Calculate per-bin mean coverage"

import argparse
from itertools import groupby
from operator import itemgetter
from pathlib import Path
import sys

from . import OmicsArgParser


# suffix of fasta files in bin directory
BIN_FASTA_SUFFIX = 'fa'


def main(argv=None, namespace=None):
    args, cov_files = get_args(argv, namespace)
    bins = get_bins(args.bins)
    for row in compile_mean_coverage(bins, *cov_files):
        print(*row, sep='\t', file=args.outfile)


def load_cov_data(bins, file):
    """
    Load coverage file (genomeCoverageBed output)

    Will filter and sort by contig
    """
    data = []
    with file.open() as f:
        row0 = f.readline()
        if not row0.startswith('Contig\tDepth\t'):
            # omics-mapping actually sets this header
            raise RuntimeError(
                'Unexpected first line in {}, should be genomeCovBed '
                'header but got: {}'.format(f.name, row0)
            )

        for line in f:
            row = line.strip().split('\t')
            if len(row) != 5:
                raise RuntimeError('Expected 5 tab-separated columns in '
                                   '{} but got: {}'.format(f.name, line))
            if '.' in row[0]:
                contig, _, chunk = row[0].rpartition('.')
            else:
                contig = row[0]
                chunk = None

            if contig not in bins:
                continue

            if chunk is not None:
                try:
                    int(chunk)
                except ValueError:
                    raise RuntimeError(
                        'Failed parsing contig id, found chunk but failed '
                        'to parse it as int: {}'.format(line)
                    )
            try:
                depth = int(row[1])
            except ValueError:
                raise RuntimeError('Failed parsing depth as int in {}, '
                                   'line: {}'.format(f.name, line))

            try:
                bases = int(row[2])
            except ValueError:
                raise RuntimeError('Failed parsing bases_with_depth as int'
                                   'in {}, line: {}'.format(f.name, line))

            try:
                length = int(row[3])
            except ValueError:
                raise RuntimeError('Failed parsing contig length as int'
                                   'in {}, line: {}'.format(f.name, line))

            data.append((contig, chunk, depth, bases, length))

    return sorted(data, key=itemgetter(0))


def compile_mean_coverage(bins, *cov_files, header=True):
    """ compile per bin and sample mean coverage """

    def bin_id_to_num(bin_id):
        """ key function for bin id sorting """
        digits = [0]
        for i in bin_id:
            try:
                digits.append(int(i))
            except ValueError:
                pass
        return int(''.join(map(str, digits)))

    bin_list = sorted(set(bins.values()), key=bin_id_to_num)

    print('Got {} bins with total {} contigs'
          ''.format(len(bin_list), len(bins.keys())), file=sys.stderr)

    header = ['sample'] + bin_list
    if header:
        yield header

    for sample, cov_file in cov_files:

        # initialize zero bases and length for each bin
        bin_data = {i: [0, 0] for i in bin_list}

        # get filtered and sorted per-chunk and depth coverage data
        cov_data = load_cov_data(bins, cov_file)

        for contig, contig_rows in groupby(cov_data, key=itemgetter(0)):
            total_bases = 0
            length = 0
            chunk_count = 0
            for chunk, chunk_rows in groupby(contig_rows, key=itemgetter(1)):
                chunk_count += 1
                first = True
                for row in chunk_rows:
                    total_bases += row[2] * row[3]  # depth * bases
                    if first:
                        # add length to contig (only once per chunk)
                        length += row[4]
                        first = False

                # if chunk_count > 1:
                #    print(contig, 'multiple chunks!', file=sys.stderr)

            bin_data[bins[contig]][0] += total_bases
            bin_data[bins[contig]][1] += length
            # print(sample, contig, total_bases, length, total_bases / length,
            #      sep='\t', file=sys.stderr)

        # yield row with coverages in consistent order
        yield [sample] + [bin_data[i][0] / bin_data[i][1] for i in bin_list]


def get_bins(path):
    """ Compile contig to bin mapping """
    bins = {}
    for i in Path(path).glob('*.{}'.format(BIN_FASTA_SUFFIX)):
        with i.open() as f:
            for line in f:
                if line.startswith('>'):
                    contig = line[1:].split()[0]
                    bins[contig] = i.stem

    if not bins:
        raise RuntimeError('Failed to find bin fasta files')

    return bins


def get_args(argv=None, namespace=None):
    prog = __loader__.name.replace('.', ' ').replace('_', '-')
    argp = OmicsArgParser(prog=prog, description=__doc__, threads=False)
    argp.add_argument(
        'sample',
        nargs='+',
        help='Sample directory names',
    )
    argp.add_argument(
        '--coverage-path',
        help='Path to the sample\'s coverage file relative to the sample '
             'directory',
    )
    argp.add_argument(
        '--bins',
        help='Path to the bins',
    )
    argp.add_argument(
        '-o', '--outfile',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='File to write table with per sample per bin coverage data',
    )
    args = argp.parse_args(args=argv, namespace=namespace)

    # make pairs of sample id matching coverage file
    cov_files = [
        (Path(i).name, Path(i) / args.coverage_path)
        for i in args.sample
    ]
    for _, i in cov_files:
        if not i.is_file():
            argp.error('File not found: {}'.format(i))

    return args, cov_files


if __name__ == '__main__':
    main()
