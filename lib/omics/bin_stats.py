"""
Module to compile statistics for bins
"""
import argparse
from pathlib import Path
import sys


def make_statistics(checkm_table, bin_stats_analyse, bin_dir, out):
    ...


def bin_stats_convert(infile):
    """
    Convert the bin_stats.analyse.tsv file into a real tab-delimited table
    """
    is_first_row = True
    for line in infile:
        line = line.strip()
        try:
            bin_id, data = line.split('\t')
        except Exception:
            raise RuntimeError('Failed parsing {}: expected two elements '
                               'separated by a tab:\n{}'
                               ''.format(infile.name, line))

        data = data.lstrip('{').rstrip('}').split(', ')
        if is_first_row:
            header = ['bin id']

        row = [bin_id]
        for i in data:
            col_id, value = i.split(': ')
            col_id = col_id.strip("'")

            if 'scaffold' in col_id:
                continue

            if is_first_row:
                header.append(col_id)

            try:
                value = int(value)
            except Exception:
                try:
                    value = float(value)
                except Exception:
                    raise RuntimeError(
                        'Failed to parse value {} in line:\n{}\nin file {}'
                        ''.format(value, line, infile.name)
                    )
            row.append(value)
        if is_first_row:
            is_first_row = False
            yield header

        yield row


class DirectoryAction(argparse.Action):
    def __init__(self, option_strings, dest, default=Path.cwd(),
                 metavar='PATH', **kwargs):
        super().__init__(option_strings, dest, default=default,
                         metavar=metavar, **kwargs)

    def __call__(self, _, np, vals, opts=None):
        setattr(np, self.dest, Path(vals))


def get_args():
    argp = argparse.ArgumentParser(description=__doc__)
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
    return argp.parse_args()


def main():
    args = get_args()
    print(args, file=sys.stderr)
    if args.bin_stats_analyse is not None:
        s = bin_stats_convert(args.bin_stats_analyse)
        for row in s:
            for i in range(len(row)):
                if isinstance(row[i], float):
                    row[i] = '{:.4f}'.format(row[i])
            print(*row, sep='\t')
    return
    make_statistics(
        checkm_table=args.checkm_table,
        bin_stats=args.bin_stats_analyse,
        bin_dir=args.bin_dir,
        out=args.out,
    )
if __name__ == '__main__':
    main()
