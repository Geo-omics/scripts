"""
Module to compile statistics for bins
"""
import argparse
import csv
from pathlib import Path
import sys


def make_statistics(checkm_table, bin_stats_analyse, bin_dir, out):
    ...


def bin_stats_convert(infile):
    """
    Convert the bin_stats.analyse.tsv file into a proper dict
    """
    data = {}

    for line in infile:
        line = line.strip()
        try:
            bin_id, row = line.split('\t')
        except Exception:
            raise RuntimeError('Failed parsing {}: expected two elements '
                               'separated by a tab:\n{}'
                               ''.format(infile.name, line))

        row = row.lstrip('{').rstrip('}').split(', ')

        data[bin_id] = {}
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
        reader_func = csv.DictReader
        data = {}
    else:
        reader_func = csv.reader
        data = []

    for row in reader_func(file, delimiter=sep):
        if to_dict:
            for k, v in row.items():
                # NOTE: also converts row ids
                row[k] = type_conv(v)

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
