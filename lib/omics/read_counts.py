import argparse
from itertools import zip_longest
from operator import itemgetter
from pathlib import Path
import sys

import matplotlib

from . import get_argparser


def main(argv=None):
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__
    )
    argp.add_argument(
        'source',
        nargs='*',
        default=['.'],
        help='List of per-sample directories of fasta/q files for input.',
    )
    argp.add_argument(
        '-f', '--file-name',
        metavar='FILE',
        default='fwd.fastq',
        help='Filename of fasta/q input files that any directories in '
             'source will be searched for.',
    )
    argp.add_argument(
        '-o', '--out-txt',
        metavar='PATH',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='File to write text output to.  Stdout by default',
    )
    argp.add_argument(
        '--pdf',
        metavar='PATH',
        type=argparse.FileType('w'),
        help='Save a read-count plot to a PDF file of given name',
    )
    args = argp.parse_args(argv)

    files = []  # list of tuples (sample id, path)
    for i in map(Path, args.source):
        if i.is_dir():
            file = i / args.file_name
            if file.is_file():
                files.append((i.name, file))
            else:
                print('Warning: File not found in directory {}: {}'
                      ''.format(i, file), file=sys.stderr)
        elif i.is_file():
            files.append((i.name, i))
        else:
            argp.error('File or directory not found: {}'.format(i))

    files = sorted(set(files))

    if not files:
        argp.error('No files found')

    files = list(set(files))
    read_counts = {}

    for sample, path in files:
        if sample in read_counts:
            print('[WARNING] Ambigious sample ids: {} -- keeping only one '
                  'read count', file=sys.stderr)
        else:
            read_counts[sample] = count_fastaq_reads(path)


def count_fastaq_reads(path):
    """
    Count number of reads in fasta/q file
    """
    file = path.open('rb')
    marker = file.peek()[:1]
    if marker == b'@':
        # FASTQ
        read_lines = 4
    elif marker == b'>':
        # FASTA
        read_lines = 2
    else:
        raise RuntimeError(
            'Failed to detect fileformat: {} is neither valid FASTA nor FASTQ:'
            ' file starts: "{}"'.format(path, file.peek()[:10].decode())
        )

    count = 0
    args = [iter(file)] * read_lines
    for read in zip_longest(*args):
        if read[-1] is None:
            raise RuntimeError(
                'Line count is not a multiple of {}: {}, read count at {}\n'
                'Last lines are:\n{}'
                ''.format(path, count, '\n'.join([str(i) for i in read]))
            )
        count += 1
    return count


def make_output(read_counts, outfile, pdfout=None):
    """
    Write counts to file(s)

    :param dict read_counts: Dict mapping sample id to read count
    :param file outfile: file-like object
    :param pdfout: String of PDF file or file-like object to save
                   PDF figure or None to disable PDF output
    """
    for sample, count in sorted(read_counts.items()):
        out = '{}\t{}\n'.format(sample, count)
        outfile.write(out)

    if pdfout is not None:
        plot(read_counts, pdfout)


def plot(read_counts, outfile):
    """
    Plot counts

    :param dict read_counts: Dict mapping sample id to read count
    :param outfile: String containing name of output file or a file-like
                    object
    """
    matplotlib.use('pdf')
    from matplotlib import pyplot

    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    data = sorted(read_counts.items(), key=itemgetter(1), reverse=True)
    samples = [i[0] for i in data]
    counts = [i[1] for i in data]
    idx = range(len(data))
    ax.bar(idx, counts)
    pyplot.title('Paired-read count per sample')
    pyplot.ylabel('read count')
    pyplot.xticks(idx, samples)
    fig.savefig('read_counts.pdf')
    pyplot.close()


if __name__ == '__main__':
    main()
