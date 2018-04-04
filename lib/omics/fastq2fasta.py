"""
Convert fastq into fasta
"""
import argparse
import sys

from . import get_argparser


def convert(data, output, check=True):
    """
    Convert data from FASTQ into FASTA format

    :param data: File-like object with input data
    :param output: File-like object for output
    """
    state = 0
    for line in data:
        if check:
            if state is 0 and not line.startswith('@'):
                raise RuntimeError('Input not in FASTQ format? Expected line '
                                   'to start with @: {}'.format(line))
            elif state is 2 and not line == '+\n':
                raise RuntimeError('Input not in FASTQ format? Expected + '
                                   'separator line: {}'.format(line))

        if state is 0:
            output.write('>' + line[1:])
        elif state is 1:
            output.write(line)

        state = (state + 1) % 4

    if state is not 0:
        raise RuntimeError('Input not in FASTQ format? Expected total number '
                           'of lines to be multiple of 4, last line: {}'
                           ''.format(line))


def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
        project_home=False,
        threads=False,
    )
    argp.add_argument(
        'inputfile',
        metavar='FILE',
        nargs='?',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help='Fastq file to be converted, by default data is read from stdin.'
    )
    argp.add_argument(
        '-o', '--output',
        metavar='FILE',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Name of output filie.  Write to stdout by default.'
    )
    argp.add_argument(
        '--force', '-f',
        action='store_true',
        help='Overwrite existing files',
    )
    argp.add_argument(
        '--no-check',
        action='store_false',
        dest='check',
        help='Skip sanity check on input data.  By default it is checked that '
             'the input is indeed in fastq format.',
    )
    args = argp.parse_args()
    convert(args.inputfile, args.output, check=args.check)


if __name__ == '__main__':
    main()
