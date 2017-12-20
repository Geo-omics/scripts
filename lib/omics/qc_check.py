"""
Run QC on metagenomic reads from multiple samples
"""
from pathlib import Path
import sys
import zipfile

from omics import get_argparser

DEFAULT_FWD_PREFIX = 'fwd'
DEFAULT_REV_PREFIX = 'rev'
DEFAULT_POST_QC_INFIX = 'derep_scythe_sickle'


def get_fail(*paths, fastqc_dir='FASTQC', prefixes=None, infixes=None,
             warnings=False):
    """
    Check results of FASTQC

    :param paths: List of pathlib Path pointing to directory for ech sample
    :param str fastqc_dir: Name of sub-directory with FASTQC results
    :param list prefixes: List of prefixes to build analysed files' name
    :param list infixes: List of infixes to build analysed files' name
    :param bool warning: Report warnings, too.  By default only fails are
                         reported
    """
    if prefixes is None:
        # what omics-qc created by default
        prefixes = ['fwd', 'rev']

    if infixes is None:
        # what omics-qc creates by default
        infixes = ['derep_scythe_sickle']  # post-qc files only

    query = ['FAIL']
    if warnings:
        query.append('WARN')

    if not prefixes:
        # if empty list was passed
        prefixes = ['']

    # compile names of result files
    for p in paths:
        for pref in prefixes:
            for inf in infixes:
                if pref and inf:
                    basename = pref + '_' + inf
                else:
                    # one is empty (at least) so no _ separator
                    basename = pref + inf
                basename += '_fastqc'

                zipf = p / fastqc_dir / '{}.zip'.format(basename)
                with zipfile.ZipFile(str(zipf)) as zf:
                    summary_name = '{}/summary.txt'.format(basename)
                    with zf.open(summary_name) as summary:
                        for line in summary:
                            mark, test_name, _ = line.decode().split('\t')
                            if mark in query:
                                yield p, pref, inf, mark, test_name


def get_argp():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
    )
    argp.add_argument(
        'samples',
        nargs='*',
        default=['.'],
        help='List of directories, one per sample that contain the sample\'s '
             'reads. The default is to take the current directory and '
             'process a single sample.  The names of the reads files must be '
             'fwd.fastq and rev.fastq, currently this can not be set manually.'
             ' Use the omics-qc-sample script directly to specify filenames, '
             'omics-qc is just a wrapper after all.'
    )
    argp.add_argument(
        '-w', '--warnings',
        action='store_true',
        help='Also report warnings',
    )
    argp.add_argument(
        '-a', '--all',
        action='store_true',
        help='Report both, pre- and post qc results.  By default only report '
             'post-qc results',
    )
    argp.add_argument(
        '--fwd',
        default=DEFAULT_FWD_PREFIX,
        help=('Basename of forward reads file, by default this is \'{}\''
              ''.format(DEFAULT_FWD_PREFIX))
    )
    argp.add_argument(
        '--rev',
        default=DEFAULT_REV_PREFIX,
        help=('Basename of reverse reads file, by default this is \'{}\''
              ''.format(DEFAULT_REV_PREFIX))
    )
    argp.add_argument(
        '--post-qc-infix',
        default=DEFAULT_POST_QC_INFIX,
        help=('The infix part of the post-qc reads files, default is {}'
              ''.format(DEFAULT_POST_QC_INFIX))
    )

    return argp


def main():
    argp = get_argp()
    args = argp.parse_args()
    args.samples = [Path(i) for i in args.samples]
    for i in args.samples:
        if not i.is_dir():
            argp.error('Directory not found: {}'.format(i))

    infixes = [args.post_qc_infix]
    if args.all:
        infixes = [''] + infixes

    try:
        fails = get_fail(
            *args.samples,
            infixes=infixes,
            warnings=args.warnings,
        )
        for i in fails:
            print(*i)
    except Exception as e:
        if args.traceback:
            raise
        else:
            print('[qc-check] (error) {}: {}'.format(e.__class__.__name__, e),
                  file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
