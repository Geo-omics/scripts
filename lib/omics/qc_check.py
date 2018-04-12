"""
Run QC on metagenomic reads from multiple samples
"""
from collections import namedtuple
from itertools import groupby
from pathlib import Path
import sys
import zipfile

from omics import get_argparser, DEFAULT_VERBOSITY

DEFAULT_FASTQC_DIR = 'FASTQC'
DEFAULT_FWD_PREFIX = 'fwd'
DEFAULT_REV_PREFIX = 'rev'
DEFAULT_POST_QC_INFIX = 'good'
BASENAME_SEPARATOR = '.'  # As used by omics-qc; Note: FASTQC also uses _ sep


def get_test_marks(*paths, fastqc_dir=DEFAULT_FASTQC_DIR, prefixes=None,
                   infixes=None, warnings=False, passes=False):
    """
    Check results of FASTQC

    :param paths: List of pathlib Path pointing to directory for ech sample
    :param str fastqc_dir: Name of sub-directory with FASTQC results
    :param list prefixes: List of prefixes to build analysed files' name
    :param list infixes: List of infixes to build analysed files' name
    :param bool warnings: Report warnings, too.  By default only fails are
                         reported
    :param bool passes: Report passes, too.  By default only fails are
                         reported
    """
    if prefixes is None:
        # what omics-qc created by default
        prefixes = ['fwd', 'rev']

    if infixes is None:
        # what omics-qc creates by default
        infixes = [DEFAULT_POST_QC_INFIX]  # post-qc files only

    query = ['FAIL']
    if warnings:
        query.append('WARN')

    if passes:
        query.append('PASS')

    if not prefixes:
        # if empty list was passed
        prefixes = ['']

    result = namedtuple('result', ['path', 'pref', 'inf', 'test', 'mark'])
    # compile names of result files
    for p in paths:
        for inf in infixes:
            for pref in prefixes:
                basename = get_basename(pref, inf)
                zipf = p / fastqc_dir / '{}.zip'.format(basename)
                with zipfile.ZipFile(str(zipf)) as zf:
                    summary_name = '{}/summary.txt'.format(basename)
                    with zf.open(summary_name) as summary:
                        for line in summary:
                            mark, test_name, _ = line.decode().split('\t')
                            if mark in query:
                                yield result(p, pref, inf, test_name, mark)


def get_basename(prefix, infix, sep=BASENAME_SEPARATOR):
    """
    Build base of relevant filenames
    """
    parts = [prefix, infix]
    parts = filter(None, parts)
    basename = sep.join(parts)
    basename += '_fastqc'
    return basename


def get_details(path, fastqc_dir, prefix, infix, test):
    """
    Retrieve a single detailed test result
    """
    basename = get_basename(prefix, infix)
    zipf = path / fastqc_dir / '{}.zip'.format(basename)
    with zipfile.ZipFile(str(zipf)) as zf:
        data_file_name = '{}/fastqc_data.txt'.format(basename)
        with zf.open(data_file_name) as data:
            ret = ''
            in_section = False
            for line in data:
                line = line.decode()
                if line.startswith('>>{}\t'.format(test)):
                    in_section = True
                    continue

                if in_section and line == '>>END_MODULE\n':
                    break

                if in_section:
                    ret += line

            return ret


def get_argp():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
        project_home=False,
        threads=False,
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
        help='Also report warnings. By default only results involving a FAIL'
             'get reported',
    )
    excl_grp = argp.add_mutually_exclusive_group()
    excl_grp.add_argument(
        '-d', '--diff',
        action='store_true',
        help='Show pre- vs. post-qc difference.  ',
    )
    excl_grp.add_argument(
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
    argp.add_argument(
        '-t', '--test',
        action='append',
        help='Restrict output to given tests, this option can be given '
             'multiple times.  It is okay to only give the beginning of'
             'the test\'s name.'
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
    if args.all or args.diff:
        # pre-qc files have empty infix
        # order such that before-qc comes before after-qc
        infixes = [''] + infixes

    try:
        results = get_test_marks(
            *args.samples,
            infixes=infixes,
            warnings=args.warnings or args.diff,
            passes=args.diff,
        )

        if args.test:
            # show all tests unless restircted by command line
            results = [
                i for i in results
                if any([i.test.startswith(j) for j in args.test])
            ]

        if args.diff:
            def diffsort(x):
                return (x.path, x.pref, x.test)

            results = list(sorted(results, key=diffsort))

            for (path, pref, test), results in groupby(results, diffsort):
                try:
                    before, after = list(results)
                except ValueError:
                    raise RuntimeError(
                        'Unexpected number of results for {}/{}/{}: {}'
                        ''.format(path, pref, test, results)
                    )
                else:
                    if before.mark != after.mark:
                        if args.warnings \
                                or 'FAIL' in [before.mark, after.mark]:
                            print(
                                '[{} {}]'.format(path, pref),
                                '{} --> {} ({})'
                                ''.format(before.mark, after.mark, test),
                            )
        else:
            for i in results:
                if args.all:
                    if i.inf == infixes[0]:
                        order = ' before'
                    elif i.inf == infixes[1]:
                        order = '  after'
                    else:
                        raise ValueError('bad infix value: {}'.format(i.inf))
                else:
                    order = ''

                if args.verbosity > DEFAULT_VERBOSITY:
                    print('[{} {}{}] {} {} -- Detailed report:'
                          ''.format(i.path, i.pref, order, i.mark, i.test))
                    print(get_details(i.path, DEFAULT_FASTQC_DIR, i.pref,
                                      i.inf, i.test))
                else:
                    print(
                        '[{} {}{}]'.format(i.path, i.pref, order),
                        '{} ({})'.format(i.mark, i.test),
                    )

    except Exception as e:
        if args.traceback:
            raise
        else:
            print('[qc-check] (error) {}: {}'.format(e.__class__.__name__, e),
                  file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
