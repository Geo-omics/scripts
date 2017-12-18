"""
Run QC on metagenomic reads from multiple samples
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys
import zipfile

from omics import get_argparser, DEFAULT_THREADS, DEFAULT_VERBOSITY


def qc_sample(path, *, clean_only=False, adapters=None, keep_all=False,
              less_mem=False, no_dereplicate=False, no_interleave=False,
              no_fasta_interleave=False, verbosity=DEFAULT_VERBOSITY,
              project=None):
    """
    Do QC for a single sample

    This is a wrapper for the omics-qc-sample script
    """
    script = 'omics-qc-sample'

    args = []
    not clean_only or args.append('--clean-only')
    not keep_all or args.append('--keep-all')
    if adapters:
        args += ['--adapters', adapters]
    not less_mem or args.append('--less-mem')
    not no_dereplicate or args.append('--no-dereplicate')
    not no_interleave or args.append('--no-interleave')
    not no_fasta_interleave or args.append('--no-fasta-interleave')
    if verbosity > DEFAULT_VERBOSITY:
        args.append('--verbosity={}'.format(verbosity))

    p = subprocess.run(
        [script] + args,
        cwd=str(path),
    )
    if p.check_returncode():
        raise RuntimeError('Failed to run qc-sample: sample: {}: exit status: '
                           '{}'.format(path, p.returncode))


def qc(samples, *, clean_only=False, adapters=None, keep_all=False,
       less_mem=False, no_dereplicate=False, no_interleave=False,
       no_fasta_interleave=False, verbosity=DEFAULT_VERBOSITY,
       threads=DEFAULT_THREADS, project=None):
    """
    Do quality control on multiple samples

    :param samples: List of pathlib.Path containing read data, one per sample
    :param kwargs: Options, see omics.qc.main for argsparse options.
    """
    if less_mem:
        # qc-sample script uses two threads internally
        threads = int(threads / 2)
    else:
        # qc-sample script uses four threads internally
        threads = int(threads / 4)
    threads = max(1, threads)

    errors = []

    with ThreadPoolExecutor(max_workers=threads) as e:
        futures = {}
        for path in samples:
            futures[e.submit(
                qc_sample,
                path,
                clean_only=clean_only,
                adapters=adapters,
                keep_all=keep_all,
                less_mem=less_mem,
                no_dereplicate=no_dereplicate,
                no_interleave=no_interleave,
                no_fasta_interleave=no_fasta_interleave,
                verbosity=verbosity,
                project=project,
            )] = path

        for fut in as_completed(futures.keys()):
            sample_path = futures[fut]
            e = fut.exception()
            if e is None:
                if verbosity >= DEFAULT_VERBOSITY + 2:
                    print('Done: {}'.format(sample_path.name))
            else:
                errors.append(
                    '[qc] (error): {}: {}: {}'
                    ''.format(sample_path, e.__class__.__name__, e)
                )

    if errors:
        raise(RuntimeError('\n'.join(errors)))


def check_result(*path):
    """
    Check results of FASTQC
    """
    for i in path:
        for d in ['fwd', 'rev']:
            result_file = \
                i / 'FASTQC' / '{}_derep_scythe_sickle_fastqc.zip'.format(d)
            with zipfile.ZipFile(result_file) as zf:
                summary_name = \
                    '{}_derep_scythe_sickle_fastqc/summary.txt'.format(d)
                with zf.open(summary_name) as summary:
                    for line in summary:
                        mark, test_name, _ = line.decode().split('\t')
                        if mark == 'WARN':
                            print(i, d, mark, ':', test_name)


def main():
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
        '--clean-only',
        action='store_true',
        help='Remove all files made by a previously run of qc and exit.',
    )
    argp.add_argument(
        '-a', '--adapters',
        metavar='FILE',
        help='Specify the adapters file used in the adpater trimming step.  By'
             ' default the Illumina adapter file TruSeq3-PE-2.fa as '
             'distributed by the Trimmomatic project will be used.',
    )
    argp.add_argument(
        '--keep-all',
        action='store_true',
        help='Keep all intermediate files, by default some not-so-important '
             'intermediate results will be deleted to save disk space.',
    )
    argp.add_argument(
        '--less-mem',
        action='store_true',
        help='This option will reduce the dominating memory requirements for '
             'the de-replication step by half, typically, and double the '
             'computation time.',
    )
    argp.add_argument(
        '--no-dereplicate',
        action='store_true',
        help='Option to skip the de-replication step',
    )
    argp.add_argument(
        '--no-interleave',
        action='store_true',
        help='Option to skip building the interleaved reads file',
    )
    argp.add_argument(
        '--no-fasta-interleave',
        action='store_true',
        help='Skip building the interleaved fasta file, interleaved fastq '
             'files will still be build.',
    )
    argp.add_argument(
        '--check',
        action='store_true',
        help='Check fastqc results',
    )
    args = argp.parse_args()
    args.samples = [Path(i) for i in args.samples]
    for i in args.samples:
        if not i.is_dir():
            argp.error('Directory not found: {}'.format(i))

    if args.check:
        check_result(*args.samples)
        return

    try:
        qc(
            args.samples,
            clean_only=args.clean_only,
            adapters=args.adapters,
            keep_all=args.keep_all,
            less_mem=args.less_mem,
            no_dereplicate=args.no_dereplicate,
            no_interleave=args.no_interleave,
            no_fasta_interleave=args.no_fasta_interleave,
            verbosity=args.verbosity,
            threads=args.threads,
            project=args.project,
        )
    except Exception as e:
        if args.traceback:
            raise
        else:
            print('[qc] (error) {}: {}'.format(e.__class__.__name__, e),
                  file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
