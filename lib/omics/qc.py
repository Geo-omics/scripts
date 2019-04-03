"""
Run QC on metagenomic reads from multiple samples
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys

from omics import get_argparser, DEFAULT_VERBOSITY

FILTER_CHOICES = ['trimmomatic', 'scythe', 'bbtools']
QC_BINARY_NAME = 'omics-qc-sample'
DEFAULT_RQCFILTERDATA = '/reference-data/bbtools/RQCFilterData'


def qc_sample(path, **kwargs):
    """
    Do QC for a single sample

    This is a wrapper for the omics-qc-sample script
    """
    for k, v in vars(get_args(argv=[])).items():
        kwargs.setdefault(k, v)

    args = []
    if kwargs['fwd']:
        args += ['--fwd', kwargs['fwd']]
    if kwargs['rev']:
        args += ['--rev', kwargs['rev']]
    not kwargs['clean_only'] or args.append('--clean-only')
    not kwargs['keep_all'] or args.append('--keep-all')
    if kwargs['adapters']:
        args += ['--adapters', kwargs['adapters']]
    not kwargs['no_dereplicate'] or args.append('--no-dereplicate')
    not kwargs['no_fasta_interleave'] or args.append('--no-fasta-interleave')
    if kwargs['filter']:
        args += ['--filter', kwargs['filter']]
    if kwargs['rqcfilterdata'] is not None:
        args += ['--rqcfilterdata', kwargs['rqcfilterdata']]
    if kwargs['threads'] is not None:
        args += ['--threads', str(kwargs['threads'])]
    if kwargs['verbosity'] > DEFAULT_VERBOSITY:
        args.append('--verbosity={}'.format(kwargs['verbosity']))
        print('[qc] Calling qc-sample with arguments: {}'.format(args))

    p = subprocess.run(
        [QC_BINARY_NAME] + args,
        cwd=str(path),
    )
    if p.check_returncode():
        raise RuntimeError('Failed to run qc-sample: sample: {}: exit status: '
                           '{}'.format(path, p.returncode))


def qc(**kwargs):
    """
    Do quality control on multiple samples

    :param kwargs: Options, see omics.qc.main for argsparse options.
    """
    for k, v in vars(get_args(argv=[])).items():
        kwargs.setdefault(k, v)

    # number of workers: how many samples are processed in parallel
    # HOWTO distribute CPUs among workers:
    #   1. allow at most (hard-coded) 6 workers
    #   2. # of workers is the least of 6, # of CPUs, and # of samples
    #   3. # of CPUs/worker at least 1
    MAX_WORKERS = 1
    num_workers = min(MAX_WORKERS, len(kwargs['samples']), kwargs['threads'])
    threads_per_worker = max(1, int(kwargs['threads'] / num_workers))

    if kwargs['verbosity'] > DEFAULT_VERBOSITY:
        print('[qc] Processing {} samples in parallel, using {} threads each'
              ''.format(num_workers, threads_per_worker))

    kwargs['threads'] = threads_per_worker

    errors = []
    with ThreadPoolExecutor(max_workers=num_workers) as e:
        futures = {}
        for path in kwargs['samples']:
            futures[e.submit(qc_sample, path, **kwargs)] = path

        for fut in as_completed(futures.keys()):
            sample_path = futures[fut]
            e = fut.exception()
            if e is None:
                if kwargs['verbosity'] >= DEFAULT_VERBOSITY + 2:
                    print('Done: {}'.format(sample_path.name))
            else:
                errors.append(
                    '[qc] (error): {}: {}: {}'
                    ''.format(sample_path, e.__class__.__name__, e)
                )

    if errors:
        raise(RuntimeError('\n'.join(errors)))


def get_args(argv=None, namespace=None):
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
        '-f', '--fwd',
        help='Fastq file name with forward reads, default is fwd.fastq.  A '
             'file with this name needs to be present in each sample '
             'directory.',
    )
    argp.add_argument(
        '-r', '--rev',
        help='Fastq file name with reverse reads, default is rev.fastq.  A '
             'file with this name needs to be present in each sample '
             'directory.',
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
             'distributed by the Trimmomatic project will be used.  Ignored '
             'when using the bbtools filter',
    )
    argp.add_argument(
        '--keep-all',
        action='store_true',
        help='Keep all intermediate files, by default some not-so-important '
             'intermediate results will be deleted to save disk space.',
    )
    argp.add_argument(
        '--no-dereplicate',
        action='store_true',
        help='Option to skip the de-replication step',
    )
    argp.add_argument(
        '--no-fasta-interleave',
        action='store_true',
        help='Skip building the interleaved fasta file, interleaved fastq '
             'files will still be build.',
    )
    argp.add_argument(
        '-F', '--filter',
        choices=FILTER_CHOICES,
        default=FILTER_CHOICES[0],
        help='Which quality filter to use.  Trimmomatic is used by default. '
             'Scythe/sickle based filtering the the "old" way.  The rqcfilter2'
             ' pipeline from BBTools is an alternative',
    )
    argp.add_argument(
        '--rqcfilterdata',
        metavar='PATH',
        help='Path to directory containing the reference data for BBTool\'s '
             'rqcfilter2. The default is ' + DEFAULT_RQCFILTERDATA + ' and '
             'is good for running omics qc inside the omics container.'
    )
    args = argp.parse_args(args=argv, namespace=namespace)
    args.samples = [Path(i) for i in args.samples]
    for i in args.samples:
        if not i.is_dir():
            argp.error('Directory not found: {}'.format(i))
    return args


def main(argv=None, namespace=None):
    args = get_args(argv, namespace)
    try:
        qc(**vars(args))
    except Exception as e:
        if args.traceback:
            raise
        else:
            print('[qc] (error) {}: {}'.format(e.__class__.__name__, e),
                  file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
