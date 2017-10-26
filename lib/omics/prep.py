"""
Prepare fastq files for processing with Geomicro Illumina Reads Pipeline
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
from itertools import groupby
from pathlib import Path
import re
import shutil
import sys

from . import get_argparser, get_project, DEFAULT_VERBOSITY

raw_reads_file_pat = re.compile(
    r'(?P<sampleid>[^_]+)_((?P<index>[a-zA-Z]+)_)?S(?P<snum>\d+)_'
    r'L(?P<lane>\d+)_R(?P<dir>[12])_(?P<fnum>\d+)\.(?P<suffix>.*)'
)


def group(files, keep_lanes=False):
    """
    Generate groups of raw reads files per sample

    :param list files: Files, list of pathlib.Path objects
    :param bool keep_lanes: If true, then data from one sample that was spread
                            across several lanes are concatenated.

    :return: Iterator over groups (which are also iterators)

    Filenames are assumed to adhere to the usual scheme:

        <identifier>[_<index>]_S<nnn>_L<nnn>_R<1|2>_<nnn>.fastq[.gz]

        e.g. 66145_CATTGAC_S1_L007_R1_001.fastq.gz

    but graceful degradation is attempted in case the files do not correspond
    to this scheme.
    """
    files = sorted(files, key=lambda x: x.name)

    def keyfunc(x):
        """ key function for grouping """
        # get sample id, all up to first underscore
        key = x.name.partition('_')[0]
        if keep_lanes:
            # try to obtain the lane
            m = re.search('_L\d+', x.name)
            try:
                lane = m.group()
            except AttributeError:
                raise RuntimeError(
                    'Failed extract lane info from filename: {}'.format(x)
                )
            else:
                key += lane
        return key

    for key, group in groupby(files, keyfunc):
        yield (key, group)


def without_filenumber(path):
    """
    Helper function to removethe file number from a path

    Does this:
        id_Snnn_Lnnn_Rn_nnn.suf --> id_Snnn_Lnnn_Rn

    This leave the stem, identifying a series of read file that where split up
    Use this as a key function for grouping.  If paring fails the original
    path is returned.
    """
    pre_stem = path.name.partition('.')[0]  # rm suffices
    stem, _, fnum = pre_stem.rpartition('_')
    try:
        int(fnum)
    except ValueError:
        # it's not the file number -> degrade gracefully
        return path

    if stem:
        # fnum is a number and stem is what we think it is
        return stem + '_'  # typically the _ is still needed to parse direction

    # rpartition failed, stem is empty -> degrade gracefully
    return path


def get_direction(filename):
    """
    Extract the read direction from filename

    :param filename: Path or (partial) filename, can be str or Path
    :return: Integer 1 or 2
    :raise: AttributeError if parsing fails
    """
    filename = Path(filename).name
    m = re.search(r'_R(?P<dir>[12])_', filename)
    return int(m.groupdict()['dir'])


def prep(sample, files, dest=Path.cwd(), force=False, verbosity=1,
         executor=None):
    """
    Decompress and copy files into sample directory

    :param str sample: Sample identifier
    :param list files: List of pathlib.Path objects, must be sorted by filename
                       such that files belonging to the same series (i.e.
                       differing by 'filenumber' are grouped together.
    :param Path dest: Project directory, destination / output directory
    :param bool force: If true, skip any safety check and overwrite existing
                       files.
    :param executor: A concurrent.futures.Executor object.  If executor is None
                     then all will be done single-threaded.

    :return: Dictionary of futures

    :raise: May raise OSError and relatives, in particular when destination
            files exist already or can not be written.

    """
    futures = {}

    destdir = dest / sample
    destdir.mkdir(exist_ok=True, parents=True)

    for stem, series in groupby(files, without_filenumber):
        # a 'series' is a bunch of files from the same lane/sample that got
        # split up, and that we need to put back together.  It's not clear if
        # such data is still produced in the wild.
        try:
            direction = get_direction(stem)
        except Exception as e:
            raise RuntimeError(
                'Failed to parse read direction from filename: {} in {}'
                ''.format(stem, files)
            )
        if direction == 1:
            outname = 'fwd.fastq'
        elif direction == 2:
            outname = 'rev.fastq'
        else:
            raise ValueError('Illegal value for direction: {}'
                             ''.format(direction))
        series = list(series)
        outfile = destdir / outname
        if outfile.is_file() and not force:
            raise FileExistsError(outfile)

        if executor is None:
            _do_extract_and_copy(outfile, series)
        else:
            futures[
                executor.submit(
                    _do_extract_and_copy,
                    sample,
                    outfile,
                    series,
                    verbosity
                )
            ] = (sample, direction)

    return futures


def _do_extract_and_copy(sample, outfile, series, verbosity):
    """
    Helper function doing all the parallelizable IO work

    :param Path outfile: Output file
    :param series: List of input files
    """
    with outfile.open('wb') as outf:
        for i in series:
            if i.suffix == '.gz':
                infile = gzip.open(str(i), 'r')
                action = 'extr'
            else:
                infile = i.open('rb')
                action = 'copy'

            if verbosity > DEFAULT_VERBOSITY:
                print('{}: {} {} >> {}'.format(sample, action, i,
                      outfile.relative_to(Path.cwd())))

            try:
                shutil.copyfileobj(infile, outf, 4 * 1024 * 1024)
            finally:
                infile.close()


def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__
    )
    argp.add_argument(
        'rawreads',
        nargs='*',
        default=['.'],
        help='List of fastq.gz files or path to directory containing them. '
             'By default this is the current directory.',
    )
    argp.add_argument(
        '--force', '-f',
        action='store_true',
        help='Overwrite existing files',
    )
    argp.add_argument(
        '--keep-lanes-separate',
        action='store_true',
        dest='keep_lanes',
        help='Keep data from different lanes separate.  The default is to '
             'collect reads originating from the same physical sample if '
             'sequencing was done using several lanes.',
    )
    argp.add_argument(
        '--suffix',
        metavar='LIST',
        default='fastq,fastq.gz',
        help='Comma-separated list of valid file suffices used for raw reads.'
             ' This is used to find files when a directory is given as '
             'positional argument.  By default .fastq and .fastq.gz files '
             'are considered.',
    )
    argp.add_argument(
        '--cpus', '--threads', '-t',
        type=int,
        dest='threads',
        default=1,
        help='Number of threads / CPUs to employ',
    )
    args = argp.parse_args()

    project = get_project(args.project_home)

    if args.verbosity == DEFAULT_VERBOSITY:
        verbosity = project['verbosity']
    else:
        verbosity = args.verbosity

    suffices = args.suffix.split(',')
    files = []
    for i in map(Path, args.rawreads):
        if i.is_dir():
            for suf in suffices:
                suf = suf.lstrip('.')
                files += i.glob('*.' + suf)
        elif i.is_file():
            files.append(i)
        else:
            argp.error('File or directory not found: {}'.format(i))

    if files:
        if verbosity >= DEFAULT_VERBOSITY:
            print('Found {} read files.'.format(len(files)))
        if verbosity >= DEFAULT_VERBOSITY + 2:
            for i in files:
                print('  ->', i)
    else:
        argp.error('No files found.')

    files = list(set(files))
    samp_count = 0

    try:
        with ThreadPoolExecutor(max_workers=args.threads) as e:
            futures = {}
            for sample, sample_grp in group(files, keep_lanes=args.keep_lanes):
                samp_count += 1
                futures.update(
                    prep(
                        sample,
                        list(sample_grp),
                        dest=project['project_home'],
                        force=args.force,
                        verbosity=verbosity,
                        executor=e,
                    )
                )

            for fut in as_completed(futures.keys()):
                if fut.exception() is not None:
                    print(fut.result())

    except Exception as e:
        if args.traceback:
            raise
        else:
            print('{}: {}'.format(e.__class__.__name__, e), file=sys.stderr)
            sys.exit(1)

    if verbosity >= 1:
        print('Processed {} samples'.format(samp_count))


if __name__ == '__main__':
    main()
