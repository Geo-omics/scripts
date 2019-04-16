"""
Prepare fastq files for processing with Geomicro Illumina Reads Pipeline
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
from itertools import groupby, zip_longest
from pathlib import Path
import re
import shutil
import sys

from . import get_argparser, DEFAULT_VERBOSITY

READ_COUNT_FILE_NAME = 'read_count.tsv'

FORWARD_READS_FILE = 'fwd.fastq'
REVERSE_READS_FILE = 'rev.fastq'

raw_reads_file_pat = re.compile(
    r'(?P<sampleid>[^_]+)(_(?P<index>[-a-zA-Z]+))?(_S(?P<snum>\d+))?'
    r'_L(?P<lane>\d+)_R(?P<dir>[12])_(?P<fnum>\d+)\.(?P<suffix>.*)'
)

# Multi-run file names, based on example from Matt
# no index here
multi_run_file_pat = re.compile(
    r'(?P<runid>[^_]+)_(?P<sampleid>[^_]+)(_(?P<index>[-a-zA-Z]+))?'
    r'(_S(?P<snum>\d+))?'
    r'_L(?P<lane>\d+)_R(?P<dir>[12])_(?P<fnum>\d+)\.(?P<suffix>.*)'
)


class FileNameDoesNotMatch(Exception):
    pass


def group(files, keep_lanes=False, multi_run=False, verbosity=1):
    """
    Generate groups of raw reads files per sample

    :param list files: Files, list of pathlib.Path objects
    :param bool keep_lanes: If true, then data from one sample that was spread
                            across several lanes are not concatenated and kept
                            separate.
    :param bool multi_run:  Use multi-run file patters, and consider runid in
                            sorting.

    :return: Iterator over tuples containing the sort key and the group.  The
             key comes as a tuple of sub-keys, i.e. the sample id and if
             ``keep_lanes == True`` also the lane.
             The group is an iterator over the files with the given keys.

    Filenames are assumed to adhere to the usual scheme:

        <identifier>[_<index>][_S<nnn>]_L<nnn>_R<1|2>_<nnn>.fastq[.gz]

        e.g. 66145_CATTGAC_S1_L007_R1_001.fastq.gz

    but graceful degradation is attempted in case the files do not correspond
    to this scheme.
    """
    if multi_run:
        pattern = multi_run_file_pat
    else:
        pattern = raw_reads_file_pat

    def id_lane_dir(x):
        """
        key function for initial sorting

        Primary key is sample name.
        If keep_lanes then lanes come before direction, otherwise they're last.
        In multi_run mode the run id is treated as (the primary) part of the
        lane.

        Degrade to identity if filename can't be parsed.
        """
        m = re.match(pattern, x.name)
        if m is None:
            key = (x.name, 0, '')
        else:
            m = m.groupdict()
            key = (m['sampleid'],)

            if multi_run:
                runid = m['runid']
            else:
                runid = None

            if keep_lanes:
                key += (runid, int(m['lane']), m['dir'])
            else:
                key += (m['dir'], runid, int(m['lane']))

        if verbosity >= DEFAULT_VERBOSITY + 2:
            print('parsed: {} -> {}'.format(x, key))

        return key

    files = sorted(files, key=id_lane_dir)

    def id_lane(x):
        """
        key function for grouping

        Key is sample name and if we keep lanes also the lane
        """
        m = re.match(pattern, x.name)
        try:
            m = m.groupdict()
        except AttributeError:
            raise FileNameDoesNotMatch(x)

        key = (m['sampleid'],)
        if keep_lanes:
            if multi_run:
                # need runid to identify lanes since runs may use same lane
                key += (m['runid'],)
            key += (m['lane'],)
        return key

    for key, group in groupby(files, key=id_lane):
        yield (key, group)


def without_filenumber(path):
    """
    Helper function to remove the file number from a path

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


def sample_direction(filename):
    """
    Extract the read direction from filename

    :param filename: Path or (partial) filename, can be str or Path
    :return: Integer 1 or 2
    :raise: AttributeError if parsing fails

    This is somewhat more permissive than the other filename parsing in prep
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

    fwd_outfile = destdir / FORWARD_READS_FILE
    rev_outfile = destdir / REVERSE_READS_FILE
    for i in [fwd_outfile, rev_outfile]:
        if i.is_file():
            if force:
                with i.open('wb'):
                    # truncating
                    pass
            else:
                raise FileExistsError(i)

    for direction, series in groupby(files, key=sample_direction):
        # a 'series' is a bunch of files that got split up, and that we need
        # to put back together in a consistent order
        if direction == 1:
            outfile = fwd_outfile
        elif direction == 2:
            outfile = rev_outfile
        else:
            raise ValueError('Illegal value for direction: {}'
                             ''.format(direction))
        series = list(series)

        args = (sample, outfile, series, verbosity)
        if executor is None:
            _do_extract_and_copy(*args)
        else:
            futures[
                executor.submit(_do_extract_and_copy, *args)
            ] = (sample, direction)

    return futures


def _do_extract_and_copy(sample, outfile, series, verbosity):
    """
    Helper function doing all the parallelizable IO work

    :param Path outfile: Output file
    :param series: List of input files

    Returns name of output file
    """
    with outfile.open('ab') as outf:
        for i in series:
            if i.suffix == '.gz':
                infile = gzip.open(str(i), 'r')
                action = 'extr'
            else:
                infile = i.open('rb')
                action = 'copy'

            if verbosity > DEFAULT_VERBOSITY:
                try:
                    dest = outfile.resolve().relative_to(Path.cwd())
                except:
                    dest = outfile
                print('{}: {} {} >> {}'.format(sample, action, i,
                      dest))

            try:
                shutil.copyfileobj(infile, outf, 4 * 1024 * 1024)

            finally:
                infile.close()

    return outf.name


def count_fastq_reads(path, verbose=False):
    """
    Count number of reads in fastq file
    """
    count = 0
    args = [iter(path.open('rb'))] * 4
    if verbose:
        print('Start counting reads for {}...'.format(path))
    for read in zip_longest(*args):
        if read[-1] is None:
            raise RuntimeError(
                'Line count is not a multiple of 4: {}, read count at {}\n'
                'Last lines are:\n{}'
                ''.format(path, count, '\n'.join([str(i) for i in read]))
            )
        count += 1
    return count


def main(argv=None):
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__
    )
    argp.add_argument(
        'rawreads',
        nargs='*',
        default=['.'],
        help='Your raw fastq or fastq.gz files or path to directory containing'
             ' them. By default this is the current directory.  If this is a '
             'directory the script will respect the --suffix option and try '
             'to copy all the files with the given suffices.  If any of those '
             'files\'s name does not follow the support Illumina pattern '
             'then the script will abort. In this case you need to list the '
             'validly named file explicitly on the command line.',
    )
    argp.add_argument(
        '-d', '--dest',
        metavar='PATH',
        default='.',
        help='Destination directory, by default the current working directory',
    )
    argp.add_argument(
        '--count-reads',
        action='store_true',
        help='Make simple read-count statistics',
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
        '--multi-run',
        action='store_true',
        help='Assume data is from multiple sequencing runs.  File name are '
             'interpreted such that the first two underscore-separated parts '
             'are run and sample identifiers.  '
             'Obviously the different runs must use the same sample ids for '
             'this to work.'
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
    args = argp.parse_args(argv)

    dest = Path(args.dest)
    if not dest.is_dir():
        argp.exit('Not a directory: {}'.format(args.dest))

    verbosity = args.verbosity
    quiet = verbosity < DEFAULT_VERBOSITY
    verbose = verbosity > DEFAULT_VERBOSITY
    very_verbose = verbosity > DEFAULT_VERBOSITY + 1

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
        if not quiet:
            print('Found {} read files.'.format(len(files)))
        if very_verbose:
            for i in files:
                print('  ->', i)
    else:
        argp.error('No files found.')

    if very_verbose:
        print('Using {} threads.'.format(args.threads))

    files = list(set(files))
    samp_count = 0
    read_counts = {}

    try:
        with ThreadPoolExecutor(max_workers=args.threads) as e:
            futures = {}
            file_groupings = group(files, keep_lanes=args.keep_lanes,
                                   multi_run=args.multi_run,
                                   verbosity=verbosity)
            for samp_key, samp_grp in file_groupings:
                samp_count += 1

                samp_key = list(samp_key)
                if len(samp_key) > 1:
                    try:
                        # remove leading zeros
                        samp_key[1] = str(int(samp_key[1]))
                    except Exception:
                        pass

                futures.update(
                    prep(
                        '_'.join(samp_key),
                        list(samp_grp),
                        dest=dest,
                        force=args.force,
                        verbosity=verbosity,
                        executor=e,
                    )
                )

            while futures:
                for fut in as_completed(futures.keys()):
                    job_meta = futures[fut]
                    if len(job_meta) == 1:
                        # read counting job, job_meta is 1-tuple
                        sample = futures[fut][0]
                        count = fut.result()
                        if sample in read_counts:
                            raise RuntimeError(
                                'Already counted reads for {} (internal error)'
                                ''.format(sample)
                            )
                        else:
                            read_counts[sample] = count
                            if verbose:
                                print('{}: {} reads'.format(sample, count))
                    elif len(job_meta) == 2:
                        # was an extract/copy job, job_meta is 2-tuple
                        sample, direction = futures[fut]
                        outfile = Path(fut.result())
                        if very_verbose:
                            print(
                                'Done: {} {}'.format(
                                    sample,
                                    'fwd' if direction == 1 else 'rev'
                                ),
                                end=''
                            )
                        if args.count_reads and direction == 1:
                            futures[e.submit(
                                count_fastq_reads,
                                outfile,
                                verbose=verbose,
                            )] = (sample, )
                        if fut.exception() is not None:
                            print('Failed to write: {}: {}: {}'
                                  ''.format(sample, direction, outfile),
                                  file=sys.stderr)
                    else:
                        raise RuntimeError(
                            '(internal error) invalid futures structure: '
                            '{}: {}: {}'
                            ''.format(e.__class__.__name__, e, job_meta)
                        )
                    del futures[fut]

    except FileNameDoesNotMatch as e:
        print('The name of file {} does not follow the supported pattern:\n'
              '  <sample>_<index>_S<n>_L<nnn>_R(1|2)_<nnn>.fastq[.gz]\n'
              'or with --multi-run\n'
              '  <run>_<sample>_<index>_S<n>_L<nnn>_R(1|2)_<nnn>.fastq[.gz]\n'
              'You may need to manually set up your sample '
              'directories.'.format(e), file=sys.stderr)
        sys.exit(1)

    except Exception as e:
        if args.traceback:
            raise
        else:
            print('{}: {}'.format(e.__class__.__name__, e), file=sys.stderr)
            sys.exit(1)

    with open(READ_COUNT_FILE_NAME, 'w') as f:
        for sample, count in sorted(read_counts.items()):
            out = '{}\t{}\n'.format(sample, count)
            f.write(out)
            if verbose:
                print(out, end='')

    if verbose:
        print('read counts written to', READ_COUNT_FILE_NAME)

    if not quiet:
        print('Processed {} samples'.format(samp_count))


if __name__ == '__main__':
    main()
