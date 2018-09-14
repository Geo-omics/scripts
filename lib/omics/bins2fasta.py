"""
Copy contigs into per-bin FASTA files

Takes bin output table from CONCOCT or MetaBAT
"""
import argparse
from collections import Counter
import csv
from itertools import combinations
from pathlib import Path
from warnings import warn


# modes of operation
CONCOCT = 'CONCOCT'
METABAT = 'MetaBAT'


def get_argp():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'bintable',
        metavar='bin-table',
        type=argparse.FileType(),
        help='tab- or comma-separated, two-column table, mapping contigs to '
             'bins.  For CONCOCT this is a file usually named '
             'clustering_gt1000.csv and for MetaBAT with the --saveCls option '
             'a file named bins',
    )
    argp.add_argument(
        'real_assembly',
        metavar='real-assembly',
        type=argparse.FileType(),
        help='Filename of original, non-chopped assembly',
    )
    argp.add_argument(
        'chopped_assembly',
        nargs='?',
        metavar='chopped-assembly',
        type=argparse.FileType(),
        help='Filename of chopped assembly',
    )

    argp.add_argument(
        '-s', '--suffix',
        default='fna',
        help='Suffix for output files, default is .fna'
    )
    argp.add_argument(
        '-o', '--output-dir',
        default='.',
        help='Output directory, default is the current directory.'
    )
    argp.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print more info in case of errors'
    )
    return argp


def main():
    args = get_argp().parse_args()
    outdir = Path(args.output_dir)

    if args.verbose:
        print('Loading chunk info... ', end='', flush=True)

    cdata = load_contig_chunk_info(args.chopped_assembly or args.real_assembly)

    if args.verbose:
        print('done')
        print('Reading contig-to-bin mapping... ', end='', flush=True)

    mode = load_bins(cdata, args.bintable)

    if args.verbose:
        print('({} mode) done'.format(mode))
        print('Processing data... ', end='', flush=True)

    good_data, unbinned, ambiguous, cross = apply_policy(cdata)

    if args.verbose:
        print('done')
        print('Summary:')
        print('  total  contigs:', len(cdata.keys()))
        print('  binned contigs:', len(good_data.keys()))
        print('  ambiguous binnings:', len(ambiguous.keys()))
        print('  cross binnings:', len(cross))

    if args.verbose:
        print('Writing fasta files for bins... ', end='', flush=True)

    if args.real_assembly.closed:
        # when chopped asm was not given, file would be used and closed earlier
        args.real_assembly = open(args.real_assembly.name)

    files_written = bin_fasta(
        outdir=outdir,
        mapping=good_data,
        input_contigs=args.real_assembly,
        suffix=args.suffix,
        verbose=args.verbose,
    )

    if args.verbose:
        print(files_written, 'files written')


def parse_contig_chunk(s):
    """
    parse k255_99.7 -> ('k255_99', 7)
    """
    cont_id = s.strip().rsplit('.', 1)
    contig = cont_id[0]
    if len(cont_id) == 1:
        # no chunk info, keep as-is
        chunk = None
    elif len(cont_id) == 2:
        try:
            chunk = int(cont_id[1])
        except ValueError:
            raise RuntimeError('Failed to parse int chunk id '
                               '{}'.format(s))
    else:
        raise RuntimeError('Failed to parse contig + chunk in '
                           '{}'.format(s))
    return contig, chunk


def load_contig_chunk_info(file):
    """
    Read chop/chunk-aware contig ids from a fasta file

    :param Path file: Path to fasta formatted file with chopped contigs

    Returns empty contig-chunks-to-bin mapping, i.e. a dict of contig id to a
    list of bins, on bin per chunk, all initialized to None since the binning
    info has not been loaded yet.
    """
    data = {}
    for line in file:
        if line.startswith('>'):
            line = line.lstrip('>')
            contig_id = line.split()[0]
            try:
                contig, chunk = parse_contig_chunk(contig_id)
            except Exception as e:
                raise RuntimeError(
                    'Failed to parse contig id in {}: {}: {}: {}'
                    ''.format(file.name, line, e.__class__.__name__, e)
                )

            if contig in data:
                data[contig]
                data[contig].append(chunk)
            else:
                # initialize
                data[contig] = [chunk]

    file.close()

    # check for inconsistencies
    for contig, chunks in data.items():
        if chunks == [None]:
            continue
        elif sorted(chunks) == list(range(len(chunks))):
            continue
        else:
            raise RuntimeError('Non-consequtive chunks for contig {} found in '
                               '{}: {}'.format(contig, file, chunks))

    # return data in contig-to-chunk-bin-list format
    return {contig: [None] * len(chunks) for contig, chunks in data.items()}


def load_bins(cdata, file):
    """
    Load data from bin table

    :param dict cdata:  Initialized dictionary, mapping expected contig
                        chunks to bins
    :param file: file object with Binning table

    The input file must be a two-column table mapping contigs to bins
    The separator will be auto-detected.

    Changes cdata in-place.
    """

    # detect file format
    peeked = file.buffer.peek()
    if b',' in peeked:
        mode = CONCOCT
        sep = ','
    elif b'\t' in peeked:
        mode = METABAT
        sep = '\t'
    else:
        raise RuntimeError('Failed to detect file format of {}'
                           ''.format(file.name))

    try:
        for c, b in csv.reader(file, delimiter=sep):
            contig, chunk = parse_contig_chunk(c)
            try:
                b = int(b)
            except ValueError:
                raise RuntimeError('Integer expected for bin: {}{}{}'
                                   ''.format(c, sep, b))

            if mode == METABAT and b == 0:
                # skip, bin 0 in MetaBAT means contig is not binned
                continue

            # assign bin to chunk, if chunk is None also use 0 index
            cdata[contig][chunk or 0] = b

    except Exception as e:
        raise RuntimeError('Failed to parse binning result file {}: {}: {}'
                           ''.format(file, e.__class__.__name__, e))

    return mode


def apply_policy(cdata):
    """
    Determine final binning according to policy.
    """
    good = {}  # maps good contigs to bin
    unbinned = []  # list of unbinned / rejected contigs
    bad_ambiguous = {}  # maps unbinned / rejected contig's chunks to bins
    cross_binning = []  # list of edges of cross-binning graph

    for contig, bins in cdata.items():
        # POLICY: this and the following test implements the following policy:
        #
        #   A contig is binned to the bin that more than half of
        #   its chunks (by chunk count, not by sequence length)
        #   are binned to.  If no such absolute majority exists,
        #   then the contig remains unbinned.
        #

        # get most common bin and count
        _bin, count = sorted(Counter(bins).items(), key=lambda x: x[1])[-1]

        if _bin is not None and count / len(bins) > 0.5:
            # the is an absolute majority bin
            is_good = True
        else:
            is_good = False

        # compile outputs

        # get sorted list of bins for this bad contig
        chunk_bins = sorted({i for i in bins if i is not None})

        if is_good:
            good[contig] = _bin
        else:
            unbinned.append(contig)

            if chunk_bins:
                # some chunks were binned
                bad_ambiguous[contig] = bins

        if len(chunk_bins) > 1:
            # it's an edge in the cross-binning-graph
            for b1, b2 in combinations(chunk_bins, 2):
                cross_binning.append((b1, contig, b2))

    return good, unbinned, bad_ambiguous, cross_binning


def get_fa_sequences(file, contigs):
    """
    Generate fasta sequences from a fasta formatted file

    :param file:  File-like object, fasta-formatted text
    :param contig: Collection of names of contigs that should be yielded
    """
    name = None
    seq = None
    for line in file:
        if line.startswith('>'):
            if seq is not None:
                yield name, seq

            # start new sequence
            # sequence/contig name is text between > and first whitespace
            name = line.split()[0][1:]
            if name in contigs:
                seq = line
            else:
                # skip this sequence
                name = None
                seq = None
        else:
            if seq is not None:
                seq += line

    # last sequence or empty file
    if seq is not None:
        yield name, seq


def bin_fasta(outdir, mapping, input_contigs, suffix='fna', verbose=False):
    """
    Split fasta file into binned output
    """
    outdir.mkdir(parents=True, exist_ok=True)

    error = None
    files = {}  # map of bin to output file handle
    try:
        for contig, seq_txt in get_fa_sequences(input_contigs, mapping.keys()):
            _bin = mapping[contig]
            if _bin not in files:
                outfile = outdir / 'bin_{}.{}'.format(_bin, suffix)
                files[_bin] = outfile.open('w')

            files[_bin].write(seq_txt)
    except Exception as e:
        # save exception for after closing files
        error = e
    finally:
        for i in files.values():
            try:
                i.close()
            except Exception:
                warn('Failed to close file {}'.format(i.name))
        if error:
            raise error

    return len(files)


if __name__ == '__main__':
    main()
