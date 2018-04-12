"""
Find and remove replicated reads from fastq files.
"""
import argparse
from binascii import hexlify
from pathlib import Path

from . import get_argparser, DEFAULT_VERBOSITY

ST_HEAD = 1
ST_SEQ = 2
ST_PLUS = 3
ST_SCORE = 4


def find_duplicates(fwd_in, rev_in=None, *, check=False):
    """
    Find duplicate reads in given fastq files

    Uses hash of sequence with newlines to detect replicated reads.  The
    file offsets of each header-@ of forward reads is recorded and semi-ordered
    such that the first offset listed belongs to the forward read of the
    highest quality paired duplicate read.
    """
    if rev_in is None:
        raise NotImplemented('Single reads processing not implemented')

    state = ST_HEAD
    data = {}
    fwd_read_pos = 0
    paired_hash = None
    total_count = 0
    for fwd_line, rev_line in zip(fwd_in, rev_in):
        if state is ST_HEAD:
            if check:
                if ord('>') in [fwd_line[0], rev_line[0]]:
                    raise NotImplementedError('Fasta support not implemented')
                if fwd_line[0] != ord('@'):
                    raise RuntimeError('Expected fastq header in {}: {}'
                                       ''.format(fwd_in.name, fwd_line))
                if rev_line[0] != ord('@'):
                    raise RuntimeError('Expected fastq header in {}: {}'
                                       ''.format(rev_in.name, rev_line))
            state = ST_SEQ
        elif state is ST_SEQ:
            paired_hash = (fwd_line + rev_line).__hash__()
            state = ST_PLUS
        elif state is ST_PLUS:
            if check and not fwd_line == rev_line == b'+\n':
                raise RuntimeError('Expected "+":\n{}\n'
                                   ''.format(fwd_line, rev_line))
            state = ST_SCORE
        elif state is ST_SCORE:
            # get mean score of concatenated quality with newlines
            cur_mean = mean_quality_score(fwd_line + rev_line)
            if paired_hash in data:
                best_mean, pos_list = data[paired_hash]
                if cur_mean > best_mean:
                    # new best read goes first
                    data[paired_hash] = (
                        cur_mean,
                        [fwd_read_pos] + pos_list
                    )
                else:
                    # read to be deleted, goes at end
                        data[paired_hash] = (
                            best_mean,
                            pos_list + [fwd_read_pos]
                        )
            else:
                data[paired_hash] = (
                    cur_mean,
                    [fwd_read_pos]
                )

            state = ST_HEAD
            # set positions for next read
            fwd_read_pos = fwd_in.tell()
            total_count += 1

        else:
            raise RuntimeError('Illegal state reached: {}'.format(state))

    return data, total_count


def mean_quality_score(score):
    """
    :param bytes score: Byte string encoding the score
    """
    return sum(score) / len(score)


def build_filter(data):
    """
    Transform data from hash-indexed to set of read positions
    """
    refuse = set()
    for _, pos_list in data.values():
        for i in pos_list[1:]:
            refuse.add(i)
    return refuse


def filter_write(refuse, fwd_in, rev_in, fwd_out, rev_out, check=False,
                 dupe_file=None):
    """
    Write out filtered data

    :param set refuse: Set of file offset positions of the read headers @ of
                       duplicated reads in the forward reads file.
    """
    if rev_in is None or rev_out is None:
        raise NotImplemented('Single reads processing not implemented')

    state = ST_HEAD
    pos = fwd_in.tell()
    for fwd_line, rev_line in zip(fwd_in, rev_in):
        if pos in refuse:
            if state is ST_HEAD and dupe_file is not None:
                dupe_file.write(fwd_line.rstrip() + b'\t')
            if state is ST_SEQ and dupe_file is not None:
                hash_ = (fwd_line + rev_line).__hash__()
                hash_ = hash_.to_bytes(length=8, byteorder='big', signed=True)
                hash_ = hexlify(hash_)
                dupe_file.write(hash_ + b'\n')
        else:
            fwd_out.write(fwd_line)
            rev_out.write(rev_line)

        if state is ST_HEAD:
            if check and not fwd_line[0] == rev_line[0] == ord('@'):
                raise RuntimeError('Fastq header expected but found:\n{}\n{}'
                                   ''.format(fwd_line, rev_line))
            state = ST_SEQ
        elif state is ST_SEQ:
            state = ST_PLUS
        elif state is ST_PLUS:
            if check and not fwd_line == rev_line == b'+\n':
                raise RuntimeError('Plus separator expected but found:\n{}\n{}'
                                   ''.format(fwd_line, rev_line))
            state = ST_SCORE
        elif state is ST_SCORE:
            state = ST_HEAD
            pos = fwd_in.tell()
        else:
            raise RuntimeError('Illegal internal state: {}'.format(state))


def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
        threads=False,
    )
    argp.add_argument('forward_reads', type=argparse.FileType())
    argp.add_argument('reverse_reads', type=argparse.FileType())
    argp.add_argument(
        '-c', '--check',
        action='store_true',
        help='Run sanity checks on input'
    )
    argp.add_argument(
        '-i', '--infix',
        default='.derep',
        help='Infix to construct output filenames',
    )
    argp.add_argument(
        '-o', '--out-dir',
        default='.',
        help='Output directory, by default this is the current working '
             'directory. The directory must exist.',
    )
    argp.add_argument(
        '--replicates-list',
        default=None,
        metavar='FILE',
        type=argparse.FileType('wb'),
        help='If provided, the list of replicated reads is written to the '
             'given file.',
    )
    args = argp.parse_args()

    out_dir = Path(args.out_dir)
    if not out_dir.is_dir():
        argp.error('Directory does not exist: {}'.out_dir)

    args.forward_reads.close()
    args.reverse_reads.close()

    fwd_path = Path(args.forward_reads.name)
    rev_path = Path(args.reverse_reads.name)

    fwd_in = fwd_path.open('rb')
    rev_in = rev_path.open('rb')

    data, total_reads = find_duplicates(fwd_in, rev_in, check=args.check)

    if args.verbosity > DEFAULT_VERBOSITY:
        print('total paired-read count: {}'.format(total_reads))

    refuse = build_filter(data)

    if args.verbosity > DEFAULT_VERBOSITY:
        print('replicated paired-reads:', len(refuse))

    fwd_in.seek(0)
    rev_in.seek(0)

    fwd_out_path = out_dir / (fwd_path.stem + args.infix + fwd_path.suffix)
    rev_out_path = out_dir / (rev_path.stem + args.infix + rev_path.suffix)

    fwd_out = fwd_out_path.open('wb')
    rev_out = rev_out_path.open('wb')

    if args.verbosity > DEFAULT_VERBOSITY:
        print('writing dereplicated output to {} and {} ...'
              ''.format(fwd_out_path, rev_out_path), end='', flush=True)

    filter_write(refuse, fwd_in, rev_in, fwd_out, rev_out,
                 dupe_file=args.replicates_list)
    fwd_out.close()
    rev_out.close()

    if args.verbosity > DEFAULT_VERBOSITY:
        print('done')


if __name__ == '__main__':
    main()
