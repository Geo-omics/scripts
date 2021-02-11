# Copyright 2019 Regents of The University of Michigan.

# This file is part of geo-omics-scripts.

# Geo-omics-scripts is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# Geo-omics-scripts is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with Geo-omics-scripts.  If not, see <https://www.gnu.org/licenses/>.

"""
Find and remove replicated reads from fastq files.
"""

import argparse
from binascii import hexlify
from itertools import zip_longest
from pathlib import Path

from . import get_argparser, DEFAULT_VERBOSITY

ST_HEAD = 1
ST_SEQ = 2
ST_PLUS = 3
ST_SCORE = 4


def hash_read_pair(fwd_head, fwd_seq, rev_head='', rev_seq=''):
    """
    Hash function for reads, reverse read may be omitted

    Uses built-in hash of flowcell + lane + sequence
    """
    return b':'.join(
        fwd_head.strip().split(b':')[:4] +
        [fwd_seq.strip()] +
        rev_head.strip().split(b':')[:4] +
        [rev_seq.strip()]
    ).__hash__()


def read_groups(fwd, rev):
    """
    Iterate over reads

    :param fwd: File handle to forward reads
    :param rev: File handle to reverse reads

    This implements the grouper recipe on a zipped input.  The output
    iterates over 4-tuples of 2-tuples, i.e. the forward, reverse pairs
    of the four lines of a read:

      (@fwd_head, @rev_head), (fwd_seq, rev_seq), (+, +), (fwd_qual, rev_qual)
    """
    paired = zip_longest(fwd, rev)
    args = [iter(paired)] * 4
    return zip_longest(*args)


def find_duplicates(fwd_in, rev_in=None, *, check=False):
    """
    Find duplicate reads in given fastq files

    Calculates hash of concatenation of flowcell + lane part of header with
    whitespace-trimmed sequence string to detect replicated reads.  The
    file offsets of each header-@ of forward reads is recorded and semi-ordered
    such that the first offset listed belongs to the forward read of the
    highest quality paired duplicate read.
    """
    if rev_in is None:
        raise NotImplementedError('Single reads processing not implemented')

    data = {}
    fwd_read_pos = 0
    total_count = 0

    for (fh, rh), (fs, rs), (fp, rp), (fq, rq) in read_groups(fwd_in, rev_in):
        if check:
            if ord('>') in [fh[0], rh[0]]:
                raise NotImplementedError('Fasta support not implemented')
            if fh[0] != ord('@'):
                raise RuntimeError('Expected fastq header in {}: {}'
                                   ''.format(fwd_in.name, fh))
            if rh[0] != ord('@'):
                raise RuntimeError('Expected fastq header in {}: {}'
                                   ''.format(rev_in.name, rh))
            # TODO: check if headers match
            if not fp == rp == b'+\n':
                raise RuntimeError('Expected two + lines:\n{}\n{}'
                                   ''.format(fp, rp))

        paired_hash = hash_read_pair(fh, fs, rh, rs)

        # get mean score of concatenated quality with newlines
        cur_mean = mean_quality_score(fq + rq)

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

        # set positions for next read
        fwd_read_pos = fwd_in.tell()
        total_count += 1

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
        raise NotImplementedError('Single reads processing not implemented')

    pos = fwd_in.tell()
    for (fh, rh), (fs, rs), (fp, rp), (fq, rq) in read_groups(fwd_in, rev_in):
        if pos in refuse:
            if dupe_file is not None:
                hash_ = hash_read_pair(fh, fs, rh, rs)
                hash_ = hash_.to_bytes(length=8, byteorder='big', signed=True)
                hash_ = hexlify(hash_)
                dupe_file.write(fh.rstrip() + b'\t' + hash_ + b'\n')
        else:
            fwd_out.write(fh + fs + fp + fq)
            rev_out.write(rh + rs + rp + rq)

        if check:
            if not fh[0] == rh[0] == ord('@'):
                raise RuntimeError('Fastq header expected but found:\n{}\n{}'
                                   ''.format(fh, rh))

            if not fp == rp == b'+\n':
                raise RuntimeError('Plus separator expected but found:\n{}\n{}'
                                   ''.format(fp, rp))

        pos = fwd_in.tell()


def main(argv=None, namespace=None):
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
    args = argp.parse_args(args=argv, namespace=namespace)

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
