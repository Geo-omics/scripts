import argparse
from pathlib import Path

from . import get_argparser

ST_HEAD = 1
ST_SEQ = 2
ST_PLUS = 3
ST_SCORE = 4


def find_duplicates(fwd_path, rev_path=None, *, check=False):
    """
    Find duplicate reads in given fastq files
    """
    if rev_path is None:
        raise NotImplemented('Single reads processing not implemented')

    fwd_f = fwd_path.open('rb')
    rev_f = rev_path.open('rb')

    state = ST_HEAD
    data = {}
    fwd_read_pos = 0
    rev_read_pos = 0
    paired_hash = None
    total_count = 0
    for fwd_line, rev_line in zip(fwd_f, rev_f):
        if state is ST_HEAD:
            if check:
                if ord('>') in [fwd_line[0], rev_line[0]]:
                    raise NotImplementedError('Fasta support not implemented')
                if fwd_line[0] != ord('@'):
                    raise RuntimeError('Expected fastq header in {}: {}'
                                       ''.format(fwd_path, fwd_line))
                if fwd_line[0] != ord('@'):
                    raise RuntimeError('Expected fastq header in {}: {}'
                                       ''.format(fwd_path, fwd_line))
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
            cur_mean = mean_quality_score(fwd_line + rev_line)
            if paired_hash in data:
                best_mean, pos_list = data[paired_hash]
                if cur_mean > best_mean:
                    # new best read goes first
                    data[paired_hash] = (
                        cur_mean,
                        [(fwd_read_pos, rev_read_pos)] + pos_list
                    )
                else:
                    # read to be deleted, goes at end
                        data[paired_hash] = (
                            best_mean,
                            pos_list + [(fwd_read_pos, rev_read_pos)]
                        )
            else:
                data[paired_hash] = (
                    cur_mean,
                    [(fwd_read_pos, rev_read_pos)]
                )

            state = ST_HEAD
            # set positions for next read
            fwd_read_pos = fwd_f.tell()
            rev_read_pos = rev_f.tell()
            total_count += 1

        else:
            raise RuntimeError('Illegal state reached: {}'.format(state))

    return data, total_count


def mean_quality_score(score):
    """
    :param bytes score: Byte string encoding the score
    """
    return sum(score) / len(score)


def unique_pairs(fwd_data, rev_data):
    if len(fwd_data) != len(rev_data):
        raise RuntimeError('Sanity check failed, fwd vs. rev length difference:'
                           '{} != {}'.format(len(fwd_data), len(rev_data)))

    for (fwd_seq, fwd_reads), (rev_seq, rev_reads) \
            in zip(fwd_data.items(), rev_data.items()):

        if len(fwd_reads) != len(rev_reads):
            raise RuntimeError(
                'Sanity check failed: different number of replicates fwd vs. '
                'rev:\n{}: {}\n{}: {}'
                ''.format(fwd_seq, fwd_reads, rev_seq, rev_reads)
            )

        if len(fwd_reads) > 1:
            best_mean = -1
            best_pair = ()
            for f, r in zip(fwd_reads, rev_reads):
                mean = mean_quality_score(f) + mean_quality_score(r)
                if best_mean < mean:
                    best_mean = mean
                    best_pair = (f, r)
            yield f, r
        else:
            yield fwd_reads[0], rev_reads[0]


def write_read(read, file):
    file.write(
        read.header + '\n'
        + read.seq + '\n+\n'
        + read.score + '\n'
    )


def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__
    )
    argp.add_argument('forward_reads', type=argparse.FileType())
    argp.add_argument('reverse_reads', type=argparse.FileType())
    argp.add_argument(
        '-c', '--check',
        action='store_true',
        help='Run sanity checks on input'
    )
    argp.add_argument(
        '-o', '--out-prefix',
        default='derep_',
        help='Prefix attached to output files',
    )
    args = argp.parse_args()

    args.forward_reads.close()
    args.reverse_reads.close()

    fwd_path = Path(args.forward_reads.name)
    rev_path = Path(args.reverse_reads.name)

    data, total_reads = find_duplicates(fwd_path, rev_path, check=args.check)

    print('total read count: {}, duplicates: {}'
          ''.format(total_reads, total_reads - len(data)))

    fwd_out = fwd_reads.parent / (args.out_prefix + fwd_reads.name)
    rev_out = rev_reads.parent / (args.out_prefix + rev_reads.name)

    with fwd_out.open('w') as fout, rev_out.open('w') as rout:

        for fwd_read, rev_read in unique_pairs(fwd_data, rev_data):
            write_read(fwd_read, fout)
            write_read(rev_read, rout)

    #print(fwd, rev)


if __name__ == '__main__':
    main()
