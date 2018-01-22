import argparse
from collections import OrderedDict
from enum import Enum
from pathlib import Path

from . import get_argparser

ST_HEAD = 1
ST_SEQ = 2
ST_PLUS = 3
ST_SCORE = 4

def load(path):
    state = ST_HEAD
    data = OrderedDict()
    cur = {}
    total_count = 0
    for line in path.open('rb'):
        if state is ST_HEAD:
            if not line[0] == ord('@'):
                if line[0] == ord('>'):
                    raise NotImplementedError('Fasta support not implemented')
                raise RuntimeError('Expected fastq header: {}'.format(line))
            cur['header'] = line
            state = ST_SEQ
        elif state is ST_SEQ:
            cur['seq'] = line
            state = ST_PLUS
        elif state is ST_PLUS:
            if not line == b'+\n':
                raise RuntimeError('Expected "+": {}'.format(line))
            state = ST_SCORE
        elif state is ST_SCORE:
            cur['score'] = line
            if not cur:
                raise RuntimeError('Bad internal state: {}:\n{}'
                                   ''.format(cur, line))
            if len(cur['seq']) != len(cur['score']):
                raise RuntimeError('Sanity check failed: sequence and quality '
                                   'score length differ: {}'.format(cur))

            read = (cur['header'], cur['seq'], cur['score'])
            seq_hash = cur['seq'].__hash__()
            if seq_hash in data:
                data[seq_hash].append(read)
            else:
                data[seq_hash] = [read]

            cur = {}
            state = ST_HEAD
            total_count += 1

        else:
            raise RuntimeError('Illegal state reached: {}'.format(state))

    return data, total_count


def mean_quality_score(read):
    return sum([ord(i) for i in read.score]) / len(read.score)


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
        '-o', '--out-prefix',
        default='derep_',
        help='Prefix attached to output files',
    )
    args = argp.parse_args()

    args.forward_reads.close()
    args.reverse_reads.close()

    fwd_reads = Path(args.forward_reads.name)
    rev_reads = Path(args.reverse_reads.name)

    fwd_data, fwd_total = load(fwd_reads)
    print('{}: total read count: {}, duplicates: {}'
          ''.format(fwd_reads, fwd_total, fwd_total - len(fwd_data)))
    rev_data, rev_total = load(rev_reads)
    print('{}: total read count: {}, duplicates: {}'
          ''.format(rev_reads, rev_total, rev_total - len(rev_data)))

    fwd_out = fwd_reads.parent / (args.out_prefix + fwd_reads.name)
    rev_out = rev_reads.parent / (args.out_prefix + rev_reads.name)

    with fwd_out.open('w') as fout, rev_out.open('w') as rout:

        for fwd_read, rev_read in unique_pairs(fwd_data, rev_data):
            write_read(fwd_read, fout)
            write_read(rev_read, rout)

    #print(fwd, rev)


if __name__ == '__main__':
    main()
