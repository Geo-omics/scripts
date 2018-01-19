import argparse
from collections import OrderedDict, namedtuple
from enum import Enum
from pathlib import Path

from . import get_argparser


def load(path):
    Read = namedtuple('Read', ['header', 'seq', 'score'])
    State = Enum('State', 'head seq plus score')
    state = State.head
    data = OrderedDict()
    cur = {}
    total_count = 0
    for line in path.open():
        line = line.strip()
        if state is State.head:
            if not line[0] == '@':
                raise RuntimeError('Expected fastq header: {}'.format(line))
            cur['header'] = line
            state = State.seq
        elif state is State.seq:
            cur['seq'] = line
            state = State.plus
        elif state is State.plus:
            if not line == '+':
                raise RuntimeError('Expected "+": {}'.format(line))
            state = State.score
        elif state is State.score:
            cur['score'] = line
            if not cur:
                raise RuntimeError('Bad internal state: {}:\n{}'
                                   ''.format(cur, line))
            if len(cur['seq']) != len(cur['score']):
                raise RuntimeError('Sanity check failed: sequence and quality '
                                   'score length differ: {}'.format(cur))

            read = Read(cur['header'], cur['seq'], cur['score'])
            if cur['seq'] in data:
                data[cur['seq']].append(read)
            else:
                data[cur['seq']] = [read]

            cur = {}
            state = State.head
            total_count += 1

        else:
            raise RuntimeError('Illegal state reached: {}'.format(state))

    return data, total_count




def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__
    )
    argp.add_argument('forward_reads', type=argparse.FileType())
    argp.add_argument('reverse_reads', type=argparse.FileType())
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
    dupes = get_duplicates(fwd)
    print('Found', len(dupes), '(unique) sequences with multiple occurrences')

    #print(fwd, rev)


if __name__ == '__main__':
    main()
