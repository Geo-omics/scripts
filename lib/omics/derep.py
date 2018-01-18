import argparse
from collections import namedtuple
from enum import Enum

from . import get_argparser


def load(file):
    Read = namedtuple('Read', ['header', 'seq', 'score'])
    State = Enum('State', 'head seq plus score')
    state = State.head
    data = {}
    cur = {}
    cur_id = None
    for line in file:
        line = line.strip()
        if state is State.head:
            if not line[0] == '@':
                raise RuntimeError('Expected fastq header: {}'.format(line))
            cur['header'] = line
            try:
                cur_id = tuple((
                    int(i) for i in
                    line[1:].split()[0].split(':')[4:7]
                ))
            except Exception:
                print('Failed to parse header:', line)
                raise

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
            if cur_id is None or not cur:
                raise RuntimeError('Bad internal state: {} {}\n{}'
                                   ''.format(cur_id, cur, line))
            try:
                data[cur_id] = Read(**cur)
            except Exception:
                print('Failed to constract Read object:', cur_id, '\n', cur)
                raise
            else:
                cur = {}
                cur_id = None
                state = State.head

        else:
            raise RuntimeError('Illegal state reached: {}'.format(state))

    return data


def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__
    )
    argp.add_argument('forward_reads', type=argparse.FileType())
    argp.add_argument('reverse_reads', type=argparse.FileType())
    args = argp.parse_args()

    fwd = load(args.forward_reads)
    rev = load(args.reverse_reads)
    print(fwd, rev)


if __name__ == '__main__':
    main()
