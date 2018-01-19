import argparse
from collections import defaultdict, namedtuple
from enum import Enum
from pathlib import Path

from . import get_argparser


def load(path):
    Read = namedtuple('Read', ['header', 'seq', 'score'])
    State = Enum('State', 'head seq plus score')
    state = State.head
    data = {}
    cur = {}
    cur_id = None
    for line in path.open():
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
            if len(cur['seq']) != len(cur['score']):
                raise RuntimeError('Sanity check failed: sequence and quality '
                                   'score length differ: {}'.format(cur))
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


def get_duplicates(data):
    ID_mean = namedtuple('ID_mean', ['id', 'mean_score'])
    seqs = defaultdict(list)
    print('Searching for duplicates...', end='', flush=True)
    for id, read in data.items():
        mean_score = sum([ord(i) for i in read.score]) / len(read.score)
        seqs[read.seq].append(ID_mean(id, mean_score))

    print(' done')
    dups = {}
    for seq, id_means in seqs.items():
        if len(id_means) > 1:
            dups.update(seq=id_means)
            
    return dups


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

    fwd = load(fwd_reads)
    print('Loaded', len(fwd), 'reads from', fwd_reads)
    dupes = get_duplicates(fwd)
    print('Found', len(dupes), '(unique) sequences with multiple occurrences')
    for k,v in list(dupes.items())[:10]:
        print(k, v)

    #rev = load(args.reverse_reads)
    #print(fwd, rev)


if __name__ == '__main__':
    main()
