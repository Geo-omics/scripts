"""
Module dealing with chopped bins.
"""
import argparse
from itertools import groupby
from pathlib import Path
import sys


# Data structures:
#
# contig indexed:
# dict contig_info = contig_id --> dict (length, chunks = dict chunk->bin)
#
# bin indext:
# dict bin_info = bin_id --> (file, dict contig)
def load_contig_ids(file, no_chunk=False):
    """
    Generate chop/chunk-aware contig ids from a fasta file

    :param bool no_chunk:  If true, then skip parsing for chunk info and
                           return a simple list of contig ids

    Yields tuples of contig id and chunk number.
    """
    with file.open() as f:
        for line in f:
            if line.startswith('>'):
                if no_chunk:
                    yield line.strip()
                else:
                    # parse >k255_99.7 -> ('k255_99', 7)
                    cont_id = line.strip().lstrip('>').rsplit('.', 1)
                    if len(cont_id) == 1:
                        # no chunk info, keep as-is
                        yield cont_id[0], None
                    elif len(cont_id) == 2:
                        try:
                            yield cont_id[0], int(cont_id[1])
                        except ValueError:
                            raise RuntimeError('Failed to parse int chunk id '
                                               'in {}: {}'.format(file, line))
                    else:
                        raise RuntimeError('Failed to parse contig id in '
                                           '{}: {}'.format(file, line))


def read_chopping(assembly):
    """
    Compile chopping info from (chopped-up) assembly fasta file

    :param pathlib.Path assembly:  Path to input file

    Assumes that the chop-contig tool of the goe-omics-scripts
    package was used for the copping and the resulting chunk names.

    Returns a dict mapping contig id to chunk count
    """
    data = {}
    contig_info = sorted(load_contig_ids(assembly))
    for contig_id, items in groupby(contig_info, lambda x: x[0]):
        chunks = [i[1] for i in items]  # strip out contig_id
        if chunks == [None]:
            count = 1
        else:
            # sanity check on chunks
            if chunks != list(range(len(chunks))):
                raise RuntimeError('Chunk check failed: {} / {}'
                                   ''.format(contig_id, chunks))
            count = len(chunks)

        data[contig_id] = {'length': count}
    return data


def read_bins_metablast(bin_path, contig_info):
    """
    Compile binning info from Metablast run

    :param pathlib.Path bin_path: Path were matablast stored chopped bins
    :param dict contig_info:  Contig-indexed data
    """
    # compile binning from Metablast output, which is a directory full of fasta
    # files, one per bin, into a data structure
    # bins: dict int -> (file, bin_contig_list)
    bins = {}
    for i in bin_path.glob('bin.*.fa'):
        try:
            bin_num = int(i.name.split('.')[1])  # expect name to be bin.NNN.fa
        except Exception as e:
            raise RuntimeError('Failed to parse bin file name: {} '.format(i))

        bin_data = {}  # to collect info on contigs

        # like read_chopping() but we need more details
        binning = sorted(load_contig_ids(i))
        for contig_id, items in groupby(binning, lambda x: x[0]):
            chunks = [i[1] for i in items]  # strip out contig_id

            # compute missing chunk info
            try:
                total = contig_info[contig_id]['length']
            except Exception as e:
                print(contig_info[contig_id])
                raise RuntimeError(
                    'BORKED: {}: {} / contig_id: {} / file: {}'
                    ''.format(e.__class__.__name__, e, contig_id, i)
                )

            if total == len(chunks):
                missing = []
            else:
                should = list(range(total))
                missing = list(sorted(set(should) - set(chunks)))

            bin_data[contig_id] = {
                'chunks': chunks,
                'total': total,
                'missing': missing,
            }

            if contig_id not in contig_info:
                raise RuntimeError('Missing data: Expected contig {} from {} '
                                   'in contig info'.format(contig_id, i))

            # update contig index
            if 'binning' not in contig_info[contig_id]:
                contig_info[contig_id]['binning'] = [None] * total
            for j in chunks:
                # assumes not double binning
                contig_info[contig_id]['binning'][j] = bin_num

        bins[bin_num] = (i, bin_data)
    return bins


def print_split_contigs(contig_info, file=sys.stdout):
    """
    Print contigs with split binning
    """
    for contig_id, info in contig_info.items():
        if 'binning' in info:
            bins = {k: [] for k in set(info['binning'])}
            if len(bins) > 1:
                for i in bins:
                    bins[i] = [
                        j for j
                        in range(info['length'])
                        if info['binning'][j] == i
                    ]
                print(contig_id, bins, sep='\t', file=file)


def main(asm, chopped_asm, bin_path):
    print('Reading chop info... ', end='', flush=True)
    contig_info = read_chopping(chopped_asm)
    print('[done]')

    print('Reading binning info... ', end='', flush=True)
    bin_info = read_bins_metablast(bin_path, contig_info)
    print('[done]')

    return contig_info, bin_info


if __name__ == '__main__':
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        '-a', '--assembly',
        type=argparse.FileType(),
        help='The (presumely) chopped assembly in fasta format',
    )
    argp.add_argument(
        '-b', '--bin-path',
        help='Directory were bins are stored',
    )
    args = argp.parse_args()
    args.bin_path = Path(args.bin_path)
    args.assembly.close()
    args.assembly = Path(args.assembly.name)

    c, _ = main(None, args.assembly, args.bin_path)
    print_split_contigs(c)
