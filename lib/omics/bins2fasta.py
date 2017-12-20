"""
Copy contigs into per-bin FASTA files

Takes output from CONCOCT, the result can be fed into CheckM.
"""
import argparse
from pathlib import Path
import sys


def get_argp():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'clustering',
        type=argparse.FileType(),
        help='CONCOCT output file assigning contigs to clusters'
    )
    argp.add_argument(
        'contigs',
        type=argparse.FileType(),
        help='Fasta formatted contigs file'
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
    to_fasta(
        outdir=outdir,
        clustering=args.clustering,
        input_contigs=args.contigs,
        suffix=args.suffix,
        verbose=args.verbose,
    )


def to_fasta(outdir, clustering, input_contigs, suffix='fna', verbose=False):
    outdir.mkdir(parents=True, exist_ok=True)

    # load bin assignments
    contigs = {}
    bins = set()
    print('Reading bin assignment data... ', end='', flush=True)
    for line in clustering:
        line = line.strip()

        try:
            contig_id, xbin = line.split(',')
            xbin = int(xbin)
        except Exception as e:
            print('Failed to parse clustering file {}, offending line: {}: '
                  'Error: {}: {}'
                  ''.format(clustering, line, e.__class__.__name__, e))
            sys.exit(1)

        contigs[contig_id] = xbin
        bins.add(xbin)
    print('done')

    print('  Found {} contigs with bin assignments'.format(len(contigs)))
    print('  Found {} bins'.format(len(bins)))

    print('Copying contigs...', end='', flush=True)
    # get output file handles
    files = {}
    for i in bins:
        suffix = suffix.lstrip('.')
        fname = 'bin_{}.{}'.format(i, suffix)
        files[i] = (outdir / fname).open('w')

    # write contigs
    contig_id = None
    num_not_binned = 0
    for line in input_contigs:
        is_header = line.startswith('>')
        if is_header:
            # account previous contig
            try:
                del contigs[contig_id]
            except KeyError:
                pass
            contig_id, _, _ = line.strip().lstrip('>').partition(' ')

        try:
            bin = contigs[contig_id]
        except KeyError:
            # contig was not binned for whatever reason
            if is_header:
                num_not_binned += 1
            continue
        else:
            files[bin].write(line)

    # account for last contig
    try:
        del contigs[contig_id]
    except KeyError:
        pass

    for i in files.values():
        if i.tell() == 0:
            print('WARNING: {} is empty'.format(i.name))
        i.close()
    print('done')

    if contigs:
        print('WARNING: {} contigs were clustered but not written out (not '
              'found in contigs file?)'.format(len(contigs)))
        if verbose:
            for i, b in contigs.items():
                print(i, b)

    if num_not_binned > 0:
        print('Found {} contigs without bin assignment'.format(num_not_binned))


if __name__ == '__main__':
    main()
