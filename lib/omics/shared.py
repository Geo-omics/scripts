#!/usr/bin/env python3

# Copyright 2020 Regents of The University of Michigan.

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
Implementation of mothur shared file format
"""
from array import array
import argparse
from collections import Counter
from collections import namedtuple
from concurrent.futures import (ProcessPoolExecutor as PoolExecutor,
                                as_completed)
from mmap import mmap, ACCESS_READ
import os.path
from pathlib import Path
import sys

import numpy
import pandas

DEFAULT_THREADS = 1

# file format id names
AUTO_DETECT = 'autodetect'
SHARED = 'mothur shared'
COUNT_TABLE = 'mothur count table'
COMPRESSED_COUNT_TABLE = 'mothur compressed count table'

SPEC_FIELDS = ['index_cols', 'num_index_cols', 'transposed',
               'compressed', 'magic', 'skip_initial_comments']
SPEC = {
    SHARED: {
        'index_cols': ['label', 'Group', 'numOtus'],
        'num_index_cols': 3,
    },
    COUNT_TABLE: {
        'index_cols': ['Representative_Sequence', 'total'],
        'num_index_cols': 2,
        'transposed': True,
    },
    COMPRESSED_COUNT_TABLE: {
        'magic': '#Compressed Format: groupIndex,abundance',
        'skip_initial_comments': True,
        'index_cols': ['Representative_Sequence', 'total'],
        'num_index_cols': 2,
        'transposed': True,
        'compressed': True,
    },
}

# Required fields for the specs are index_cols and num_index_cols
# Defaults for optional fields are:
#  transposed: False
#  compressed: False
#  magic: None  # i.e. derive magic from index column names
#  skip-initial_comments: False
SPEC_DEFAULTS = (False, False, None, False)
Spec = namedtuple('Spec', ['name'] + SPEC_FIELDS, defaults=SPEC_DEFAULTS)


class MothurShared():
    def __init__(self, file_arg, verbose=True, threads=DEFAULT_THREADS,
                 file_format=AUTO_DETECT):
        self.verbose = verbose
        self.threads = threads
        if isinstance(file_arg, str):
            file = open(file_arg, 'r')
        elif isinstance(file_arg, Path):
            file = file_arg.open('r')
        else:
            # assume a file-like object
            file = file_arg

        file.seek(0)
        head = file.readline().strip()

        if file_format == AUTO_DETECT:
            for i, testspec in SPEC.items():
                magic = testspec.get('magic')
                if magic is None:
                    magic = '\t'.join(testspec['index_cols'])
                if head.startswith(magic):
                    file_format = i
                    break
            else:
                # no for loop break
                raise RuntimeError('Fileformat auto-detection failed')

        try:
            self.spec = Spec(name=file_format, **SPEC[file_format])
        except KeyError:
            raise ValueError('Invalid file format identifier')

        if self.spec.skip_initial_comments:
            while head.startswith('#'):
                head = file.readline().strip()

        if not head.startswith('\t'.join(self.spec.index_cols)):
            raise RuntimeError(
                'input file does not seem to be a {} file, table column header'
                'starts with:\n{}'.format(self.spec.name, head[:100])
            )

        input_cols = head.split('\t')[self.spec.num_index_cols:]

        if self.spec.transposed:
            self.samples = input_cols
            self.info('total samples:  ', len(self.samples))
        else:
            self.otus = input_cols
            self.info('total OTUs:  ', len(self.otus))

        self.label = None

        row_ids, counts = self.get_counts(file)

        if len(row_ids) != len(set(row_ids)):
            raise RuntimeError('Duplicate row ids in intput file')

        if self.spec.transposed:
            self.otus = row_ids
            counts = counts.T
        else:
            self.samples = row_ids

        self.info('Making data frame...', end='', flush=True)
        self.counts = pandas.DataFrame(counts, index=self.samples,
                                       columns=self.otus)
        self.info('[ok]', inline=True)
        self._update_from_counts(trim=False)
        if self.spec.transposed:
            self.info('total otus:  ', self.ncols)
        else:
            self.info('total samples:  ', self.nrows)
        self.info('label:  ', self.label)

        if isinstance(file_arg, (str, Path)):
            file.close()
        # end init

    def _counts_per_sample(self, file):
        """
        Generates array of counts for one sample

        This is called by __init__ to import the data.  As side effect it
        collects the sample names and sizes and checks label.
        """
        self._input_row_ids = []
        if 'label' in self.spec.index_cols:
            # store sample ids as list here to make to index later
            self.label = None

        for line in file:

            row_id, row_meta, counts = self.process_line(line)
            self._input_row_ids.append(row_id)

            if 'label' in self.spec.index_cols:
                # row_meta is the label
                if row_meta != self.label:
                    if self.label is None:
                        self.label = row_meta
                    else:
                        raise RuntimeError(
                            'This shared file contained multiple label, '
                            'handling this requires implementation: {} != {}'
                            ''.format(row_meta, self.label))

            yield numpy.array(counts)

    def get_counts(self, file):
        """
        Make a numpy array, single-thread implementation
        """
        if self.threads > 1:
            return self.get_counts_mp(file)

        counts = numpy.ndarray([], numpy.int32)

        self._input_row_ids = []
        # stack over a generator is deprecated, numpy issues a warning
        # TODO: research alternative
        counts = numpy.stack(
            self._counts_per_sample(file)
        )
        return self._input_row_ids, counts

    def get_counts_mp(self, file):
        """
        Make numpy array, multi-processing+futures implementation

        Using the pre-3.8 functionality, (it seems) it's not possible to use
        shared memory to send the data back from the children.  So we use
        Pool/starmap and receive arrays through pickling, assemble them and
        wrap a ndarray around.

        TODO: for 3.8 explore the shared_memory submodule
        """
        fsize = os.path.getsize(file.name)
        self.info('filesize:', fsize)
        breaks = numpy.linspace(0, fsize, num=self.threads + 1, dtype=int)
        sargs = []
        for i in range(len(breaks) - 1):
            sargs.append((file.name, breaks[i], breaks[i + 1]))

        row_ids = []
        counts = array('i')
        with PoolExecutor(max_workers=self.threads) as pe:
            futmap = {pe.submit(self.chunk2array, *i): i[1] for i in sargs}
            for future in as_completed(futmap):
                start = futmap[future]
                try:
                    chunk_row_ids, chunk_meta, a = future.result()
                except Exception as e:
                    raise RuntimeError(
                        'Chunk starting at {} failed: {}: {}'
                        ''.format(start, e.__class__.__name__, e)
                    )
                else:
                    if 'label' in self.spec.index_cols:
                        # FIXME: last worker wins, needs a check
                        self.label = chunk_meta
                    row_ids += chunk_row_ids
                    counts.extend(a)

        counts = numpy.frombuffer(counts, dtype=numpy.dtype('int32'))
        counts = counts.reshape(len(row_ids), len(counts) // len(row_ids))
        return row_ids, counts

    def chunk2array(self, path, start, end):
        """
        Process one chunk of input file into an array

        This should be run in multi-processing.  If 'start' falls between two
        lines, the partial line is ignored and processing starts with the next
        full line.  Likewise, a line extending beyond 'end' will be fully
        processes.  If no newline falls between 'start' and 'end' then nothing
        is done. The usual range-semantic is used for 'start' and 'end'.

        Using the mmap logic seems to be much faster then regular readline over
        the file object.
        """
        a = array('i')
        chunk_meta = None
        if 'label' in self.spec.index_cols:
            label = None
        row_ids = []
        with open(path, 'rb') as f:
            mm = mmap(f.fileno(), 0, access=ACCESS_READ)
            mm.seek(start)
            if start > 0 and mm[start - 1] != ord('\n'):
                # consume partial line
                mm.readline()
            if start == 0:
                # skip header
                mm.readline()

            while mm.tell() < end:
                row_id, row_meta, counts = self.process_line(
                    mm.readline().decode()
                )
                a.extend(counts)
                row_ids.append(row_id)

                if 'label' in self.spec.index_cols:
                    # row_meta is the label
                    # FIXME: setting the label is not coordinated among workers
                    #        last worker process wins??
                    if row_meta != label:
                        if label is None:
                            label = row_meta
                        else:
                            raise RuntimeError(
                                'This shared file contained multiple label, '
                                'handling this requires implementation:'
                                '{} != {}'.format(label, row_meta))

        if 'label' in self.spec.index_cols:
            chunk_meta = label

        if len(a) > 0.99 * (2 ** 29):
            # max space for pickling is 2GB, 512x10^6 32bit ints
            # see https://github.com/python/cpython/pull/10305
            # fix is in 3.8
            self.info('[WARNING] too much data in child process for pickling')
            self.info('array size is', len(a))
        return row_ids, chunk_meta, a

    def process_line(self, line):
        """
        Get count array from line
        """
        row_meta = None
        try:
            if self.spec.name == SHARED:
                # first field is label, return as row meta data
                # second field is sample, return as row id
                # third field is numOtus? -- ignore
                row_meta, row_id, _, *counts = line.strip().split('\t')
            elif self.spec.compressed:
                raise NotImplementedError
            else:
                # count tables:
                # first field is sequence id, return as row id
                # second field is the total for the sequence, ignore
                row_id, _, *counts = line.strip().split('\t')
        except Exception:
            raise RuntimeError('Failed to parse input: offending line '
                               'is:\n'.format(line))

        return row_id, row_meta, array('i', map(int, counts))

    def rows(self, counts_only=False, as_iter=False):
        it = self.counts.iterrows()
        if counts_only:
            it = (i[1] for i in it)
            if as_iter:
                return it
            else:
                return list(map(list, it))
        else:
            if as_iter:
                return it
            else:
                return [(i, list(row)) for i, row in it]

    def cols(self, counts_only=False, as_iter=False):
        it = self.counts.iteritems()
        if counts_only:
            it = (i[1] for i in it)
            if as_iter:
                return it
            else:
                return list(map(list, it))
        else:
            if as_iter:
                return it
            else:
                return [(i, list(col)) for i, col in it]

    def get_row(self, sample):
        """
        Return list of counts form single sample
        """
        return self.counts.loc[sample]

    def pick_otus(self, otus):
        """
        Pick the given OTUs from the data set
        """
        self.counts = self.counts.loc[:, otus]
        self._update_from_counts(trim=False)

    def pick_samples(self, samples, trim=True):
        """
        Pick the given samples from the data set
        """
        self.counts = self.counts.loc[samples]
        self._update_from_counts(trim=trim)

    def _update_from_counts(self, trim=True):
        """
        updates self after changes to counts

        :param bool trim: Remove zero OTUs
        """
        if trim:
            old_ncols = len(self.counts.columns)
            self.counts = self.counts.loc[:, self.counts.sum() > 0]
            new_ncols = len(self.counts.columns)
            if old_ncols != new_ncols:
                self.info('Removed OTUs with all-zero count:',
                          old_ncols - new_ncols)

        self.sample_sizes = self.counts.sum(axis=1)
        self.samples = self.counts.index
        self.otus = self.counts.columns
        self.nrows = len(self.samples)
        self.ncols = len(self.otus)

    def remove_otus(self, otus):
        """
        Remove the given OTUs from the data set
        """
        self.counts.drop(otus, axis=1, inplace=True)
        self._update_from_counts(trim=False)

    def remove_samples(self, samples, trim=True):
        """
        Remove the given samples from the data set
        """
        self.counts.drop(samples, axis=0, inplace=True)
        self._update_from_counts(trim=trim)

    def save(self, filename):
        """
        Save as mothur shared file under given name or file handle
        """
        with open(filename, 'w') as f:
            row = ['label', 'Group', 'numOtus'] + list(self.otus)
            num_otus = str(len(self.otus))
            f.write('\t'.join(row) + '\n')
            for sample, counts in self.rows():
                row = [self.label, sample, num_otus]
                row += list(map(str, counts))
                f.write('\t'.join(row) + '\n')

    def info(self, *args, inline=False, **kwargs):
        if self.verbose:
            if inline:
                print(*args, file=sys.stderr, **kwargs)
            else:
                print('[{}]'.format(Path(sys.argv[0]).name), *args,
                      file=sys.stderr, **kwargs)


class Groups():
    """
    Implements support for mothur groups files
    """
    def __init__(self, file):
        if isinstance(file, str):
            file = open(file)
        elif isinstance(file, Path):
            file = file.open()
        else:
            # assume it's file-like
            pass

        self.counts = pandas.Series(Counter(
            (i.strip().split('\t')[1] for i in file)
        ))


def __main__():
    # run test
    fmt_args = {
        'auto': AUTO_DETECT,
        'shared': SHARED,
        'count-table': COUNT_TABLE,
        'compressed': COMPRESSED_COUNT_TABLE,
    }
    argp = argparse.ArgumentParser(description='Test importing a shared file '
                                   'via command line')
    argp.add_argument('shared_file', help='Name of shared file')
    argp.add_argument(
        '-f', '--format',
        choices=fmt_args.keys(),
        default='shared',
        help='Input file format',
    )
    argp.add_argument('-t', type=int, default=DEFAULT_THREADS,
                      help='number of threads')
    args = argp.parse_args()
    try:
        file_format = fmt_args[args.format]
    except KeyError:
        argp.error('Invalid value for --format option')
    MothurShared(args.shared_file, threads=args.t, file_format=file_format)


if __name__ == '__main__':
    __main__()
