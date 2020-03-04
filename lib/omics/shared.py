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
from concurrent.futures import (ProcessPoolExecutor as PoolExecutor,
                                as_completed)
from mmap import mmap, ACCESS_READ
import os.path
from pathlib import Path
import sys

import numpy
import pandas

DEFAULT_THREADS = 1


class MothurShared():
    def __init__(self, file_arg=None, verbose=True, threads=DEFAULT_THREADS):
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
        if not head.startswith('\t'.join(['label', 'Group', 'numOtus'])):
            raise RuntimeError(
                'input file does not seem to be a mothur shared file, first '
                'line starts with:\n{}'.format(head[:100])
            )

        otus = head.split('\t')[3:]
        self.ncols = len(otus)
        self.info('total OTUs:  ', self.ncols)

        self.label = None
        self.samples = []

        counts = self.get_counts(file)

        self.info('Making data frame...', end='', flush=True)
        self.counts = pandas.DataFrame(counts, index=self.samples,
                                       columns=otus)
        self.info('[ok]', inline=True)
        self._update_from_counts(trim=False)
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
        # store sample ids as list here to make to index later
        self.samples = []
        self.label = None
        for line in file:
            label, sample, counts = self.process_line(line)

            yield numpy.array(counts)

            self.samples.append(sample)
            if label != self.label:
                if self.label is None:
                    self.label = label
                else:
                    raise RuntimeError(
                        'This shared file contained multiple label, handling '
                        'this requires implementation: {} != {}'
                        ''.format(label, self.label))

    def get_counts(self, file):
        """
        Make a numpy array, single-thread implementation
        """
        if self.threads > 1:
            return self.get_counts_mp(file)

        counts = numpy.ndarray([], numpy.int32)

        # stack over a generator is deprecated, numpy issues a warning
        # TODO: research alternative
        counts = numpy.stack(
            self._counts_per_sample(file)
        )
        return counts

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

        counts = array('i')
        with PoolExecutor(max_workers=self.threads) as pe:
            futmap = {pe.submit(self.chunk2array, *i): i[1] for i in sargs}
            for future in as_completed(futmap):
                start = futmap[future]
                try:
                    label, s, a = future.result()
                except Exception as e:
                    raise RuntimeError(
                        'Chunk starting at {} failed: {}: {}'
                        ''.format(start, e.__class__.__name__, e)
                    )
                else:
                    self.label = label
                    self.samples += s
                    counts.extend(a)

        counts = numpy.frombuffer(counts, dtype=numpy.dtype('int32'))
        counts = counts.reshape(len(self.samples), self.ncols)
        return counts

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
        label = None
        samples = []
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
                label0, sample, counts = self.process_line(
                    mm.readline().decode()
                )
                a.extend(counts)
                samples.append(sample)

                if label0 != label:
                    if label is None:
                        label = label0
                    else:
                        raise RuntimeError(
                            'This shared file contained multiple label, '
                            'handling this requires implementation: {} != {}'
                            ''.format(label, label0))

        if len(a) > 0.99 * (2 ** 29):
            # max space for pickling is 2GB, 512x10^6 32bit ints
            # see https://github.com/python/cpython/pull/10305
            # fix is in 3.8
            self.info('[WARNING] too much data in child process for pickling')
            self.info('array size is', len(a))
        return label, samples, a

    def process_line(self, line):
        """
        Get count array from line
        """
        try:
            label, sample, _, *counts = line.strip().split('\t')
        except Exception:
            raise RuntimeError('Failed to parse input: offending line '
                               'is:\n'.format(line))

        if sample in self.samples:
            raise RuntimeError('got sample {} a second time'
                               ''.format(sample))

        if len(counts) != self.ncols:
            raise RuntimeError('Not all rows have equal number of columns')

        return label, sample, array('i', map(int, counts))

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


def __main__():
    # run test
    argp = argparse.ArgumentParser(description='Test importing a shared file '
                                   'via command line')
    argp.add_argument('shared_file', help='Name of shared file')
    argp.add_argument('-t', type=int, help='number of threads')
    args = argp.parse_args()
    MothurShared(args.shared_file, threads=args.threads)


if __name__ == '__main__':
    __main__()
