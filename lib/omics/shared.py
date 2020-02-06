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
from sys import stderr

import numpy
import pandas

DEFAULT_MIN_SAMPLE_SIZE = 0


class MothurShared():
    def __init__(self, file=None, min_sample_size=DEFAULT_MIN_SAMPLE_SIZE):
        file.seek(0)

        head = file.readline().strip()
        if not head.startswith('\t'.join(['label', 'Group', 'numOtus'])):
            raise RuntimeError(
                'input file does not seem to be a mothur shared file, first '
                'line starts with:\n{}'.format(head[:100])
            )

        self.otus = head.split('\t')[3:]
        self.ncols = len(self.otus)
        print('total OTUs:  ', self.ncols, file=stderr)

        counts = numpy.ndarray([], numpy.int32)
        self.samples = []
        self.sample_sizes = []

        self.label = None

        # stack over a generator is deprecated, numpy issues a warning
        # TODO: research alternative
        counts = numpy.stack(
            self._counts_per_sample(file, min_sample_size)
        )

        self.nrows = len(self.samples)
        print('total samples:  ', self.nrows, file=stderr)

        self.counts = pandas.DataFrame(counts, index=self.samples,
                                       columns=self.otus)
        # end init

    def _counts_per_sample(self, file, min_sample_size):
        """
        Generates array of counts for one sample

        This is called by __init__ to import the data.  As side effect it
        collects the sample names and sizes and checks label.
        """
        for line in file:
            try:
                label, sample, _, *counts = line.strip().split('\t')
            except Exception:
                raise RuntimeError('Failed to parse input: offending line '
                                   'is:\n'.format(line))

            if self.label is None:
                self.label = label
                print('label:  ', label, file=stderr)
            elif self.label != label:
                raise RuntimeError('This shared file contained multiple label,'
                                   ' handling this requires implementation')

            if sample in self.samples:
                raise RuntimeError('got sample {} a second time'
                                   ''.format(sample))

            if len(counts) != self.ncols:
                raise RuntimeError('Not all rows have equal number of columns')

            counts = numpy.fromiter(map(int, counts), numpy.int32,
                                    count=self.ncols)
            size = counts.sum()

            if size < min_sample_size:
                print('Sample {} too small, size {}, ignoring'
                      ''.format(sample, size), file=stderr)
                continue

            yield counts

            self.samples.append(sample)
            self.sample_sizes.append(size)

    def rows(self, counts_only=False):
        if counts_only:
            return (list(i[1]) for i in self.counts.iterrows())
        else:
            return ((i, list(row)) for i, row in self.counts.iterrows())

    def cols(self, counts_only=False):
        if counts_only:
            return (list(i[1]) for i in self.counts.iteritems())
        else:
            return ((i, list(col)) for i, col in self.counts.iteritems())

    def get_row(self, sample):
        """
        Return list of counts form single sample
        """
        if type(sample) == str:
            idx = self.samples.index(sample)

        return self.counts[idx, :]

    def pick(self, samples=None, otus=None):
        """
        Keep given samples and OTUs
        """
        if samples is None:
            samples = set(self.samples)
        else:
            samples = set(samples)
            if samples - set(self.samples):
                raise ValueError('invalid sample(s)')

        if otus is None:
            otus = set(self.otus)
        else:
            otus = set(otus)
            if otus - set(self.otus):
                raise ValueError('invalid otu(s)')

        # get indices for array access
        sis = [i for i, j in enumerate(self.samples) if j in samples]
        ois = [i for i, j in enumerate(self.otus) if j in otus]

        # update self
        self.counts = self.counts.take(sis, axis=0).take(ois, axis=1)
        self.sample_sizes = list(self.counts.sum(axis=1))
        self.samples = [i for i in self.samples if i in samples]
        self.otus = [i for i in self.otus if i in otus]
        self.nrows = len(samples)
        self.ncols = len(otus)

    def remove(self, samples=[], otus=[]):
        """
        Remove the given samples and OTUs from the data set
        """
        for i in samples:
            if i not in self.samples:
                raise ValueError('not a sample: {}'.format(i))

        for i in otus:
            if i not in self.otus:
                raise ValueError('not an otu: {}'.format(i))

        self.pick(
            samples=(set(self.samples) - set(samples)),
            otus=(set(self.otus) - set(otus)),
        )

    def save(self, filename):
        """
        Save as mothur shared file under given name or file handle
        """
        with open(filename, 'w') as f:
            row = ['label', 'Group', 'numOtus'] + self.otus
            f.write('\t'.join(row) + '\n')
            for sample, counts in self.rows():
                row = [self.label, sample, str(len(self.otus))]
                row += list(map(str, counts))
                f.write('\t'.join(row) + '\n')
