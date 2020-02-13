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
        self._update_from_counts()

    def pick_samples(self, samples):
        """
        Pick the given samples from the data set
        """
        self.counts = self.counts.loc[samples]
        self._update_from_counts()

    def _update_from_counts(self):
        """
        updates self after changes to counts
        """
        self.sample_sizes = list(self.counts.sum(axis=1))
        self.samples = list(self.counts.index)
        self.otus = list(self.counts.columns)
        self.nrows = len(self.samples)
        self.ncols = len(self.otus)

    def remove_otus(self, otus):
        """
        Remove the given OTUs from the data set
        """
        self.counts.drop(otus, axis=1, inplace=True)
        self._update_from_counts()

    def remove_samples(self, samples):
        """
        Remove the given samples from the data set
        """
        self.counts.drop(samples, axis=0, inplace=True)
        self._update_from_counts()

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
