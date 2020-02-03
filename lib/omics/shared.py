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

"""
from array import array
from sys import stderr


DEFAULT_MIN_SAMPLE_SIZE = 0


class MothurShared():
    def __init__(self, file, min_sample_size=DEFAULT_MIN_SAMPLE_SIZE):
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

        self.counts = array('i')
        self.samples = []
        self.sample_sizes = []

        label_ref = None
        for line in file:
            try:
                label, sample, _, *counts = line.strip().split('\t')
            except Exception:
                raise RuntimeError('Failed to parse input: offending line '
                                   'is:\n'.format(line))

            if label_ref is None:
                label_ref = label
                print('label:  ', label, file=stderr)
            elif label_ref != label:
                raise RuntimeError('This shared file contained multiple label,'
                                   ' handling this requires implementation')

            if sample in self.samples:
                raise RuntimeError('got sample {} a second time'
                                   ''.format(sample))

            if len(counts) != self.ncols:
                raise RuntimeError('Not all rows have equal number of columns')

            counts = array('i', map(int, counts))
            size = sum(counts)

            if size < min_sample_size:
                print('Sample {} too small, size {}, ignoring'
                      ''.format(sample, size), file=stderr)
                continue

            self.counts.extend(counts)
            self.samples.append(sample)
            self.sample_sizes.append(size)

        self.nrows = len(self.samples)
        print('total samples:  ', self.nrows, file=stderr)

        # end init

    def rows(self, sample=False):
        for i in range(self.nrows):
            offs = i * self.ncols
            if sample:
                yield self.samples[i], self.counts[offs:offs + self.ncols]
            else:
                yield self.counts[offs:offs + self.ncols]

    def cols(self):
        for i in range(self.ncols):
            yield (self.counts[j * self.ncols + i] for j in range(self.nrows))

    def get_row(self, sample):
        """
        Return list of counts form single sample
        """
        if type(sample) == str:
            offs = self.samples.index(sample) * self.ncols
        else:
            # assume integer index given
            offs = sample * self.ncols

        return self.counts[offs:offs + self.ncols]
