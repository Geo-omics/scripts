# Copyright 2019 Regents of The University of Michigan.

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
Omics utilities collection
"""

from collections import defaultdict
from matplotlib.pyplot import subplots
from pathlib import Path


def load_read_coordinates(reads_file, file_format='fastq'):
    """
    Load read coordinates into dict of lists of points

    :param reads_file: File object or str or Path of inout file
    :param str file_format: Either 'fasta' or 'fastq'
    """
    if file_format == 'fastq':
        mark = '@'
    elif file_format == 'fasta':
        mark = '>'
    else:
        raise ValueError('Illegal file format specified: {}'
                         ''.format(file_format))

    # handle input parameter type
    keep_open = False
    if isinstance(reads_file, str):
        reads_file = open(reads_file)
    elif isinstance(reads_file, Path):
        reads_file = reads_file.open()
    else:
        keep_open = True

    try:
        data = defaultdict(list)
        for line in reads_file:
            if line.startswith(mark):
                try:
                    point = line.split(':')[4:7]
                except Exception:
                    raise RuntimeError('Failed parsing sequence header '
                                       '(split by :): {}'.format(line))
                point[2] = point[2].split()[0]

                try:
                    data[int(point[0])].append((int(point[1]), int(point[2])))
                except Exception:
                    raise RuntimeError('Failed parsing sequence header '
                                       '(int conversion): {}'.format(line))
    except Exception:
        raise
    else:
        return data
    finally:
        if not keep_open:
            reads_file.close()


def scatter(points):
    """
    Diplay scatterplot
    """
    fig, ax = subplots()
    ax.scatter(
        [i[0] for i in points],
        [i[1] for i in points],
        marker='.',
    )
    fig.show()


def hist(data):
    """
    Diplay histogram
    """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(
        data,
        bins='auto',
    )
    plt.show()
    plt.close()
