"""
Omics utilities collection
"""
from collections import defaultdict
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
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(
        [i[0] for i in points],
        [i[1] for i in points],
        marker='.',
    )
    plt.show()
    plt.close()


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
