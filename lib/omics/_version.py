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

import os.path
import subprocess

# Set to real version when distribute outside of git vcs
VERSION = None


def get_version(version=VERSION, raise_on_error=False):
    """
    Get the version string

    Get the hard-coded version if possible, then fall back to ask git.  If that
    fails raise an execption or return an 'unknown' depending on the
    raise_on_error flag.
    """
    if version:
        return version

    try:
        p = subprocess.run(
            ['git', 'describe'],
            cwd=os.path.dirname(__file__),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=raise_on_error,
        )
    except Exception as e:
        out = e.stdout.decode()
        err = e.stderr.decode()
        raise RuntimeError(
            'Failed to get version info from git: {}: {}\n{}{}'
            ''.format(e.__class__.__name__, e, out, err))
    else:
        version = p.stdout.decode().strip()
        # version should be like 1.0.134-42-gd3adb33f
        # make this a PEP440 local version like 1.0.134+42-gd3adb33f
        version = version.replace('-', '+', 1)
        if version:
            return version
    return 'unknown'
