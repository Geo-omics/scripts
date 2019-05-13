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
Creates a new omics project in the given directory by setting up a
configuration file in a hidden sub-directory ``.omics`` to primarily store
options common to various command in the omics tool chain.  This is purely a
convenience feature and not necessary to run any of the omics tools.
"""

from pathlib import Path
from string import Template

from omics import db
from omics import OMICS_DIR, CONFIG_FILE, CONF_SECTION_PROJECT, get_argparser

empty_conf_template = """\
[$section_name]

"""


def make_minimal_config_file(**options):
    """
    Make up a minimal configuration file

    :param options: Optional key value pairs to store in the configuration
    :return str: Contents of config file
    """
    ret = Template(empty_conf_template).substitute(
        section_name=CONF_SECTION_PROJECT
    )
    for k, v in options.items():
        ret += '{} = {}\n'.format(k, v)
    return ret


def create_config_file(path, **options):
    """
    Write a minimal configuration file to disk

    :param Path path: Full path to the file
    :param options: Optional key value pairs to store in the configuration
    """
    cfg_file = path / CONFIG_FILE
    cfg_file.touch()
    cfg_file.write_text(make_minimal_config_file(**options))


def init(path=Path.cwd(), name=None):
    if isinstance(path, str):
        path = Path(path)

    path = path.resolve()
    omics_dir = path / OMICS_DIR

    if name is None:
        name = path.name

    if not name:
        # maybe path is root
        name = None

    if omics_dir.is_dir():
        # TODO: check config file
        return omics_dir, True
    else:
        omics_dir.mkdir(parents=True)
        create_config_file(omics_dir, name=name)
        return omics_dir, False


def get_argp():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
        project_home=False,
        threads=False,
    )
    argp.add_argument(
        'directory',
        default='.',
        nargs='?',
        help='Path to project directory, by default the current directory.',
    )
    argp.add_argument(
        '--name',
        default=None,
        help='Optional project name, by default, a project name will be '
             'derived from the project directory',
    )
    return argp


def main(argv=None):
    args = get_argp().parse_args(argv)
    path, exists = init(path=args.directory, name=args.name)
    db.setup(db_path=path)
    if exists:
        print('Reinitialized existing omics project in {}'.format(path))
    else:
        print('Initialized omics project in {}'.format(path))


if __name__ == '__main__':
    main()
