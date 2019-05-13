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

from io import StringIO
from pathlib import Path
import sys

import django
from django.conf import settings
from django.core.management import call_command, execute_from_command_line

from omics import OmicsArgParser, OMICS_DIR


DEFAULT_DB_FILE = 'omics.db'
PROG_NAME = 'omics db'  # for argparse and django management cmd benefit

db_settings = {
    'INSTALLED_APPS': ['omics.db.apps.OmicsDBConfig'],
    'DATABASES': {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': None,
        },
    },
}


argp = OmicsArgParser(prog=PROG_NAME, add_help=False)
argp.add_argument(
    '--db',
    default=None,
    help='Path and name of database file. If this is not provided, then the '
         'database is looked up in the {} project directory, if one exists, '
         'and otherwise the default {} in the current directory is used.'
         ''.format(OMICS_DIR, DEFAULT_DB_FILE),
)


def setup(db_path):
    """
    Setup and configure django framework

    This is public API exposed at the package level to be used by other code
    outside of omics.db but do not call this from within omics.db since
    django.setup() is called here only for the benefit of external scripts and
    per Django documentation should not be call explicitly in code using the
    framework, i.e. everything under omics.db.
    """
    if db_path.is_dir():
        db_file_name = str(db_path / DEFAULT_DB_FILE)
    else:
        db_file_name = str(db_path)

    configure(db_file_name)
    django.setup()
    out = StringIO()
    try:
        call_command('migrate', stdout=out)
    except Exception as e:
        print('Database setup / migration failed: {}: {}, output of the '
              'migrate command was:'.format(e.__class__.__name__, e),
              file=sys.stderr)
        print(out.getvalue(), file=sys.stderr)


def configure(db_file_name=DEFAULT_DB_FILE):
    """
    Do the Django settings configuration dance

    This needs to be called once whenever we want to access the database
    """
    db_settings['DATABASES']['default']['NAME'] = str(db_file_name)
    settings.configure(**db_settings)


def main(argv=None):
    """
    Run the 'omics db' command
    """
    args, rest_argv = argp.parse_known_args(argv)

    if args.db is None:
        if hasattr(args, 'project'):
            if (args.project['project_home'] / OMICS_DIR).is_dir():
                # omics init was called
                db_path = \
                    args.project['project_home'] / OMICS_DIR / DEFAULT_DB_FILE
            else:
                db_path = args.project['project_home'] / DEFAULT_DB_FILE
        else:
            # FIXME: when will this path be taken?
            db_path = Path.cwd() / DEFAULT_DB_FILE
    else:
        db_path = Path(args.db)

    if not Path(db_path).is_file():
        print('No database present, file not found: {}\nTo properly initiate '
              'a project run "omics init" or run "omics db migrate" '
              'to just set up a database in the current directory.'
              ''.format(args.db), file=sys.stderr)

    configure(db_path)
    execute_from_command_line([PROG_NAME] + rest_argv)


if __name__ == "__main__":
    main()
