import argparse
from io import StringIO
from pathlib import Path
import sys

import django
from django.conf import settings
from django.core.management import call_command, execute_from_command_line


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


argp = argparse.ArgumentParser(add_help=False, prog=PROG_NAME)
argp.add_argument(
    '--db',
    default=DEFAULT_DB_FILE,
    help='Path and name of database file',
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
    db_settings['DATABASES']['default']['NAME'] = db_file_name
    settings.configure(**db_settings)


def main(argv=None):
    """
    Run the 'omics db' command
    """
    args, rest_argv = argp.parse_known_args(argv)
    if not Path(args.db).is_file():
        print('No database present, file not found: {}\nRun "omics db migrate"'
              ' to set up a database in the current directory.'
              ''.format(args.db), file=sys.stderr)

    configure(args.db)
    execute_from_command_line([PROG_NAME] + rest_argv)


if __name__ == "__main__":
    main()
