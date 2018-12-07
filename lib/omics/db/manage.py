import argparse
from pathlib import Path
import sys

from django.conf import settings
from django.core.management import execute_from_command_line


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


def configure(db_file_name=DEFAULT_DB_FILE):
    """
    Do the Django configuration dance

    This needs to be called once whenever we want to access the database
    """
    db_settings['DATABASES']['default']['NAME'] = db_file_name
    settings.configure(**db_settings)
    print('BORK: configured!')


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
