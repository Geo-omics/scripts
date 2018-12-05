import sys

from django.conf import settings
from django.core.management import execute_from_command_line

from . import config as db_settings


if __name__ == "__main__":
    settings.configure(**db_settings)
    execute_from_command_line(sys.argv)
