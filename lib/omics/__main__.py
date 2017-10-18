"""
Invoke omics scripts
"""
import argparse
import subprocess
import sys

from . import process_command_line


def main():
    argp = argparse.ArgumentParser(__package__, description=__doc__)

    argp.add_argument(
        '-n', '--dry-run',
        action='store_true',
        help='Do not actually run command, just print the full command line '
             'that would have been run.',
    )
    argp.add_argument(
        '--traceback',
        action='store_true',
        help='Show python stack trace in case of some internal errors for '
             'debugging.',
    )
    argp.add_argument(
        'command',
        nargs=argparse.REMAINDER,
        help='The command to run.',
    )
    args = argp.parse_args()
    cmd = process_command_line(args.command)
    if args.dry_run:
        print(*cmd)
    else:
        try:
            p = subprocess.run(cmd)
        except Exception as e:
            if args.traceback:
                raise
            else:
                print('{}: {}'.format(e.__class__.__name__, e),
                      file=sys.stderr)
                sys.exit(1)
        else:
            sys.exit(p.returncode)


main()
