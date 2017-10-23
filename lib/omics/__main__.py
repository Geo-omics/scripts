"""
Invoke omics scripts
"""
import argparse
from pathlib import Path
import subprocess
import sys

from . import SCRIPT_PREFIX, process_command_line


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
        '--script-dir',
        metavar='PATH',
        default=Path(sys.argv[0]).parent,
        help='Path to directory containing the scripts that implement the '
             'omics commands.  This can be used to override the default, '
             'which is: {}'.format(Path(sys.argv[0]).parent),
    )
    argp.add_argument(
        'command',
        nargs=argparse.REMAINDER,
        help='The command to run.',
    )
    args = argp.parse_args()
    args.script_dir = Path(args.script_dir)
    if not args.script_dir.is_dir():
        argp.error('Not a directory: {}'.format(args.script_dir))

    if args.command:
        cmd = args.command[0]
        cmd_opts = args.command[1:]
        cmdline = process_command_line(
            cmd,
            cmd_opts,
            script_dir=args.script_dir
        )
        if args.dry_run:
            print(*cmdline)
        else:
            try:
                p = subprocess.run(cmdline)
            except FileNotFoundError as e:
                if args.traceback:
                    raise
                else:
                    argp.error('Not a valid omics command: {}\n({})'
                               ''.format(cmd, e))
            except Exception as e:
                if args.traceback:
                    raise
                else:
                    argp.error('{}: {}'.format(e.__class__.__name__, e))
            else:
                argp.exit(status=p.returncode)
    else:
        subcmds = get_available_commands()
        if subcmds:
            print('Available commands:')
            for i in subcmds:
                print('  ', i)
            print('Type `omics -h` or `omics <cmd> -h` to get help.')
        else:
            argp.print_help()


def get_available_scripts():
    """
    Try to get the available omics commands

    :return: List of paths for subcommand scripts.
    :raises: In case of errors
    """
    p = subprocess.run(['which', 'omics'], stdout=subprocess.PIPE)
    p.check_returncode()
    path = Path(p.stdout.decode().strip()).parent
    if path.is_dir():
        return list(path.glob(SCRIPT_PREFIX + '*'))
    else:
        raise RuntimeError('Failed to determine directory containing omics '
                           'executable: {}'.format(path))


def get_available_commands():
    """
    Get list of available sub-commands

    :return list: List of str names of sub-commands.
                  List is empty in case of errors.
    """
    ret = []
    try:
        commands = get_available_scripts()
    except:
        pass
    else:
        for i in commands:
            _, _, subcmd = i.name.partition('-')
            if subcmd:
                ret.append(subcmd)
    return ret


if __name__ == '__main__':
    main()
