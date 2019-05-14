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
Invoke omics scripts and sub-commands
"""

from pathlib import Path
import os

from . import process_command_line, get_main_arg_parser, get_available_commands
from . import launch_cmd_as_sub_module


def main():
    argp = get_main_arg_parser(description=__doc__)
    args = argp.parse_args()
    args.script_dir = Path(args.script_dir)
    if not args.script_dir.is_dir():
        argp.error('Not a directory: {}'.format(args.script_dir))

    if args.command:
        import_err = launch_cmd_as_sub_module(args, argp)

        # try calling as shellscript
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
                p = os.execv(cmdline[0], cmdline)
            except FileNotFoundError as e:
                if args.traceback:
                    raise
                else:
                    msg1 = '\n  {}: {}'.format(import_err.__class__.__name__,
                                               import_err)
                    msg2 = '  {}: {}'.format(e.__class__.__name__, e)
                    msg3 = '  ==> Not a valid omics command: {}'.format(cmd)
                    argp.error('\n'.join([msg1, msg2, msg3]))
            except Exception as e2:
                if args.traceback:
                    raise
                else:
                    argp.error('Command "{}" failed: {}: {}'
                               ''.format(cmdline, e2.__class__.__name__, e2))
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


if __name__ == '__main__':
    main()
