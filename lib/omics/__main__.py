"""
Invoke omics scripts and sub-commands
"""
from pathlib import Path
import subprocess

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
                p = subprocess.run(cmdline)
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
