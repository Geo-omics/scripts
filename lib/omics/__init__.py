"""
Package to support geo-omics-scripts
"""
import argparse
import configparser
from importlib import import_module
import inspect
from os import environ
from pathlib import Path
from pkgutil import iter_modules
import re
import subprocess
import sys


OMICS_DIR = '.omics'
CONFIG_FILE = 'config'
CONF_SECTION_PROJECT = 'project'
SCRIPT_PREFIX = 'omics-'
DEFAULT_THREADS = 1
DEFAULT_VERBOSITY = 1


class OmicsArgParser(argparse.ArgumentParser):
    """
    Implements extension to the standard parser
    """
    def __init__(self, *args, project_home=True, threads=True, add_help=True,
                 is_main_omics_parser=False,
                 auto_complete=True, **kwargs):
        """
        Provide the canonical omics argparse argument parser

        :param bool project_home: Include a --project-home option.
        :param bool is_main_omics_parser:  Should be set True for the parser of
                                           the "omics" executable
        :param bool auto_complete:  Weather to allow bash auto completion, this
                                    is set to false by the auto completion code
                                    to avoid infinite recursion
        :param: *args and **kwargs are handed over to argparse.ArgumentParser()

        :return: A new argparse.ArgumentParser object

        Does not use the normal help option since out help option shall go into
        the common group.  --project-dir is made optional as it does not go
        into the init script.  With add_help=False no --help option will be
        parsed, this can be used with partial parsing, and passing --help to a
        subsequent parser, e.g. for the db command which wraps the django
        management commands that have their own arg parsers.
        """
        super().__init__(*args, add_help=False, **kwargs)
        self.auto_complete = 'OMICS_AUTO_COMPLETE' in environ and auto_complete
        self.is_main_omics_parser = is_main_omics_parser
        common = self.add_argument_group('common omics options')

        # help option is inspired by argparse.py
        if add_help:
            common.add_argument(
                '-h', '--help',
                action='help', default=argparse.SUPPRESS,
                help='show this help and exit',
            )
        if project_home:
            common.add_argument(
                '--project-home',
                metavar='PATH',
                help='Omics project directory, by default, this is the '
                     'current directory.',
            )
        if threads:
            common.add_argument(
                '--cpus', '--threads', '-t',
                type=int,
                metavar='N',
                dest='threads',
                default=None,  # None signifies option not given on cmd line
                help='Number of threads / CPUs to employ',
            )
        common.add_argument(
            '-v', '--verbose',
            action='count',
            default=DEFAULT_VERBOSITY,
            dest='verbosity',
            help='Show increased diagnostic output.',
        )
        common.add_argument(
            '--traceback',
            action='store_true',
            help='Show python stack trace in case of some internal errors for '
                 'debugging.',
        )

    def format_help(self):
        """
        Like parent method but abbreviate common options in usage text
        """
        usage = super().format_help()
        # assume -h was defined first and --traceback is last common option
        usage = re.sub(r'\[-h\].*\[--traceback\]', '[OPTIONS...]', usage)
        return usage

    def parse_known_args(self, args=None, namespace=None):
        """
        Parse options and substitue missing options with configured values

         * Missing values for --verbose and --threads are substituted.
         * Remove project_home and add project as object.

        Note: There is some redundancy between the arguments and the project.
        """
        if 'OMICS_AUTO_COMPLETE' in environ:
            try:
                self.bash_complete()
            except:
                if 'OMICS_AUTO_COMPLETE_DEBUG' in environ:
                    raise
            sys.exit()

        args, argv = super().parse_known_args(*args, **kwargs)

        try:
            project = get_project(args.project_home)
        except AttributeError:
            # if project_home is missing, i.e. omics init
            project = None
        else:
            del args.project_home

        args.project = project

        if args.verbosity == DEFAULT_VERBOSITY:
            try:
                args.verbosity = project['verbosity']
            except TypeError:
                pass

        try:
            args.threads
        except AttributeError:
            # constructed with threads=False
            pass
        else:
            if args.threads is None:
                try:
                    args.threads = project['threads']
                except (TypeError, AttributeError, KeyError):
                    args.threads = DEFAULT_THREADS
            else:
                if args.threads <= 0:
                    self.error('The number of threads given via --threads '
                               'must be >=1')

        return args, argv

    def bash_complete(self):
        # get cursor pos
        index = int(environ['OMICS_AUTO_COMPLETE'])

        word = sys.argv[index]

        if 'OMICS_AUTO_COMPLETE_DEBUG' in environ:
            print(file=sys.stderr)
            print('argv:', sys.argv, file=sys.stderr)
            print('cur:', index, '"{}"'.format(word), file=sys.stderr)
        # print('argv:', sys.argv, file=sys.stderr)
        # print('cursor at:', cword, file=sys.stderr)

        # iterate through args left of cursor to get current word type:

        # states:
        #  0 - common optional parameters (cop)
        #  1 - sub command (sub)
        #  2 - sub optional or positionals
        # +nargs
        state = 0
        nargs = 0
        action = None
        for arg in sys.argv[1:index + 1]:
            if 'OMICS_AUTO_COMPLETE_DEBUG' in environ:
                print(state, nargs, '->', arg, file=sys.stderr)

            if not arg:
                # if last word is empty, do not interpret, keep current state
                continue

            if nargs > 0:
                nargs -= 1
                continue
            else:
                # reset action
                if action is not None:
                    action = None

            if state == 0:
                if arg.startswith('-'):
                    for i in self._actions:
                        if arg in i.option_strings:
                            action = i
                            break
                    if action is not None:
                        nargs = action.nargs
                else:
                    # suggest sub-command unless user started a - option
                    state = 1
            elif state == 1:
                # TODO: switch to subcommand parser and completion
                state = 2
            elif state == 2:
                if arg.startswith('-'):
                    for i in self._actions:
                        if arg in i.option_strings:
                            action = i
                            break
                    # assuming arg is valid option
                    if action is not None:
                        nargs = action.nargs

            # fix None nargs, interpret None as 1 arg
            # TODO: handle ? and *
            if nargs is None:
                nargs = 1

        if 'OMICS_AUTO_COMPLETE_DEBUG' in environ:
            print('end state:', state, action, nargs, file=sys.stderr)

        allowed = []
        if nargs > 0:
            # current word is an argument to option
            # TODO: handle files
            pass
        elif state == 0:
            # current word is option string
            for i in self._actions:
                allowed += i.option_strings
        elif state == 1:
            allowed = get_available_commands()

        compl_words = [i for i in allowed if i.startswith(word)]
        if compl_words:
            print(*compl_words, sep='\n')
        sys.exit()


def get_num_cpus():
    """
    Get the number of available CPUs

    :return type: int

    If run in the context of a PBS job, the number of requested CPUs is used.
    On a non-PBS system, i.e. some shared machine, 1/2 the available CPUs are
    taken.
    """
    if 'PBS_ENVIRONMENT' in environ:
        try:
            num_cpus = int(environ.get('PBS_NP'))
        except Exception as e:
            print('omics [WARNING]: Failed to read PBS_NP variable: {}: {}'
                  ''.format(e.__class__.__name__, e), file=sys.stderr)
            num_cpus = 1
    else:
        # `lscpu -p=cpu` prints list of cpu ids, starting with 0, so take last
        # line, divide by 2 (sharing the system) add 1 is the number of CPUs
        # used
        try:
            p = subprocess.run(['lscpu', '-p=cpu'], stdout=subprocess.PIPE)
            p.check_returncode()
            num_cpus = p.stdout.decode().splitlines()[-1].strip()
            num_cpus = int(int(num_cpus) / 2) + 1
        except Exception as e:
            print('omics [WARNING]: Failed to obtain number of CPUs: {}: {}'
                  ''.format(e.__class__.__name__, e), file=sys.stdout)
            num_cpus = 1

    return num_cpus


def get_argparser(*args, **kwargs):
    """
    Provide a canonical omics argparse argument parser

    This function is deprecated and for compatibility only.  Use the
    OmicsArgParser contructor directly.
    """
    return OmicsArgParser(*args, **kwargs)


def process_command_line(command, options, script_dir=Path()):
    """
    Process command line given via the omics command

    :param str command: User-supplied command name
    :param list options: User-supplied options for command
    :param Path script_dir: Directory with scripts

    :return list: Processed command line ready to hand off to subprocess.run
    """
    script = script_dir / (SCRIPT_PREFIX + command)
    return [str(script)] + options


def get_project(path=None):
    """
    Retrieve the current project
    """
    if path is None:
        path = Path.cwd()
    else:
        # allow str input
        path = Path(path)
    return OmicsProject.from_directory(path)


def get_main_arg_parser(*args, **kwargs):
    """
    Get the argument parser for the main omics command
    """
    argp = OmicsArgParser(
        prog=__package__,
        is_main_omics_parser=True,
        *args,
        **kwargs
    )

    argp.add_argument(
        '-n', '--dry-run',
        action='store_true',
        help='Do not actually run command, just print the full command line '
             'that would have been run.',
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
    return argp


def launch_cmd_as_sub_module(args, argp=None):
    """
    Try running the command as a omics.* submodule

    Will return an ImportError if no corresponding module exists.  In all other
    cases this function should not return but handle all exceptions and
    eventually call sys.exit().
    """
    import_err = None
    try:
        cmd_module = import_module('.' + args.command[0],
                                   package=__package__)
    except Exception as e:
        import_err = e
        if args.dry_run or args.verbosity > 1:
            print('{}: {}'.format(e.__class__.__name__, e))
    else:
        try:
            try:
                cmd_module.main(args.command[1:])
            except (AttributeError, TypeError) as e:
                # May happen when there is no function main() or
                # main() does not take a positional arg
                # so we need to distinguish the cases:
                if hasattr(cmd_module, 'main') and \
                  inspect.isfunction(cmd_module.main) and \
                  len(inspect.signature(cmd_module.main).parameters) >= 1:
                    # We assume that call to main() itself succeeded,
                    # so normal error while running the command
                    # to be caught in outer try block
                    raise
                else:
                    # like import error, try calling script later
                    import_err = e
            else:
                # successful run of command
                sys.exit()

        except Exception as e:
            if args.traceback:
                raise
            else:
                argp.error(
                    'Command {} failed: {}: {}'
                    ''.format(args.command[0], e.__class__.__name__, e)
                )
    return import_err


class OmicsProjectNotFound(FileNotFoundError):
    """
    Raised when an operation requires an omics project but none is found
    """
    pass


class OmicsProject(dict):
    """
    Dict-like container to hold an omics project's configuration data

    This is intended to be used like an environment / execution context for
    omics sub-commands.
    """

    default = {
        'project_home': Path.cwd(),  # internal use, not found in conf file
        'name': None,
        'threads': get_num_cpus(),
        'verbosity': DEFAULT_VERBOSITY,
    }
    """ Default settings """

    @classmethod
    def from_directory(cls, path):
        """
        Get the local omics project

        :param Path path: The directory for which the project should be
                          retrieved.
        :return: The omics project object.
        :raise NoOmicsContextFound: If no OMICS_DIR directory with a valid
                                    configuration is found in the given or a
                                    parent directory.
        """
        try:
            path = Path.resolve(path)
        except Exception as e:
            raise OmicsProjectNotFound from e

        omics_dir = None
        for i in [path] + list(path.parents):
            if (i / OMICS_DIR).is_dir():
                omics_dir = i / OMICS_DIR
                break

        if omics_dir:
            config_file = omics_dir / CONFIG_FILE
            if config_file.is_file():
                try:
                    return cls.from_file(config_file)
                except Exception as e:
                    raise OmicsProjectNotFound from e
            else:
                print('Warning: No config file found, using default '
                      'configuration.', file=sys.stderr)

        # use default settings as fallback
        return cls.from_default(project_home=path)

    @classmethod
    def from_default(cls, **kwargs):
        """
        Get a project with minimal default values

        :param dict kwargs: Key-value pairs that override the default.
        """
        p = cls(**cls.default)
        p.update(**kwargs)
        return p

    @classmethod
    def from_file(cls, config_file):
        """
        Get project from config file

        :param Path config_file: Configuration file
        :return: OmicsProject object
        """
        with config_file.open() as f:
            config_str = f.read()

        try:
            return cls._from_str(
                config_str,
                project_home=config_file.parent.parent
            )
        except configparser.MissingSectionHeaderError:
            # add project section
            config_str += '[{}]\n{}'.format(CONF_SECTION_PROJECT, config_str)
            return cls._from_str(config_str)

    @classmethod
    def _from_str(cls, config_str, **kwargs):
        """
        Helper method to parse the contents of a configuration file

        :param str config_str: Configuration file contents
        :param dict kwargs: Overriding key-value pairs

        :return: OmicsConfig object
        """
        proj = cls.from_default(**kwargs)

        parser = configparser.ConfigParser(
            inline_comment_prefixes=('#',),
        )
        parser.read_string(config_str)

        proj.update({
            k: v
            for k, v
            in parser[CONF_SECTION_PROJECT].items()
        })

        proj.fix_types(parser)
        return proj

    def fix_types(self, parser):
        """
        Cast values from str to other types where appropriate

        :param parser: ConfigParser object used to parse the configuration
        """
        get_funs = {
            int: parser.getint,
            float: parser.getfloat,
            bool: parser.getboolean,
        }

        for key in self:
            if self[key] is None:
                # None means variable is unset, no type conversion
                continue

            type_ = type(self.default[key])
            args = CONF_SECTION_PROJECT, key

            if not isinstance(self[key], type_):
                try:
                    self[key] = get_funs[type_](*args)
                except KeyError:
                    pass


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
        # For development environment filter out vim backup files
        return [
            i for i in path.glob(SCRIPT_PREFIX + '*')
            if not i.name.endswith('~')
        ]
    else:
        raise RuntimeError('Failed to determine directory containing omics '
                           'executable: {}'.format(path))


def get_modules_with_main():
    """
    Return list of module with a main function, i.e. the subcommands
    """
    mods = []
    for _, name, ispkg in iter_modules(path=__path__):
        # print(name, ispkg, sep='\n')
        if name == '__main__':
            continue
        mods.append(name)
    return mods


def get_available_commands():
    """
    Get list of available sub-commands

    :return list: List of str names of sub-commands.
                  List is empty in case of errors.
    """
    ret = set()
    try:
        commands = get_available_scripts()
    except:
        pass
    else:
        for i in commands:
            _, _, subcmd = i.name.partition('-')
            if subcmd:
                ret.add(subcmd)
    # do s/_/-/g on submodules since cmd line commands should be slugs
    ret = ret.union([i.replace('_', '-') for i in get_modules_with_main()])
    return sorted(ret)
