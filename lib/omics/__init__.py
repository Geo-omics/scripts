"""
Package to support geo-omics-scripts
"""
import argparse
import configparser
from pathlib import Path
import re
import sys


OMICS_DIR = '.omics'
CONFIG_FILE = 'config'
CONF_SECTION_PROJECT = 'project'
SCRIPT_PREFIX = 'omics-'
DEFAULT_VERBOSITY = 1


class OmicsArgParser(argparse.ArgumentParser):
    """
    Implements minor modifications to parent

    Assumes the parser is populated by get_argparser(), relevant changes there
    may need to be reflected here.
    """
    def format_help(self):
        """
        Like parent method but abbreviate common options in usage text
        """
        usage = super().format_help()
        # assume -h was defined first and --traceback is last common option
        usage = re.sub(r'\[-h\].*\[--traceback\]', '[OPTIONS...]', usage)
        return usage


def get_argparser(*args, project_dir=True, **kwargs):
    """
    Provide a canonical omics argparse argument parser

    :param bool project_dir: Include a --project-dir option.
    :param: *args and **kwargs are handed over to argparse.ArgumentParser()

    :return: A new argparse.ArgumentParser object

    Does not use the normal help option since out help option shall go into
    the common group.  --project-dir is made optional as it does not go into
    the init script.
    """
    argp = OmicsArgParser(*args, add_help=False, **kwargs)

    common = argp.add_argument_group('common omics options')

    # help option is inspired by argparse.py
    common.add_argument(
        '-h', '--help',
        action='help', default=argparse.SUPPRESS,
        help='show this help and exit',
    )
    if project_dir:
        common.add_argument(
            '--project-dir',
            metavar='PATH',
            help='Omics project directory, by default, this is the current '
                 'directory.',
        )
    common.add_argument(
        '--verbose', '-v',
        action='count',
        default=1,
        dest='verbosity',
        help='Show increased diagnostic output.',
    )
    common.add_argument(
        '--traceback',
        action='store_true',
        help='Show python stack trace in case of some internal errors for '
             'debugging.',
    )
    return argp


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


class OmicsProjectNotFound(FileNotFoundError):
    """
    Raised when an operation requires an omics project but none is found
    """
    pass


class OmicsProject(dict):
    """
    Dict-like container to hold an omics project's configuration data
    """

    default = {
        'cwd': Path.cwd(),
        'name': None,
        'verbosity': DEFAULT_VERBOSITY,
    }
    """ Default settings """

    @classmethod
    def from_directory(cls, cwd):
        """
        Get the local omics project

        :param Path cwd: The directory for which the project should be
                         retrieved.
        :return: The omics project object.
        :raise NoOmicsContextFound: If no OMICS_DIR directory with a valid
                                    configuration is found in the given or a
                                    parent directory.
        """
        try:
            cwd = Path.resolve(cwd)
        except Exception as e:
            raise OmicsProjectNotFound from e

        omics_dir = None
        for i in [cwd] + list(cwd.parents):
            if (i / OMICS_DIR).is_dir():
                omics_dir = i / OMICS_DIR
                break

        if omics_dir is None:
            raise OmicsProjectNotFound(
                'No omics project found in {} or any parent directory'
                ''.format(cwd)
            )

        config_file = omics_dir / CONFIG_FILE
        if config_file.is_file():
            try:
                return cls.from_file(config_file)
            except Exception as e:
                raise OmicsProjectNotFound from e
        else:
            print('Warning: No config file found, using default configuration.'
                  ' Empty config file created.', file=sys.stderr)
            return cls.from_default()

    @classmethod
    def from_default(cls):
        """
        Return default project
        """
        return cls(**cls.default)

    @classmethod
    def from_file(cls, config_file):
        with config_file.open() as f:
            config_str = f.read()

        try:
            return cls._from_str(config_str)
        except configparser.MissingSectionHeaderError:
            # add project section
            config_str += '[{}]\n{}'.format(CONF_SECTION_PROJECT, config_str)
            return cls._from_str(config_str)

    @classmethod
    def _from_str(cls, config_str):
        proj = cls.from_default()

        parser = configparser.ConfigParser(
            inline_comment_prefixes=('#',),
        )
        parser.read_string(config_str)

        proj.update({
            k: v
            for k, v
            in parser[CONF_SECTION_PROJECT].items()
        })
        return proj
