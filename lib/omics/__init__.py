"""
Package to support geo-omics-scripts
"""
import configparser
from pathlib import Path
import sys


OMICS_DIR = '.omics'
CONFIG_FILE = 'config'
CONF_SECTION_PROJECT = 'project'
SCRIPT_PREFIX = 'omics-'


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


def get_project(cwd=Path.cwd()):
    """
    Retrieve the current project
    """
    return OmicsProject.from_directory(cwd)


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
        'verbosity': 0,
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
