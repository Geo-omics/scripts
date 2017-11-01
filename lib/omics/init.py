"""
Initialize a new omics project directory
"""
from pathlib import Path
from string import Template

from omics import OMICS_DIR, CONFIG_FILE, CONF_SECTION_PROJECT, get_argparser

empty_conf_template = """\
[$section_name]

"""


def make_minimal_config_file(**options):
    """
    Make up a minimal configuration file

    :param options: Optional key value pairs to store in the configuration
    :return str: Contents of config file
    """
    ret = Template(empty_conf_template).substitute(
        section_name=CONF_SECTION_PROJECT
    )
    for k, v in options.items():
        ret += '{} = {}\n'.format(k, v)
    return ret


def create_config_file(path, **options):
    """
    Write a minimal configuration file to disk

    :param Path path: Full path to the file
    :param options: Optional key value pairs to store in the configuration
    """
    cfg_file = path / CONFIG_FILE
    cfg_file.touch()
    cfg_file.write_text(make_minimal_config_file(**options))


def init(path=Path.cwd(), name=None):
    if isinstance(path, str):
        path = Path(path)

    path = path.resolve()
    omics_dir = path / OMICS_DIR

    if name is None:
        name = path.name

    if not name:
        # maybe path is root
        name = None

    if omics_dir.is_dir():
        print('Reinitialized existing omics project in {}'.format(omics_dir))
        # TODO: check config file
    else:
        omics_dir.mkdir(parents=True)
        create_config_file(omics_dir, name=name)
        print('Initialized omics project in {}'.format(omics_dir))


def main():
    argp = get_argparser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
        project_home=False,
        threads=False,
    )
    argp.add_argument(
        'directory',
        default='.',
        nargs='?',
        help='Path to project directory, by default the current directory.',
    )
    argp.add_argument(
        '--name',
        default=None,
        help='Optional project name, by default, a project name will be '
             'derived from the project directory',
    )
    args = argp.parse_args()
    init(path=args.directory, name=args.name)


if __name__ == '__main__':
    main()
