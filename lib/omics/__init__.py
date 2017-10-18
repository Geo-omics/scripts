"""
Package to support geo-omics-scripts
"""
from pathlib import Path


SCRIPT_PREFIX = 'omics-'


def process_command_line(cmdline):
    """
    Process command line given via the omics command

    :param list cmdline: User-supplied command line
    :return list: Processed command line ready to hand off to subprocess.run
    """
    dev_env_root = Path(__file__).parent.parent.parent
    if dev_env_root.name == 'geo-omics-scripts':
        # assume dev environment
        script = dev_env_root.resolve() / 'scripts' / SCRIPT_PREFIX
        cmdline = [str(script) + cmdline[0]] + cmdline[1:]
    return cmdline
