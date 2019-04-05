import os.path
import subprocess

VERSION = None


def get_version(version=VERSION, raise_on_error=False):
    """
    Get the version string

    Get the hard-coded version if possible, then fall back to ask git.  If that
    fails raise an execption or return an 'unknown' depending on the
    raise_on_error flag.
    """
    if version:
        return version

    try:
        p = subprocess.run(
            ['git', 'describe'],
            cwd=os.path.dirname(__file__),
            stdout=subprocess.PIPE,
        )
    except Exception as e:
        if raise_on_error:
            raise RuntimeError('Failed to get version info from git: {}: {}'
                               ''.format(e.__class__.__name__, e))
    else:
        version = p.stdout.decode().strip()
        # version should be like 1.0.134-42-gd3adb33f
        # make this a PEP440 local version like 1.0.134+42-gd3adb33f
        version = version.replace('-', '+', 1)
        if version:
            return version
    return 'unknown'
