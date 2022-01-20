#! /usr/bin/env python3
import argparse
import json
import platform
import ssl
from urllib.request import urlopen


def update_pypi(system=None):
    '''Check if PyPI needs to be updated'''

    # Get local version tag
    with open('VERSION') as f:
        version = f.read().strip()

    # Get meta from PyPI
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        pass
    else:
        ssl._create_default_https_context = _create_unverified_https_context

    meta = json.load(urlopen('https://pypi.org/pypi/alouette/json'))

    # Look for an existing wheel
    try:
        data = meta['releases'][version]
    except KeyError:
        return True

    if system is None:
        system = platform.system()

    if system == 'Linux':
        tag = 'manylinux'
    elif system == 'Darwin':
        tag = 'macosx'
    else:
        return False

    for d in data:
        if tag in d['filename']:
            return False
    else:
        return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--system', help='target system')
    args = parser.parse_args()

    print(update_pypi(args.system))
