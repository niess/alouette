#! /usr/bin/env python3
import json
import ssl
from urllib.request import urlopen


def get_local_version():
    '''Get local version tag'''
    with open('VERSION') as f:
        return f.read().strip()


def get_pypi_version():
    '''Get PypI (remote) version tag'''
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        pass
    else:
        ssl._create_default_https_context = _create_unverified_https_context

    meta = json.load(urlopen('https://pypi.org/pypi/alouette/json'))
    return meta['info']['version']


def update_pypi():
    '''Check if PyPI needs to be updated'''
    loc = get_local_version()
    pypi = get_pypi_version()
    print(loc != pypi)


if __name__ == '__main__':
    update_pypi()
