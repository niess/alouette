import os
import subprocess
from setuptools import setup


CLASSIFIERS = '''\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)
Programming Language :: Python
Topic :: Scientific/Engineering :: Physics
Operating System :: POSIX :: Linux
'''


def get_version():
    '''Get the version tag'''
    version = os.getenv('ALOUETTE_VERSION')
    if not version:
        version = '0.0.1'

    p = subprocess.Popen(
        'git describe --match=NeVeRmAtCh --always --dirty',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    stdout, _ = p.communicate()
    try:
        stdout = stdout.decode()
    except AttributeError:
        stdout = str(stdout)
    git_revision = stdout.strip()

    with open('alouette/version.py', 'w+') as f:
        f.write(
            f'''\
# This file was generated by setup.py
VERSION = '{version:}'
GIT_REVISION = '{git_revision:}'
'''
        )

    return version


def main():

    with open('README.md') as f:
        long_description = f.read()

    setup(
        name='alouette',
        version=get_version(),
        author='Valentin Niess',
        author_email='valentin.niess@gmail.com',
        description='A TAUOLA wrapper with backward decays',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/niess/alouette',
        download_url='https://pypi.python.org/pypi/alouette',
        project_urls={
            'Bug Tracker': 'https://github.com/niess/alouette/issues',
            'Source Code': 'https://github.com/niess/alouette',
        },
        packages=['alouette'],
        classifiers=[s for s in CLASSIFIERS.split(os.linesep) if s.strip()],
        license='GPLv3',
        platforms=['Linux'],
        python_requires='>=3.6',
        setup_requires=['cffi>=1.0.0', 'pcpp>=1.22'],
        cffi_modules=['src/build-alouette.py:ffi'],
        install_requires=['cffi>=1.0.0', 'numpy'],
        package_data={'alouette': ['libalouette.*']},
    )


if __name__ == '__main__':
    main()
