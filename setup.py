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


def main():
    with open('VERSION') as f:
        version = f.read().strip()

    with open('README.md') as f:
        long_description = f.read()

    setup(
        name='alouette',
        version=version,
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
        package_data={'alouette': ['lib/libalouette.*', 'include/*.h']},
        entry_points={
            'console_scripts': [
                'alouette-config = alouette.__main__:main',
            ],
        },
    )


if __name__ == '__main__':
    main()
