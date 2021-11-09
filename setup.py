from setuptools import setup


def main():
    setup(name='alouette',
          version='0.0.1',
          description='A TAUOLA wrapper with backward decays',
          author='Valentin Niess',
          packages=['alouette'],
          setup_requires=['cffi>=1.0.0'],
          cffi_modules=['src/build-alouette.py:ffi'],
          install_requires=['cffi>=1.0.0']
    )


if __name__ == '__main__':
    main()
