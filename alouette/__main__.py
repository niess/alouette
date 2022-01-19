#! /usr/bin/env python3
import argparse
from pathlib import Path
import platform
import sys

from . import version


__all__ = ('main', 'PREFIX')


PREFIX = Path(__file__).parent.resolve()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', help='Alouette installation prefix',
                        action='store_true')
    parser.add_argument('--version', help='Alouette version',
                        action='store_true')
    parser.add_argument('--libs', help='Alouette linker flags',
                        action='store_true')
    parser.add_argument('--cflags', help='Alouette compiler flags',
                        action='store_true')

    args = parser.parse_args(argv)

    if args.prefix:
        print(PREFIX)
    elif args.version:
        print(version.VERSION)
    elif args.libs:
        print(f'-L{PREFIX}/lib -Wl,-rpath,{PREFIX}/lib -lalouette -lm')
    elif args.cflags:
        print(f'-I{PREFIX}/include')
    else:
        parser.print_usage(sys.stderr)
        sys.exit(1)


if __name__ == '__main__': # pragma: no cover
    main()
