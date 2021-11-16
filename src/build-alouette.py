from cffi import FFI
import io
import os
from pathlib import Path
from pcpp.preprocessor import Preprocessor
import re


PREFIX = Path(__file__).parent.parent.resolve()

INCLUDE_PATHS = ('include',)

HEADERS = ('include/alouette.h',)


def format_source(*paths):
    pruned = []
    for path in paths:
        for prefix in INCLUDE_PATHS:
            if path.startswith(prefix):
                pruned.append(path[len(prefix)+1:])
                break
    paths = pruned

    source = os.linesep.join([
        f'#include "{path}"' for path in paths])

    return source


def load_headers(*paths):
    '''
    Load the given header file(s) and run a C preprocessor.
    '''
    headers = []
    for path in paths:
        with open(PREFIX / path) as f:
            header_content = f.read()

        header_content = re.sub(r'(?m)^#include.*\n?', '', header_content)

        cpp = Preprocessor()
        cpp.parse(header_content)
        output = io.StringIO()
        cpp.write(output)

        headers.append(output.getvalue())
    headers = os.linesep.join(headers)

    headers += '''
        extern "Python" void _polarisation_callback(
            int pid,
            const double momentum[3],
            double * polarisation);
'''

    return headers


ffi = FFI()
ffi.set_source('alouette._core', format_source(*HEADERS),
    include_dirs=[f'{PREFIX}/{path}' for path in INCLUDE_PATHS],
    library_dirs=(str(PREFIX / 'alouette'),),
    extra_link_args=(f'-Wl,-rpath,$ORIGIN/.',),
    libraries=('alouette',))
ffi.cdef(load_headers(*HEADERS))


if __name__ == '__main__':
    ffi.compile(verbose=True)
