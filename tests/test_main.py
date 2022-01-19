import os
import pytest

from alouette.__main__ import main
from alouette import version


def test_main(capsys):
    '''Test main script'''

    # Test the no args case
    with pytest.raises(SystemExit) as excp:
        main([])
    assert excp.value.code == 1
    assert capsys.readouterr().err.startswith('usage')

    # Test the prefix case
    main(['--prefix'])
    prefix = capsys.readouterr().out.strip()
    assert prefix[0] == '/'

    # Test the version case
    main(['--version'])
    assert capsys.readouterr().out.strip() == version.VERSION

    # Test the libs case
    main(['--libs'])
    libdir = os.path.join(prefix, 'lib')
    assert libdir in capsys.readouterr().out

    # Test the cflags case
    main(['--cflags'])
    incdir = os.path.join(prefix, 'include')
    assert incdir in capsys.readouterr().out
