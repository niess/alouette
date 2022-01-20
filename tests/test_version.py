from alouette import version


def test_vesion():
    '''Test the version sub-module'''
    assert isinstance(version.VERSION, str)
    assert isinstance(version.GIT_REVISION, str)

    with open('VERSION') as f:
        v = f.read().strip()
    assert(version.VERSION == v)
