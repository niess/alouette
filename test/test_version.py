from alouette import version


def test_vesion():
    '''Test the version sub-module
    '''
    assert(isinstance(version.version, str))
    assert(isinstance(version.git_revision, str))
