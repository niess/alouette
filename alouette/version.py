from ._core import ffi, lib


__all__ = ('GIT_REVISION', 'VERSION')


VERSION = ffi.string(lib.alouette_version()).decode()
GIT_REVISION = ffi.string(lib.alouette_git_revision()).decode()
